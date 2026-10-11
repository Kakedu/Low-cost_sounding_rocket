# Main driver: writes a modified copy of the model, runs JSBSim, reads the log.
# Main modes: nominal, oat, recovery, dmap, mcw (Monte Carlo), morris. -h for the rest.

import argparse
import csv
import math
import random
import re
import shutil
import subprocess
import sys
from pathlib import Path

JSBSIM_ROOT = Path(__file__).resolve().parent
CASE_NAME = "_sweep"

EXE = JSBSIM_ROOT / "JSBSim.exe"
SRC_AC = JSBSIM_ROOT / "aircraft" / "myrocket"
SRC_ENG = JSBSIM_ROOT / "engine" / "M1670.xml"
SRC_SCRIPT = JSBSIM_ROOT / "scripts" / "launch.xml"

CASE_AC = JSBSIM_ROOT / "aircraft" / CASE_NAME
CASE_ENG = JSBSIM_ROOT / "engine" / f"{CASE_NAME}_engine.xml"
CASE_SCRIPT = JSBSIM_ROOT / "scripts" / f"{CASE_NAME}.xml"
CASE_CSV = JSBSIM_ROOT / "data_output" / f"{CASE_NAME}.csv"


def _scale_number(text, pattern, factor, group=1):
    m = re.search(pattern, text, re.S)
    if not m:
        raise KeyError(f"pattern not found: {pattern[:60]}")
    val = float(m.group(group))
    new = val * factor
    s, e = m.span(group)
    return text[:s] + f"{new:.6g}" + text[e:]


def _scale_table_column(text, fn_name, factor):
    m = re.search(rf'<function name="{re.escape(fn_name)}">.*?</function>', text, re.S)
    if not m:
        raise KeyError(f"function not found: {fn_name}")
    block = m.group(0)
    td = re.search(r'<tableData>(.*?)</tableData>', block, re.S)
    rows = []
    for line in td.group(1).splitlines():
        parts = line.split()
        if len(parts) == 2:
            rows.append(f"              {float(parts[0]):8.4f}   {float(parts[1])*factor:9.5f}")
        elif line.strip():
            rows.append(line)
    new_block = block[:td.start(1)] + "\n" + "\n".join(rows) + "\n            " + block[td.end(1):]
    return text[:m.start()] + new_block + text[m.end():]


def _scale_engine_table(text, col, factor):
    td = re.search(r'<tableData>(.*?)</tableData>', text, re.S)
    if not td:
        raise KeyError("thrust table not found in engine file")
    out = []
    for line in td.group(1).splitlines():
        m = re.match(r'(\s*)([-\d.]+)\s+([-\d.]+)(.*)$', line)
        if m:
            v = [float(m.group(2)), float(m.group(3))]
            v[col] *= factor
            out.append(f"{m.group(1)}{v[0]:9.4f}  {v[1]:9.2f}{m.group(4)}")
        else:
            out.append(line)
    return text[:td.start(1)] + "\n".join(out) + text[td.end(1):]


PARAMS = {
    "cd0": dict(
        file="ac", kind="table", fn="aero/force/Drag_mach",
        desc="zero-incidence drag coefficient (whole Mach table)"),
    "cn_alpha": dict(
        file="ac", kind="tables",
        fns=["aero/force/Lift_alpha", "aero/force/Side_beta"],
        desc="normal-force slope CNalpha (lift and side scaled together)"),
    "cp_x": dict(
        file="ac", kind="scale",
        pat=r'<location name="AERORP"[^>]*>\s*<x>\s*([\d.]+)\s*</x>',
        desc="centre of pressure station (drives static margin)"),
    "pitch_damp": dict(
        file="ac", kind="scale",
        pat=r'aero/moment/Pitch_damp.*?<value>\s*(-?[\d.]+)\s*</value>',
        desc="pitch damping derivative"),
    "yaw_damp": dict(
        file="ac", kind="scale",
        pat=r'aero/moment/Yaw_damp.*?<value>\s*(-?[\d.]+)\s*</value>',
        desc="yaw damping derivative (the axis a crosswind actually excites)"),
    "roll_damp": dict(
        file="ac", kind="scale",
        pat=r'aero/moment/Roll_damp.*?<value>\s*(-?[\d.]+)\s*</value>',
        desc="roll damping derivative"),
    "airbrake_dcd": dict(
        file="ac", kind="scale",
        pat=r'aero/force/Drag_airbrake.*?<value>\s*([\d.]+)\s*</value>',
        desc="air-brake drag increment at full deployment"),

    "dry_mass": dict(
        file="ac", kind="scale",
        pat=r'<emptywt unit="LBS"\s*>\s*([\d.]+)\s*</emptywt>',
        desc="dry mass"),
    "payload": dict(
        file="ac", kind="scale",
        pat=r'<pointmass name="Payload">.*?<weight unit="LBS">\s*([\d.]+)\s*</weight>',
        desc="payload mass"),
    "iyy": dict(
        file="ac", kind="scale",
        pat=r'<iyy unit="SLUG\*FT2">\s*([\d.]+)\s*</iyy>',
        desc="pitch moment of inertia"),
    "cg_shift": dict(
        file="ac", kind="scale",
        pat=r'<location name="CG" unit="IN">\s*<x>\s*([\d.]+)\s*</x>',
        desc="dry centre of gravity station"),

    "prop_mass": dict(
        file="ac", kind="prop",
        pats=[r'<tank type="FUEL".*?<capacity unit="LBS">\s*([\d.]+)',
              r'<tank type="FUEL".*?<contents unit="LBS">\s*([\d.]+)'],
        desc="propellant mass = total impulse (curve shape and thrust level "
             "kept, burn time scales)"),
    # impulse is scaled through Isp: JSBSim won't fill a tank above its capacity,
    # so scaling the propellant only works downwards
    "isp": dict(
        file="eng", kind="scale",
        pat=r'<isp>\s*([\d.]+)\s*</isp>',
        desc="specific impulse (thrust curve kept, propellant burns slower: "
             "impulse and burn time scale)"),
    "mdot": dict(
        file="eng", kind="eng_thrust",
        desc="thrust level at fixed isp (flow scales with it, burn shortens, "
             "impulse unchanged)"),

    "launch_elev": dict(
        file="ic", kind="absolute",
        pat=r'<theta unit="DEG">\s*([\d.]+)\s*</theta>',
        oat=(88.0, 89.9),
        desc="launch elevation angle, deg (88 / 89.9, not +/-10%)"),

    "launch_az": dict(
        file="ic", kind="absolute",
        pat=r'<psi   unit="DEG">\s*([\d.]+)\s*</psi>',
        oat=(0.0, 0.0),
        desc="launch azimuth, deg (Monte Carlo only)"),

    "rail_travel": dict(
        file="ac", kind="absolute",
        pat=r'<property value="([\d.]+)">fcs/rail-travel-ft</property>',
        oat=(12.0, 12.0),
        desc="guided rail travel, ft (absolute; anglesweep.py rail variant, not OAT)"),

    "gust": dict(
        file="script", kind="gust",
        oat=(0.0, 40.0),
        desc="discrete gust at t=2.5 s, fps (0 / 40, not +/-10%)"),

    "wind": dict(
        file="script", kind="wind",
        oat=(0.0, 30.0),
        desc="crosswind, fps (0 / 30, not +/-10%)"),

    "drogue_cda": dict(
        file="ac", kind="scale",
        pat=r'<property value="([\d.]+)">fcs/drogue-cda-sqft</property>',
        desc="drogue parachute CdA (assumed sizing, ~75 fps descent)"),
    "main_cda": dict(
        file="ac", kind="scale",
        pat=r'<property value="([\d.]+)">fcs/main-cda-sqft</property>',
        desc="main parachute CdA (assumed sizing, ~20 fps descent)"),

    "act_time": dict(
        file="ac", kind="scale",
        pat=r'<kinematic name="fcs/airbrake-actuator">.*?'
            r'<position>\s*1\.0\s*</position>\s*<time>\s*([\d.]+)\s*</time>',
        desc="actuator travel time, stowed to fully deployed"),
    "gate_alt": dict(
        file="ac", kind="scale",
        pat=r'position/h-agl-ft gt ([\d.]+)',
        desc="minimum altitude before the brake is allowed to act"),

    "kp": dict(
        file="ac", kind="scale",
        pat=r'<pure_gain name="fcs/brake-cmd-raw">.*?<gain>\s*([\d.]+)\s*</gain>',
        desc="proportional gain"),
    "target": dict(
        file="ac", kind="absolute",
        pat=r'<property value="([\d.]+)">fcs/target-apogee-ft</property>',
        oat=(10000.0, 10000.0),
        desc="target apogee, ft (absolute; used by the margin curve, not OAT)"),
    "cd_est": dict(
        file="ac", kind="scale",
        pat=r'<fcs_function name="fcs/coast-vt2">.*?'
            r'<property>metrics/Sw-sqft</property>\s*<value>\s*([\d.]+)\s*</value>',
        desc="Cd the on-board estimator assumes (a wrong belief, not a real change)"),
}


def load_templates():
    t = {}
    t["ac"] = (SRC_AC / "myrocket.xml").read_text()
    t["ic"] = (SRC_AC / "reset.xml").read_text()
    t["eng"] = SRC_ENG.read_text()
    t["script"] = SRC_SCRIPT.read_text()
    return t


def apply_param(text, name, factor):
    spec = PARAMS[name]
    kind = spec["kind"]
    if kind == "table":
        return _scale_table_column(text, spec["fn"], factor)
    if kind == "scale":
        return _scale_number(text, spec["pat"], factor)
    if kind == "tables":
        for fn in spec["fns"]:
            text = _scale_table_column(text, fn, factor)
        return text
    if kind == "eng_thrust":
        return _scale_engine_table(text, 1, factor)
    if kind in ("scale_all", "prop"):
        for p in spec["pats"]:
            text = _scale_number(text, p, factor)
        return text
    if kind == "gust":
        ev = (f'\n    <event name="Gust">\n'
              f'      <condition> simulation/sim-time-sec ge 2.5 </condition>\n'
              f'      <set name="atmosphere/wind-east-fps" value="{factor:.4g}"/>\n'
              f'    </event>\n')
        return text.replace('  </run>', ev + '\n  </run>', 1)
    if kind == "wind":
        anchor = '<set name="fcs/throttle-cmd-norm" value="1.0"/>'
        assert anchor in text, "ignition event not found in script"
        return text.replace(
            anchor,
            anchor + f'\n      <set name="atmosphere/wind-east-fps" value="{factor:.4g}"/>',
            1)
    if kind == "absolute":
        m = re.search(spec["pat"], text, re.S)
        s, e = m.span(1)
        return text[:s] + f"{factor:.6g}" + text[e:]
    raise ValueError(kind)


def check_template(open_loop):
    t = load_templates()["ac"]
    m = re.search(r'position/h-agl-ft gt ([\d.]+)', t)
    if not m:
        return
    gate = float(m.group(1))
    if not open_loop and gate > 100000:
        print("!" * 68)
        print(f"WARNING: gate altitude in the template is {gate:.0f} ft.")
        print("The brake can never act, so this 'closed loop' study is open")
        print("loop. Restore 1500.0 in the aircraft file before trusting any")
        print("result from this run.")
        print("!" * 68)
        print()
    elif open_loop and gate <= 100000:
        print("note: open-loop runs set an unreachable target internally.\n")


def build_case(factors, open_loop=False, wind=None):
    t = load_templates()
    factors = dict(factors)
    dT = factors.pop("dT_K", 0.0)
    if dT:
        anchor = '<set name="fcs/throttle-cmd-norm" value="1.0"/>'
        assert anchor in t["script"], "ignition event not found in script"
        # delta-T is in Rankine, hence the 1.8
        t["script"] = t["script"].replace(
            anchor, anchor + f'\n      <set name="atmosphere/delta-T" value="{1.8 * dT:.4f}"/>', 1)
    if wind is not None:
        spd, frm = wind
        to = math.radians(frm + 180.0)
        vn, ve = spd * math.cos(to), spd * math.sin(to)
        anchor = '<set name="fcs/throttle-cmd-norm" value="1.0"/>'
        assert anchor in t["script"], "ignition event not found in script"
        t["script"] = t["script"].replace(
            anchor,
            anchor + f'\n      <set name="atmosphere/wind-north-fps" value="{vn:.4f}"/>'
                     f'\n      <set name="atmosphere/wind-east-fps" value="{ve:.4f}"/>', 1)
    for name, f in factors.items():
        spec = PARAMS[name]
        t[spec["file"]] = apply_param(t[spec["file"]], name, f)
        if spec["kind"] == "prop":
            t["eng"] = _scale_engine_table(t["eng"], 0, f)

    assert t["ac"].count('<engine file="M1670">') == 1, "engine reference not found exactly once"
    t["ac"] = t["ac"].replace('<engine file="M1670">',
                              f'<engine file="{CASE_NAME}_engine">')
    t["eng"] = re.sub(r'<rocket_engine name="[^"]*">',
                      f'<rocket_engine name="{CASE_NAME}_engine">', t["eng"])

    t["ac"] = re.sub(r'<output name="[^"]*"',
                     f'<output name="data_output/{CASE_NAME}.csv"', t["ac"])

    if open_loop:
        t["ac"], n_sub = re.subn(r'<property value="[\d.]+">fcs/target-apogee-ft</property>',
                                 '<property value="1000000000.0">fcs/target-apogee-ft</property>', t["ac"])
        assert n_sub == 1, "target property not found exactly once"

    t["script"] = t["script"].replace('aircraft="myrocket"', f'aircraft="{CASE_NAME}"')

    CASE_AC.mkdir(parents=True, exist_ok=True)
    (CASE_AC / f"{CASE_NAME}.xml").write_text(t["ac"], newline="\r\n")
    (CASE_AC / "reset.xml").write_text(t["ic"], newline="\r\n")
    CASE_ENG.write_text(t["eng"], newline="\r\n")
    CASE_SCRIPT.write_text(t["script"], newline="\r\n")


_BINDING_RUNNER = """
import os, sys, jsbsim
f = jsbsim.FGFDMExec(os.getcwd())
f.set_debug_level(0)
f.load_script(sys.argv[1])
f.run_ic()
while f.run():
    pass
"""


def run_case():
    if CASE_CSV.exists():
        CASE_CSV.unlink()
    if EXE.exists():
        cmd = [str(EXE), f"--script=scripts/{CASE_NAME}.xml"]
    else:
        cmd = [sys.executable, "-c", _BINDING_RUNNER, f"scripts/{CASE_NAME}.xml"]
    r = subprocess.run(cmd, cwd=str(JSBSIM_ROOT), capture_output=True, text=True)
    if not CASE_CSV.exists():
        tail = "\n".join((r.stdout + r.stderr).splitlines()[-15:])
        raise RuntimeError("no CSV produced:\n" + tail)
    return read_results(CASE_CSV)


def read_results(path):
    with open(path) as f:
        rows = list(csv.reader(f))
    hdr = [h.strip() for h in rows[0]]

    def idx(*frags):
        for fr in frags:
            for i, h in enumerate(hdr):
                if fr.lower() in h.lower():
                    return i
        raise KeyError(frags)

    i_t, i_h = idx("time"), idx("altitude agl")
    i_m, i_q = idx("mach"), idx("dynamic pressure")
    i_b = idx("airbrake position")
    i_a, i_be = idx("alpha"), idx("beta")
    i_th = idx("pitch angle")
    try:
        i_d = idx("ground track")
    except KeyError:
        i_d = None
    data = [[float(x) for x in r] for r in rows[1:] if len(r) == len(hdr)]
    apo_row = max(data, key=lambda r: r[i_h])
    asc = [r for r in data if r[i_t] <= apo_row[i_t]]
    flight = [r for r in asc if r[i_t] > 1.0] or asc

    thrust_max = max(r[idx("thrust")] for r in data)
    if apo_row[i_t] < 1.0 or thrust_max < 1.0:
        raise RuntimeError(
            f"vehicle did not fly: peak altitude {apo_row[i_h]:.0f} ft at "
            f"t={apo_row[i_t]:.2f} s, peak thrust {thrust_max:.1f} lbf, "
            f"{len(data)} rows. Peak thrust near zero means the engine "
            f"produced nothing - check the engine file.")
    pos = {}
    try:
        i_la, i_lo = idx("latitude"), idx("longitude")
        R = 6378137.0
        la0, lo0 = data[0][i_la], data[0][i_lo]
        def _ne(r):
            return (math.radians(r[i_la] - la0) * R,
                    math.radians(r[i_lo] - lo0) * R * math.cos(math.radians(la0)))
        pos["apo_n"], pos["apo_e"] = _ne(apo_row)
        pos["land_n"], pos["land_e"] = _ne(data[-1])
        pos["land_dist"] = math.hypot(pos["land_n"], pos["land_e"])
        pos["t_land"] = data[-1][i_t]
    except KeyError:
        pass
    return dict(**pos, apogee=apo_row[i_h], t_apogee=apo_row[i_t],
                max_mach=max(r[i_m] for r in asc),
                max_q=max(r[i_q] for r in asc),
                max_brake=max(r[i_b] for r in data),
                max_alpha=max(abs(r[i_a]) for r in flight),
                max_beta=max(abs(r[i_be]) for r in flight),
                theta_apogee=apo_row[i_th],
                drift=(max(r[i_d] for r in data) if i_d is not None else float("nan")))


def cmd_check():
    t = load_templates()
    ok = bad = 0
    for name, spec in PARAMS.items():
        try:
            before = t[spec["file"]]
            after = apply_param(before, name, 1.1)
            if after == before:
                print(f"  NO-OP   {name:14s} {spec['desc']}")
                bad += 1
            else:
                print(f"  ok      {name:14s} {spec['desc']}")
                ok += 1
        except Exception as e:
            print(f"  FAIL    {name:14s} {type(e).__name__}: {e}")
            bad += 1
    print(f"\n{ok} ok, {bad} need attention")
    return bad == 0


def cmd_nominal(open_loop=False, csv_out=None):
    check_template(open_loop)
    build_case({}, open_loop)
    r = run_case()
    for k, v in r.items():
        print(f"  {k:12s} {v:10.2f}")
    if csv_out:
        import shutil
        shutil.copy(CASE_CSV, csv_out)
        print(f"wrote {csv_out}")
    return r


def _per_sigma(name, lo, hi, delta):
    sd = {**MC_SIGMA, **{k: a / math.sqrt(3.0) for k, a in MC_UNIFORM.items()}}.get(name)
    if sd is None or lo != lo or hi != hi:
        return None
    return (hi - lo) / (2.0 * delta) * sd


def cmd_oat(delta=0.10, open_loop=False, metric="apogee", csv_out=None):
    check_template(open_loop)
    build_case({}, open_loop)
    base = run_case()[metric]
    unit = {"apogee": "ft", "drift": "m", "max_q": "psf",
            "max_alpha": "deg", "max_beta": "deg"}.get(metric, "")
    print(f"baseline {metric} {base:.1f} {unit}   "
          f"(perturbation +/-{delta:.0%}, {'open' if open_loop else 'closed'} loop)\n")
    out = []
    for name, spec in PARAMS.items():
        pts = spec.get("oat", (1 - delta, 1 + delta))
        vals = {}
        for sgn, f in (("lo", pts[0]), ("hi", pts[1])):
            try:
                build_case({name: f}, open_loop)
                vals[sgn] = run_case()[metric] - base
            except Exception as e:
                vals[sgn] = float("nan")
                print(f"  ! {name} {sgn}: {e}")
        rng = abs(vals["hi"] - vals["lo"])
        out.append((name, vals["lo"], vals["hi"], rng, spec["desc"]))
    out.sort(key=lambda r: -r[3])
    if csv_out:
        with open(csv_out, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["parameter", "low", "high", "range", "pct_of_base",
                        "per_sigma", "baseline", "metric", "loop", "description"])
            for name, lo, hi, rng, desc in out:
                ps = _per_sigma(name, lo, hi, delta)
                w.writerow([name, f"{lo:.3f}", f"{hi:.3f}", f"{rng:.3f}",
                            f"{100*rng/abs(base):.4f}", "" if ps is None else f"{ps:.3f}",
                            f"{base:.3f}", metric, "open" if open_loop else "closed", desc])
        print(f"\nwrote {csv_out}")
    print(f"{'parameter':16s}{'low':>10s}{'high':>10s}{'range':>10s}"
          f"{'% of base':>11s}{'per 1 sd':>10s}  description")
    print("-" * 110)
    for name, lo, hi, rng, desc in out:
        ps = _per_sigma(name, lo, hi, delta)
        pss = "       -" if ps is None else f"{ps:+10.1f}"
        print(f"{name:16s}{lo:+10.1f}{hi:+10.1f}{rng:10.1f}"
              f"{100*rng/abs(base):10.2f}%{pss}  {desc}")
    return out


def cmd_mc(n=500, open_loop=False):
    check_template(open_loop)
    IMPULSE_SIGMA_RPY = 0.0249
    IMPULSE_SIGMA_CERT = 0.067
    DIST = {"prop_mass": IMPULSE_SIGMA_RPY,
            "cd0": 0.03,
            "dry_mass": 0.03,
            "launch_elev": None}
    build_case({}, open_loop)
    base = run_case()["apogee"]
    res = []
    for i in range(n):
        f = {}
        for k, s in DIST.items():
            if k == "launch_elev":
                f[k] = min(89.9, random.gauss(89.0, 1.0))
            else:
                f[k] = max(0.5, random.gauss(1.0, s))
        try:
            build_case(f, open_loop)
            res.append(run_case()["apogee"])
        except Exception:
            pass
        if (i + 1) % 25 == 0:
            print(f"  {i+1}/{n}", end="\r", flush=True)
    res.sort()
    mean = sum(res) / len(res)
    sd = math.sqrt(sum((x - mean) ** 2 for x in res) / (len(res) - 1))
    print(f"\nn = {len(res)}   baseline {base:.0f} ft")
    print(f"  mean   {mean:8.0f} ft")
    print(f"  sd     {sd:8.0f} ft  ({100*sd/mean:.2f}%)")
    for p in (5, 50, 95):
        print(f"  p{p:<2d}    {res[int(p/100*(len(res)-1))]:8.0f} ft")
    print(f"  range  {res[0]:.0f} to {res[-1]:.0f} ft")
    print(f"\n  Note the asymmetry: the brake absorbs overshoot but cannot")
    print(f"  correct undershoot, so the upper tail is clipped and the lower")
    print(f"  tail is not. Quoting a standard deviation alone implies a normal")
    print(f"  distribution and understates the downside.")
    print(f"  p50-p5 = {res[int(0.5*(len(res)-1))]-res[int(0.05*(len(res)-1))]:.0f} ft"
          f"   p95-p50 = {res[int(0.95*(len(res)-1))]-res[int(0.5*(len(res)-1))]:.0f} ft")
    return res


CP_NOM_IN = 52.55
CG_LAUNCH_IN = 44.39
DIAM_IN = 3.94
MARGIN_NOM = (CP_NOM_IN - CG_LAUNCH_IN) / DIAM_IN


def _cp_factor(margin_cal):
    return (CG_LAUNCH_IN + margin_cal * DIAM_IN) / CP_NOM_IN


ATTITUDE_PARAMS = ["cn_alpha", "yaw_damp", "iyy", "cg_shift", "pitch_damp"]


def cmd_stab(gust=40.0, open_loop=False):
    check_template(open_loop)
    """Do the attitude parameters matter?

    In the first sweep every attitude-related input returned zero sensitivity.
    Two explanations fit: they genuinely do not matter for this vehicle, or
    the model contains no disturbance for them to act against. This separates
    the two by lowering the static margin toward neutral and adding a gust,
    then re-running the same parameters.
    """
    margins = [1.2, 2.0, 3.0, MARGIN_NOM, 5.0]
    print(f"nominal static margin {MARGIN_NOM:.2f} calibres, "
          f"gust {gust:.0f} fps at t=2.5 s, "
          f"{'open loop' if open_loop else 'closed loop'}\n")

    print(f"{'margin':>8s}{'CP (in)':>9s}{'apogee':>9s}{'drift(m)':>10s}"
          f"{'max alpha':>11s}{'max beta':>10s}{'theta@apo':>11s}")
    print("-" * 68)
    bases = {}
    for mg in margins:
        f = {"cp_x": _cp_factor(mg), "gust": gust}
        build_case(f, open_loop)
        try:
            r = run_case()
        except Exception as e:
            print(f"{mg:8.2f}   run failed: {e}")
            continue
        bases[mg] = r
        print(f"{mg:8.2f}{_cp_factor(mg)*CP_NOM_IN:9.1f}{r['apogee']:9.0f}"
              f"{r['drift']:10.0f}{r['max_alpha']:11.2f}{r['max_beta']:10.2f}"
              f"{r['theta_apogee']:11.1f}")
        if r['max_beta'] > 15.0:
            print(f"{'':8s}  ^ sideslip exceeds 15 deg: the linear normal-force "
                  f"slope and the Barrowman method it came from are outside "
                  f"their validity range here")

    print(f"\nSensitivity of attitude parameters (+/-10%), apogee change in ft")
    print(f"{'margin':>8s}" + "".join(f"{p:>13s}" for p in ATTITUDE_PARAMS))
    print("-" * (8 + 13 * len(ATTITUDE_PARAMS)))
    for mg in margins:
        if mg not in bases:
            continue
        row = []
        for p in ATTITUDE_PARAMS:
            d = []
            for fac in (0.9, 1.1):
                f = {"cp_x": _cp_factor(mg), "gust": gust, p: fac}
                if p == "cg_shift":
                    f["cp_x"] = _cp_factor(mg) * fac
                build_case(f, open_loop)
                try:
                    d.append(run_case()["apogee"] - bases[mg]["apogee"])
                except Exception:
                    d.append(float("nan"))
            row.append(max(d) - min(d))
        print(f"{mg:8.2f}" + "".join(f"{v:13.1f}" for v in row))

    print(f"\nSensitivity of the same parameters, DRIFT change in m")
    print(f"{'margin':>8s}" + "".join(f"{p:>13s}" for p in ATTITUDE_PARAMS))
    print("-" * (8 + 13 * len(ATTITUDE_PARAMS)))
    for mg in margins:
        if mg not in bases:
            continue
        row = []
        for p in ATTITUDE_PARAMS:
            d = []
            for fac in (0.9, 1.1):
                f = {"cp_x": _cp_factor(mg), "gust": gust, p: fac}
                if p == "cg_shift":
                    f["cp_x"] = _cp_factor(mg) * fac
                build_case(f, open_loop)
                try:
                    d.append(run_case()["drift"] - bases[mg]["drift"])
                except Exception:
                    d.append(float("nan"))
            row.append(max(d) - min(d))
        print(f"{mg:8.2f}" + "".join(f"{v:13.1f}" for v in row))

    print("\nApogee and drift are different objectives. A parameter that is")
    print("negligible for one can dominate the other, and a study that reports")
    print("only apogee cannot say whether a coefficient is worth measuring.")
    print("\nIf the right-hand columns stay near zero at every margin, the")
    print("attitude parameters really are irrelevant for this vehicle and the")
    print("aerodynamic stability coefficients need not be measured accurately.")
    print("If they grow as the margin falls, they matter and the first sweep")
    print("was blind because the vehicle was too stable to be disturbed.")


DOWNSIDE = {
    "isp_low":       {"isp": 0.90},
    "prop_low":      {"prop_mass": 0.90},
    "dry_high":      {"dry_mass": 1.10},
    "cd_high":       {"cd0": 1.10},
    "wind":          {"wind": 30.0},
    "combined":      {"isp": 0.95, "prop_mass": 0.95, "dry_mass": 1.05},
}


def cmd_margin(target=10000.0, impulses=None, open_loop=False, csv_out=None):
    check_template(open_loop)
    if impulses is None:
        impulses = [1.00, 1.05, 1.10, 1.15, 1.20, 1.30]

    print(f"target apogee {target:.0f} ft\n")
    print("Impulse is scaled through propellant mass. 'nominal' is the apogee")
    print("on a perfect day; 'brake' is how far the brake had to open to hold")
    print("the target, which is the price paid for the margin.\n")

    hdr = f"{'impulse':>9s}{'unbraked':>10s}{'nominal':>9s}{'brake':>7s}"
    for k in DOWNSIDE:
        hdr += f"{k:>10s}"
    hdr += f"{'worst':>9s}{'score':>9s}"
    print(hdr)
    print("-" * len(hdr))

    rows = []
    for imp in impulses:
        base_f = {"prop_mass": imp}

        build_case(base_f, open_loop=True)
        try:
            unbraked = run_case()["apogee"]
        except Exception as e:
            print(f"{imp:9.2f}   unbraked run failed: {e}")
            continue

        build_case(base_f, open_loop)
        try:
            r = run_case()
        except Exception as e:
            print(f"{imp:9.2f}   nominal run failed: {e}")
            continue
        nominal, brake = r["apogee"], r["max_brake"]

        errs = []
        cells = ""
        for name, pert in DOWNSIDE.items():
            f = dict(base_f)
            for k, v in pert.items():
                if PARAMS[k]["kind"] in ("absolute", "wind", "gust"):
                    f[k] = v
                else:
                    f[k] = f.get(k, 1.0) * v
            build_case(f, open_loop)
            try:
                e = run_case()["apogee"] - target
            except Exception:
                e = float("nan")
            errs.append(e)
            cells += f"{e:+10.0f}"

        worst = max(abs(e) for e in errs if e == e)
        score = max(worst, abs(nominal - target))
        print(f"{imp:9.2f}{unbraked:10.0f}{nominal:9.0f}{brake:7.2f}"
              f"{cells}{worst:9.0f}{score:9.0f}")
        rows.append((imp, unbraked, nominal, brake, errs, worst, score))

    if csv_out:
        with open(csv_out, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["impulse", "unbraked", "nominal", "max_brake"]
                       + list(DOWNSIDE) + ["worst", "score"])
            for r in rows:
                w.writerow([r[0], f"{r[1]:.1f}", f"{r[2]:.1f}", f"{r[3]:.3f}"]
                           + [f"{e:.1f}" for e in r[4]] + [f"{r[5]:.1f}", f"{r[6]:.1f}"])
        print(f"\nwrote {csv_out}")

    if rows:
        print("\nscore = max(worst downside error, nominal-day error). The")
        print("nominal day is just another condition the vehicle has to fly,")
        print("so it is scored alongside the perturbed cases rather than used")
        print("as a pass/fail filter.")
        best = min(rows, key=lambda r: (r[6], r[0]))
        print(f"\nBest margin: impulse x{best[0]:.2f}  ->  un-braked apogee "
              f"{best[1]:.0f} ft, {best[1]-target:+.0f} ft above target "
              f"({100*(best[1]-target)/target:.0f}% of target), "
              f"score {best[6]:.0f} ft ({100*best[6]/target:.2f}% of target).")
        print("A larger margin keeps buying downward authority until the brake")
        print("saturates on a nominal day; past that the extra impulse is spent")
        print("fighting itself and the nominal case starts to drift off target.")
    return rows


def cmd_grid(target=10000.0, impulses=None, brakes=None, open_loop=False,
             fixed_loop_gain=False, csv_out=None):
    check_template(open_loop)
    """Design chart: how much margin, with how large a brake.

    Margin buys downward authority; brake size sets how fast that authority
    can be spent. Neither alone fixes the sizing problem, so the useful
    output is the pair.
    """
    if impulses is None:
        impulses = [1.00, 1.05, 1.10, 1.15, 1.20, 1.25]
    if brakes is None:
        brakes = [0.5, 1.0, 2.0, 3.0, 4.0]

    if fixed_loop_gain:
        print("kp scaled as 2.0/dCd so that the loop gain kp*dCd is held at")
        print("its nominal value. Compare against the fixed-kp grid.\n")

    print(f"target {target:.0f} ft.\n")
    print("Cell = score = max(worst downside error, nominal-day error), in ft.")
    print("Low is good. The nominal day is scored as one more condition, not")
    print("used as a pass/fail filter, so this table and 'margin' agree.\n")
    print(f"{'impulse':>9s}" + "".join(f"{'dCd '+str(b):>10s}" for b in brakes))
    print("-" * (9 + 10 * len(brakes)))

    table = {}
    for imp in impulses:
        cells = ""
        for b in brakes:
            base_f = {"prop_mass": imp, "airbrake_dcd": b / 2.0}
            if fixed_loop_gain:
                base_f["kp"] = 2.0 / b
            build_case(base_f, open_loop)
            try:
                nominal = run_case()["apogee"]
            except Exception:
                cells += f"{'err':>10s}"
                continue
            errs = []
            for pert in DOWNSIDE.values():
                f = dict(base_f)
                for k, v in pert.items():
                    if PARAMS[k]["kind"] in ("absolute", "wind", "gust"):
                        f[k] = v
                    else:
                        f[k] = f.get(k, 1.0) * v
                build_case(f, open_loop)
                try:
                    errs.append(run_case()["apogee"] - target)
                except Exception:
                    pass
            w = max(abs(e) for e in errs) if errs else float("nan")
            score = max(w, abs(nominal - target))
            table[(imp, b)] = (score, w, nominal - target)
            cells += f"{score:10.0f}"
        print(f"{imp:9.2f}{cells}")

    for label, j in (("worst downside error, ft", 1), ("nominal-day error, ft", 2)):
        print(f"\n{'impulse':>9s}" + "".join(f"{'dCd '+str(b):>10s}" for b in brakes)
              + f"   ({label})")
        for imp in impulses:
            row = ""
            for b in brakes:
                v = table.get((imp, b))
                row += f"{v[j]:+10.0f}" if v else f"{'err':>10s}"
            print(f"{imp:9.2f}{row}")

    if csv_out:
        with open(csv_out, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["impulse", "dCd", "score", "worst_downside", "nominal_error",
                        "target", "gain_schedule"])
            for (imp, b), v in sorted(table.items()):
                if v:
                    w.writerow([imp, b, f"{v[0]:.3f}", f"{v[1]:.3f}", f"{v[2]:.3f}",
                                target, "1/dCd" if fixed_loop_gain else "fixed"])
        print(f"\nwrote {csv_out}")

    good = {k: v[0] for k, v in table.items() if v and v[0] == v[0]}
    if good:
        print(f"\nBest margin for each brake size:")
        print(f"{'brake':>10s}{'impulse':>10s}{'score':>9s}")
        for b in brakes:
            col = {i: table[(i, b)][0] for i in impulses
                   if table.get((i, b)) and table[(i, b)][0] == table[(i, b)][0]}
            if col:
                i = min(col, key=lambda x: (col[x], x))
                print(f"{b:10.1f}{i:10.2f}{col[i]:9.0f}")
        k = min(good, key=lambda x: (good[x], x[1], x[0]))
        print(f"\nLowest score overall: impulse x{k[0]:.2f} with brake dCd "
              f"{k[1]:.1f}, score {good[k]:.0f} ft. Ties are broken toward the"
              f" smaller brake and then the smaller motor, since the score"
              f" ignores cost.")
        print("Read along a row to see what a bigger brake buys at fixed")
        print("margin, and down a column to see what margin buys at fixed")
        print("brake. If a row stops improving, the brake is the constraint;")
        print("if a column does, the margin is.")
    return table


def cmd_residual(target=10000.0, impulse=1.15, dcd=2.0, open_loop=False):
    check_template(open_loop)
    """What sets the residual error at the recommended design point?

    At impulse x1.15 with dCd 2.0 the score is set by the NOMINAL-day error,
    not by the downside cases, and the nominal error is positive everywhere in
    the design grid with a floor near +18 ft even at large authority. That is
    a systematic overshoot on an undisturbed flight, so the cause is in the
    loop itself rather than in disturbance rejection. Four candidates:

      kp         proportional control leaves a steady-state offset by
                 construction: as the error shrinks the command shrinks with
                 it, so the brake retracts before the error reaches zero.
      act_time   0.5 s of travel against a coast where most of the authority
                 is spent in the first four seconds.
      cd_est     the estimator assumes a stowed Cd of 0.45; while the brake is
                 out the true Cd is higher, so the estimate is optimistic.
      gate_alt   the brake is inhibited below 1500 ft AGL, which costs part of
                 the high-dynamic-pressure window just after burnout.
    """
    base = {"prop_mass": impulse, "airbrake_dcd": dcd / 2.0}
    build_case(base, open_loop)
    r0 = run_case()
    e0 = r0["apogee"] - target
    print(f"design point: impulse x{impulse:.2f}, dCd {dcd:.1f}")
    print(f"nominal-day error {e0:+.0f} ft, max deployment {r0['max_brake']:.2f}\n")

    trials = [
        ("kp",       [0.25, 0.5, 1.0, 2.0, 4.0, 10.0]),
        ("act_time", [0.1, 0.2, 0.5, 1.0, 2.0]),
        ("cd_est",   [0.9, 1.0, 1.2, 1.5, 2.0]),
        ("gate_alt", [0.05, 0.33, 1.0, 2.0]),
    ]
    for name, factors in trials:
        nom = {"kp": 0.001, "act_time": 0.5,
               "cd_est": 0.45, "gate_alt": 1500.0}[name]
        print(f"{name}   (nominal {nom:g})")
        print(f"{'value':>12s}{'nom err':>10s}{'max dep':>10s}{'apogee':>10s}")
        for f in factors:
            g = dict(base)
            g[name] = f
            build_case(g, open_loop)
            try:
                r = run_case()
                print(f"{nom*f:12g}{r['apogee']-target:+10.0f}"
                      f"{r['max_brake']:10.2f}{r['apogee']:10.0f}")
            except Exception as ex:
                print(f"{nom*f:12g}   {ex}")
        print()

    print("A term that flattens the nominal error when relaxed is the binding")
    print("one. A term that changes nothing is not the constraint, whatever")
    print("its physical plausibility.")


MC_SIGMA = {"prop_mass": 0.025, "cd0": 0.20, "cn_alpha": 0.10, "cp_x": 0.05,
            "drogue_cda": 0.10, "main_cda": 0.10}
MC_UNIFORM = {"dry_mass": 0.0075}
WIND_MEAN_FPS = 11.73
WIND_MAX_FPS = 29.0
BRAKE_DCD_SIGMA = 0.20
DT_SIGMA_K = 10.0


def _draw_wind(rng):
    sigma = WIND_MEAN_FPS / math.sqrt(math.pi / 2.0)
    while True:
        spd = sigma * math.sqrt(-2.0 * math.log(1.0 - rng.random()))
        if spd <= WIND_MAX_FPS:
            return spd, rng.uniform(0.0, 360.0)


def launch_tilt(rng, nominal_tilt=1.0):
    tn = nominal_tilt + rng.gauss(0.0, 1.0)
    te = rng.gauss(0.0, 1.0)
    return 90.0 - math.hypot(tn, te), math.degrees(math.atan2(te, tn)) % 360.0


def draw_inputs(rng, sig=None):
    sig = sig or MC_SIGMA
    f = {}
    for k, s in sig.items():
        if k == "cd0":
            sl = math.sqrt(math.log(1.0 + s * s))
            f[k] = math.exp(rng.gauss(-0.5 * sl * sl, sl))
        else:
            f[k] = max(0.5, rng.gauss(1.0, s))
    f.update({k: 1.0 + rng.uniform(-a, a) for k, a in MC_UNIFORM.items()})
    f["launch_elev"], f["launch_az"] = launch_tilt(rng)
    w = _draw_wind(rng)
    f["airbrake_dcd"] = max(0.5, rng.gauss(1.0, BRAKE_DCD_SIGMA))
    f["dT_K"] = rng.gauss(0.0, DT_SIGMA_K)
    return f, w


def cmd_mcw(n=500, open_loop=False, csv_out=None, seed=1, oracle=False,
            cd_sigma=None, target=None, imp_sigma=None):
    rng = random.Random(seed)
    check_template(open_loop)
    sig = dict(MC_SIGMA)
    if cd_sigma is not None:
        sig["cd0"] = cd_sigma
    if imp_sigma is not None:
        sig["prop_mass"] = imp_sigma
    rows = []
    keys = None
    for i in range(n):
        f, w = draw_inputs(rng, sig)
        if oracle:
            f["cd_est"] = f["cd0"]
        if target is not None:
            f["target"] = target
        try:
            build_case(f, open_loop, wind=w)
            r = run_case()
        except Exception as e:
            print(f"  ! sample {i}: {e}")
            continue
        row = dict(sample=i, wind_fps=w[0], wind_from_deg=w[1], **f, **r)
        rows.append(row)
        keys = keys or list(row.keys())
        if (i + 1) % 25 == 0:
            print(f"  {i+1}/{n}", flush=True)
    if csv_out:
        with open(csv_out, "w", newline="") as fh:
            wr = csv.DictWriter(fh, fieldnames=keys)
            wr.writeheader()
            wr.writerows(rows)
        print(f"wrote {csv_out}")
    _mc_summary(rows)
    return rows


def _mc_summary(rows):
    import statistics as st
    ap = sorted(r["apogee"] for r in rows)
    q = lambda p: ap[int(p / 100 * (len(ap) - 1))]
    print(f"\nn = {len(rows)}")
    print(f"apogee  mean {st.mean(ap):.0f}  sd {st.stdev(ap):.0f}  "
          f"p5 {q(5):.0f}  p50 {q(50):.0f}  p95 {q(95):.0f} ft")
    N = [r["land_n"] for r in rows]; E = [r["land_e"] for r in rows]
    mn, me = st.mean(N), st.mean(E)
    snn = st.variance(N); see = st.variance(E)
    sne = sum((a - mn) * (b - me) for a, b in zip(N, E)) / (len(N) - 1)
    tr, det = snn + see, snn * see - sne ** 2
    l1 = tr / 2 + math.sqrt(max(0.0, tr * tr / 4 - det))
    l2 = tr / 2 - math.sqrt(max(0.0, tr * tr / 4 - det))
    ang = math.degrees(0.5 * math.atan2(2 * sne, snn - see))
    d = sorted(r["land_dist"] for r in rows)
    print(f"landing centre  N {mn:.0f} m  E {me:.0f} m")
    print(f"1-sigma ellipse semi-axes {math.sqrt(l1):.0f} x {math.sqrt(l2):.0f} m, "
          f"major axis {ang:.0f} deg from north")
    print(f"distance from pad  p50 {d[len(d)//2]:.0f}  p95 {d[int(0.95*(len(d)-1))]:.0f}  "
          f"max {d[-1]:.0f} m")


MORRIS_INPUTS = [
    ("prop_mass", "normal", 0.025), ("cd0", "lognormal", 0.20), ("cn_alpha", "normal", 0.10),
    ("cp_x", "normal", 0.05), ("dry_mass", "uniform", 0.0075),
    ("tilt_n", "tilt", 1.0), ("tilt_e", "tilt", 0.0),
    ("wind_spd", "rayleigh", None), ("wind_dir", "dir", None),
    ("airbrake_dcd", "normal", 0.20), ("dT_K", "dT", 10.0),
]
Q_LO, Q_HI = 0.025, 0.975


def _norm_ppf(q):
    from statistics import NormalDist
    return NormalDist().inv_cdf(q)


def _rayleigh_ppf_capped(q):
    s = WIND_MEAN_FPS / math.sqrt(math.pi / 2.0)
    fmax = 1.0 - math.exp(-WIND_MAX_FPS ** 2 / (2 * s * s))
    return s * math.sqrt(-2.0 * math.log(1.0 - q * fmax))


def _morris_map(u):
    f, spd, frm, tn, te = {}, 0.0, 0.0, 1.0, 0.0
    for (name, kind, s), x in zip(MORRIS_INPUTS, u):
        q = Q_LO + (Q_HI - Q_LO) * x
        if kind == "normal":
            f[name] = 1.0 + s * _norm_ppf(q)
        elif kind == "lognormal":
            sl = math.sqrt(math.log(1.0 + s * s))
            f[name] = math.exp(-0.5 * sl * sl + sl * _norm_ppf(q))
        elif kind == "uniform":
            f[name] = 1.0 + (2 * x - 1) * s
        elif kind == "dT":
            f[name] = s * _norm_ppf(q)
        elif kind == "tilt":
            v = s + _norm_ppf(q)
            if name == "tilt_n":
                tn = v
            else:
                te = v
        elif kind == "rayleigh":
            spd = _rayleigh_ppf_capped(q)
        elif kind == "dir":
            frm = 360.0 * x
    f["launch_elev"] = 90.0 - math.hypot(tn, te)
    f["launch_az"] = math.degrees(math.atan2(te, tn)) % 360.0
    return f, ((spd, frm) if spd > 0 else None)


def cmd_morris(r=20, open_loop=False, csv_out=None, seed=1, p=4, n_boot=1000):
    rng = random.Random(seed)
    k = len(MORRIS_INPUTS)
    delta = p / (2.0 * (p - 1))
    levels = [i / (p - 1) for i in range(p)]
    check_template(open_loop)
    ee = {name: [] for name, _, _ in MORRIS_INPUTS}
    rows = []

    def y(u):
        f, w = _morris_map(u)
        build_case(f, open_loop, wind=w)
        return run_case()["apogee"]

    for t in range(r):
        base = [rng.choice([l for l in levels if l + delta <= 1.0 + 1e-9]) for _ in range(k)]
        order = list(range(k)); rng.shuffle(order)
        u = list(base); y0 = y(u)
        for i in order:
            u2 = list(u); u2[i] = u[i] + delta
            y1 = y(u2)
            e = (y1 - y0) / delta
            ee[MORRIS_INPUTS[i][0]].append(e)
            rows.append(dict(traj=t, input=MORRIS_INPUTS[i][0], ee=e, y0=y0, y1=y1))
            u, y0 = u2, y1
        print(f"  trajectory {t + 1}/{r}", flush=True)
    import statistics as st
    brng = random.Random(seed + 1000)
    print(f"{'input':12s}{'mu*':>9s}{'95% CI':>17s}{'mu':>9s}{'sigma':>9s}"
          "   (ft per unit of the 2.5-97.5 % quantile range)")
    summ = []
    for name, _, _ in MORRIS_INPUTS:
        e = ee[name]
        mus = st.mean(abs(x) for x in e); mu = st.mean(e); sd = st.stdev(e) if len(e) > 1 else 0
        bs = sorted(st.mean(abs(brng.choice(e)) for _ in e) for _ in range(n_boot))
        lo, hi = bs[int(0.025 * n_boot)], bs[int(0.975 * n_boot) - 1]
        summ.append(dict(input=name, mu_star=mus, ci_lo=lo, ci_hi=hi, mu=mu, sigma=sd))
    for d in sorted(summ, key=lambda z: -z["mu_star"]):
        print(f"{d['input']:12s}{d['mu_star']:9.0f}   [{d['ci_lo']:5.0f}, {d['ci_hi']:5.0f}]"
              f"{d['mu']:9.0f}{d['sigma']:9.0f}")
    if csv_out:
        with open(csv_out, "w", newline="") as fh:
            wr = csv.DictWriter(fh, fieldnames=list(rows[0]))
            wr.writeheader(); wr.writerows(rows)
        sp = csv_out.replace(".csv", "_summary.csv")
        with open(sp, "w", newline="") as fh:
            wr = csv.DictWriter(fh, fieldnames=list(summ[0]))
            wr.writeheader(); wr.writerows(summ)
        print(f"wrote {csv_out}, {sp}")


def cmd_recovery(target=10000.0, csv_out=None):
    check_template(False)
    facs = [0.85, 0.90, 0.925, 0.95, 0.975, 1.00, 1.05, 1.10, 1.15, 1.20, 1.25, 1.30]
    print(f"{'impulse':>8s}{'open (ft)':>11s}{'closed (ft)':>13s}{'error':>9s}"
          f"{'brake peak':>12s}")
    out = []
    for f in facs:
        build_case({"prop_mass": f}, open_loop=True)
        op = run_case()["apogee"]
        build_case({"prop_mass": f}, open_loop=False)
        r = run_case()
        out.append((f, op, r["apogee"], r["apogee"] - target, r["max_brake"]))
        print(f"{f:8.3f}{op:11.0f}{r['apogee']:13.0f}{r['apogee']-target:+9.0f}"
              f"{r['max_brake']:12.2f}", flush=True)
    if csv_out:
        with open(csv_out, "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["impulse_factor", "open_apogee_ft", "closed_apogee_ft",
                        "error_ft", "brake_peak"])
            for r in out:
                w.writerow([f"{r[0]:.3f}"] + [f"{x:.1f}" for x in r[1:4]] + [f"{r[4]:.3f}"])
        print(f"wrote {csv_out}")
    return out


MAP_CASES = {
    "nominal":    {},
    "impulse+10": {"prop_mass": 1.10},
    "impulse+20": {"prop_mass": 1.20},
    "cd-10":      {"cd0": 0.90},
    "dry-10":     {"dry_mass": 0.90},
    "cd+10":      {"cd0": 1.10},
}


def cmd_dmap(target=10000.0, csv_out=None):
    check_template(False)
    dcds = [0.5, 1.0, 2.0, 3.0, 4.0]
    kps = [0.25, 0.5, 1.0, 2.0, 4.0]
    print("cell = worst |apogee error| (ft) over " + ", ".join(MAP_CASES))
    print(f"{'dCd / kp x':>10s}" + "".join(f"{k:>8.2f}" for k in kps))
    rows = []
    for d in dcds:
        line = f"{d:10.1f}"
        for k in kps:
            errs = {}
            for cname, c in MAP_CASES.items():
                f = dict(c, airbrake_dcd=d / 2.0, kp=k)
                build_case(f, open_loop=False)
                errs[cname] = run_case()["apogee"] - target
            worst = max(abs(e) for e in errs.values())
            rows.append(dict(dcd=d, kp_factor=k, worst=worst, **errs))
            line += f"{worst:8.0f}"
        print(line, flush=True)
    if csv_out:
        with open(csv_out, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
            w.writeheader(); w.writerows(rows)
        print(f"wrote {csv_out}")
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=["check", "nominal", "oat", "mc", "mcw", "morris", "recovery", "dmap", "stab", "margin", "grid", "residual"])
    ap.add_argument("n", nargs="?", type=int, default=500)
    ap.add_argument("--open", action="store_true", dest="open_loop",
                    help="disable the air brake (open loop)")
    ap.add_argument("--fixed-loop-gain", action="store_true",
                    dest="fixed_loop_gain",
                    help="grid only: scale kp as 1/dCd so kp*dCd is constant")
    ap.add_argument("--metric", default="apogee",
                    choices=["apogee", "drift", "max_q", "max_alpha", "max_beta"],
                    help="which outcome to rank inputs against")
    ap.add_argument("--csv", default=None, dest="csv_out",
                    help="also write the study result to this CSV, for plotting")
    ap.add_argument("--root", default=None)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--oracle", action="store_true")
    ap.add_argument("--cd-sigma", type=float, default=None, dest="cd_sigma")
    ap.add_argument("--target", type=float, default=None)
    ap.add_argument("--imp-sigma", type=float, default=None, dest="imp_sigma",
                    help="mcw: total-impulse 1-sigma other than the literature 0.025")
    a = ap.parse_args()

    if a.root:
        globals().update(_reroot(Path(a.root)))

    if a.mode == "check":
        sys.exit(0 if cmd_check() else 1)
    if a.mode == "nominal":
        cmd_nominal(a.open_loop, csv_out=a.csv_out)
    elif a.mode == "oat":
        cmd_oat(open_loop=a.open_loop, metric=a.metric, csv_out=a.csv_out)
    elif a.mode == "mc":
        cmd_mc(a.n, a.open_loop)
    elif a.mode == "recovery":
        cmd_recovery(csv_out=a.csv_out)
    elif a.mode == "dmap":
        cmd_dmap(csv_out=a.csv_out)
    elif a.mode == "mcw":
        cmd_mcw(a.n, a.open_loop, csv_out=a.csv_out, seed=a.seed, oracle=a.oracle,
                cd_sigma=a.cd_sigma, target=a.target, imp_sigma=a.imp_sigma)
    elif a.mode == "morris":
        cmd_morris(a.n, a.open_loop, csv_out=a.csv_out, seed=a.seed)
    elif a.mode == "stab":
        cmd_stab(open_loop=a.open_loop)
    elif a.mode == "margin":
        cmd_margin(open_loop=a.open_loop, csv_out=a.csv_out)
    elif a.mode == "grid":
        cmd_grid(open_loop=a.open_loop, fixed_loop_gain=a.fixed_loop_gain, csv_out=a.csv_out)
    elif a.mode == "residual":
        cmd_residual(open_loop=a.open_loop)


def _reroot(root):
    return dict(
        JSBSIM_ROOT=root, EXE=root / "JSBSim.exe",
        SRC_AC=root / "aircraft" / "myrocket",
        SRC_ENG=root / "engine" / "M1670.xml",
        SRC_SCRIPT=root / "scripts" / "launch.xml",
        CASE_AC=root / "aircraft" / CASE_NAME,
        CASE_ENG=root / "engine" / f"{CASE_NAME}_engine.xml",
        CASE_SCRIPT=root / "scripts" / f"{CASE_NAME}.xml",
        CASE_CSV=root / "data_output" / f"{CASE_NAME}.csv")


if __name__ == "__main__":
    main()
