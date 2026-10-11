# Staging sequence studies: coast time, crosswind, tip-off, and the Monte Carlo
# of the ignition inhibits.

import argparse
import csv
import math
import random

import jsbsim

ROOT = "."
DT = 0.001
MPH = 1.46667
COASTS = list(range(0, 17))
GATE_FRACTION = 0.70
TILT_MAX = 20.0
WIND_MEAN_FPS, WIND_MAX_FPS = 11.73, 29.0
IMP_SIGMA = 0.025
CD_SLN = math.sqrt(math.log(1 + 0.20 ** 2))


def launch_tilt(rng, nominal_tilt=1.0):
    tn = nominal_tilt + rng.gauss(0, 1.0)
    te = rng.gauss(0, 1.0)
    return 90.0 - math.hypot(tn, te), math.degrees(math.atan2(te, tn)) % 360.0
BOOSTER_TANK, SUST_TANK = 1, 0


def fly(coast, gate_alt=0.0, wind=None, tipoff_dps=0.0, tipoff_dir=0.0,
        k=None, elev=None, use_inhibits=1.0, t_max=150.0, iyy_sep=None, psi=None, k_stack=1.0):
    k = k or {}
    f = jsbsim.FGFDMExec(ROOT)
    f.set_debug_level(0)
    f.load_model("myrocket2")
    f.set_dt(DT)
    f.load_ic("reset", True)
    if elev is not None:
        f["ic/theta-deg"] = elev
    if psi is not None:
        f["ic/psi-true-deg"] = psi
    f.run_ic()
    f["fcs/k-stack"] = k_stack
    f["fcs/target-apogee-ft"] = 1.0e9
    f["fcs/coast-sec"] = coast
    f["fcs/gate-alt-ft"] = gate_alt
    f["fcs/tilt-max-deg"] = TILT_MAX
    f["fcs/use-inhibits"] = use_inhibits
    f["fcs/k-cd"] = k.get("cd", 1.0)
    f["fcs/k-cn"] = k.get("cn", 1.0)
    for eng, kk in ((0, "imp1"), (1, "imp2")):
        f[f"propulsion/engine[{eng}]/isp"] = f[f"propulsion/engine[{eng}]/isp"] * k.get(kk, 1.0)
    if "dry" in k:
        f["inertia/pointmass-weight-lbs[0]"] = 1.8 + (k["dry"] - 1.0) * 37.0
        f["fcs/stage1-weight-lbs"] = f["fcs/stage1-weight-lbs"] * k["dry"]
    if "cp" in k:
        f["metrics/aero-rp-x-in"] = f["metrics/aero-rp-x-in"] * k["cp"]
    if tipoff_dps and iyy_sep:
        m = iyy_sep * math.radians(tipoff_dps) / 0.05
        f["fcs/tipoff-pitch-ftlb"] = m * math.cos(math.radians(tipoff_dir))
        f["fcs/tipoff-yaw-ftlb"] = m * math.sin(math.radians(tipoff_dir))
    if wind:
        spd, frm = wind
        to = math.radians(frm + 180.0)
        f["atmosphere/wind-north-fps"] = spd * math.cos(to)
        f["atmosphere/wind-east-fps"] = spd * math.sin(to)
    f["propulsion/engine[0]/set-running"] = 1
    f["fcs/throttle-cmd-norm"] = 1.0
    apo, t, tilt_max_seen, iyy = 0.0, 0.0, 0.0, None
    while t < t_max:
        f.run()
        t = f["simulation/sim-time-sec"]
        h = f["position/h-agl-ft"]
        apo = max(apo, h)
        if f["fcs/staged"] > 0.5 and iyy is None:
            iyy = f["inertia/iyy-slugs_ft2"]
        if f["fcs/staged"] > 0.5 and f["fcs/decision-now"] < 0.5:
            tilt_max_seen = max(tilt_max_seen, f["fcs/tilt-deg"])
        if t > 5 and f["velocities/h-dot-fps"] < 0:
            break
    out = dict(t_bo1=f["fcs/t-burnout"], apogee=apo, lit=int(f["fcs/sust-lit"] > 0.5), lockout=int(f["fcs/lockout"] > 0.5),
               dec_tilt=f["fcs/lit-tilt-deg"], dec_alt=f["fcs/lit-alt-ft"],
               tilt_max_coast=tilt_max_seen, iyy_sep=iyy)
    out["lock_alt"] = int(out["lockout"] and out["dec_alt"] < gate_alt)
    out["lock_tilt"] = int(out["lockout"] and out["dec_tilt"] > TILT_MAX)
    return out


def nominal_table(elev=None, k_stack=1.0):
    return {c: fly(c, elev=elev, k_stack=k_stack) for c in COASTS}


def rayleigh(rng):
    s = WIND_MEAN_FPS / math.sqrt(math.pi / 2)
    while True:
        v = s * math.sqrt(-2 * math.log(1 - rng.random()))
        if v <= WIND_MAX_FPS:
            return v, rng.uniform(0, 360)


def write(rows, path):
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print("wrote", path)


def main():
    global ROOT
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=["coast", "wind", "tipoff", "mc"])
    ap.add_argument("--n", type=int, default=300)
    ap.add_argument("--coasts", default="0,4,8")
    ap.add_argument("--root", default=ROOT)
    ap.add_argument("--csv", default=None)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--tilt", type=float, default=1.0, help="nominal launch tilt from vertical, deg")
    ap.add_argument("--k-stack", type=float, default=1.0, dest="k_stack")
    a = ap.parse_args()
    ROOT = a.root
    nom = nominal_table(90.0 - a.tilt, a.k_stack)
    gate = {c: GATE_FRACTION * nom[c]["dec_alt"] for c in COASTS}
    iyy = nom[0]["iyy_sep"]
    rows = []
    if a.mode == "coast":
        for c in COASTS:
            r = nom[c]
            rows.append(dict(coast=c, apogee=r["apogee"], dec_alt=r["dec_alt"],
                             dec_tilt=r["dec_tilt"], gate_alt=gate[c]))
            print(rows[-1], flush=True)
    elif a.mode == "wind":
        for mph in (0, 5, 10, 15, 20):
            for c in COASTS[::2]:
                r = fly(c, gate[c], wind=(mph * MPH, 270.0) if mph else None, elev=90.0 - a.tilt)
                rows.append(dict(wind_mph=mph, coast=c, **{k: r[k] for k in
                            ("apogee", "lit", "lockout", "lock_tilt", "lock_alt", "dec_tilt", "dec_alt")}))
                print(rows[-1], flush=True)
    elif a.mode == "tipoff":
        for dps in (0, 10, 20, 40):
            for c in COASTS[::2]:
                r = fly(c, gate[c], tipoff_dps=dps, tipoff_dir=0.0, iyy_sep=iyy, elev=90.0 - a.tilt)
                rows.append(dict(tipoff_dps=dps, coast=c, **{k: r[k] for k in
                            ("apogee", "lit", "lockout", "lock_tilt", "lock_alt", "dec_tilt", "dec_alt")}))
                print(rows[-1], flush=True)
    else:
        coasts = [float(x) for x in a.coasts.split(",")]
        rng = random.Random(a.seed)
        for i in range(a.n):
            k = dict(imp1=rng.gauss(1, IMP_SIGMA), imp2=rng.gauss(1, IMP_SIGMA),
                     cd=math.exp(rng.gauss(-0.5 * CD_SLN ** 2, CD_SLN)), cn=max(0.5, rng.gauss(1, 0.10)),
                     cp=rng.gauss(1, 0.05), dry=1 + rng.uniform(-0.0075, 0.0075))
            elev, psi = launch_tilt(rng, a.tilt)
            w = rayleigh(rng)
            for c in coasts:
                ci = int(round(c))
                r = fly(c, gate[ci], wind=w, k=k, elev=elev, psi=psi, k_stack=a.k_stack)
                rows.append(dict(sample=i, coast=c, wind_fps=w[0], wind_from=w[1], elev=elev, psi=psi, **k,
                                 **{kk: r[kk] for kk in ("apogee", "lit", "lockout", "lock_tilt",
                                                         "lock_alt", "dec_tilt", "dec_alt", "t_bo1")}))
            if (i + 1) % 25 == 0:
                print(f"  {i + 1}/{a.n}", flush=True)
    if a.csv:
        write(rows, a.csv)


if __name__ == "__main__":
    main()
