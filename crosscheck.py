# Independent 3-DOF point-mass model, used to check the JSBSim results.
# Inputs are typed in here on purpose, not read from the XML files.

import argparse
import csv
import math
import sys

LB2KG = 0.45359237
FT2M = 0.3048
LBF2N = 4.4482216
PA2PSF = 0.0208854
G = 9.80665

DRY_LB = 38.8
PROP_LB = 6.837
ISP_S = 198.17
BURN_S = 3.9
THRUST_CURVE = [(0.0, 0.0), (0.055, 100), (0.092, 1500), (0.1, 2000),
                (0.15, 2200), (0.2, 1800), (0.5, 1950), (1.0, 2034),
                (1.5, 2000), (2.0, 1900), (2.5, 1760), (2.9, 1700),
                (3.0, 1650), (3.3, 530), (3.4, 350), (3.9, 0.0)]


def thrust_n(t):
    tb = t - IGNITION_S
    if tb <= 0.0 or tb >= BURN_S:
        return 0.0
    for (t0, f0), (t1, f1) in zip(THRUST_CURVE, THRUST_CURVE[1:]):
        if t0 <= tb < t1:
            return f0 + (f1 - f0) * (tb - t0) / (t1 - t0)
    return 0.0
S_FT2 = 0.0845
IGNITION_S = 0.0
START_ALT_FT = 5.0
SITE_FT = 0.0
ELEV_DEG = 89.0
RAIL_TRAVEL_FT = 12.0

CD_TABLE = [(0.00, 0.45), (0.50, 0.45), (0.80, 0.48), (0.95, 0.72),
            (1.05, 0.78), (1.30, 0.70), (2.00, 0.58), (3.00, 0.48)]

BRAKE_DCD = 2.0


def rho(h):
    h = h + SITE_FT * FT2M
    if h < 11000.0:
        t = 288.15 - 0.0065 * h
        p = 101325.0 * (t / 288.15) ** 5.2559
    else:
        t = 216.65
        p = 22632.0 * math.exp(-G * (h - 11000.0) / (287.05 * t))
    return p / (287.05 * t)


def a_sound(h):
    h = h + SITE_FT * FT2M
    t = 288.15 - 0.0065 * h if h < 11000.0 else 216.65
    return math.sqrt(1.4 * 287.05 * t)


def cd_of_mach(m):
    if m <= CD_TABLE[0][0]:
        return CD_TABLE[0][1]
    if m >= CD_TABLE[-1][0]:
        return CD_TABLE[-1][1]
    for (m0, c0), (m1, c1) in zip(CD_TABLE, CD_TABLE[1:]):
        if m0 <= m <= m1:
            return c0 + (c1 - c0) * (m - m0) / (m1 - m0)


def _state0():
    return dict(x=0.0, h=START_ALT_FT * FT2M, vx=0.0, vh=0.0, s=0.0,
                on_rail=True)


def _step(st, f, cd_total, m, dt):
    area = S_FT2 * 0.09290304
    th = math.radians(ELEV_DEG)
    if st["on_rail"]:
        v = math.hypot(st["vx"], st["vh"])
        drag = 0.5 * rho(st["h"]) * v * v * area * cd_total
        a = (f - drag) / m - G * math.sin(th)
        v = max(0.0, v + a * dt)
        st["s"] += v * dt
        st["vx"], st["vh"] = v * math.cos(th), v * math.sin(th)
        st["x"] = st["s"] * math.cos(th)
        st["h"] = START_ALT_FT * FT2M + st["s"] * math.sin(th)
        if st["s"] >= RAIL_TRAVEL_FT * FT2M:
            st["on_rail"] = False
            st.setdefault("v_rail", v)
        return
    v = math.hypot(st["vx"], st["vh"])
    ux, uh = (st["vx"] / v, st["vh"] / v) if v > 0 else (math.cos(th), math.sin(th))
    drag = 0.5 * rho(st["h"]) * v * v * area * cd_total
    st["vx"] += (f - drag) / m * ux * dt
    st["vh"] += ((f - drag) / m * uh - G) * dt
    st["x"] += st["vx"] * dt
    st["h"] += st["vh"] * dt


def fly(brake_window=None, brake_level=0.0, dt=0.001):
    m = (DRY_LB + PROP_LB) * LB2KG
    st = _state0()
    t = 0.0
    burn_end = IGNITION_S + BURN_S
    out = dict(qmax=0.0, t_qmax=0.0, mach_max=0.0)

    while st["vh"] >= 0.0 or t < IGNITION_S + 0.2:
        f = thrust_n(t)
        m -= f / (ISP_S * G) * dt

        b = 0.0
        if brake_window is not None:
            tc = t - burn_end
            if brake_window[0] <= tc < brake_window[1]:
                b = brake_level

        v = math.hypot(st["vx"], st["vh"])
        mach = v / a_sound(st["h"])
        q = 0.5 * rho(st["h"]) * v * v
        if q > out["qmax"]:
            out["qmax"], out["t_qmax"] = q, t
        out["mach_max"] = max(out["mach_max"], mach)

        _step(st, f, cd_of_mach(mach) + BRAKE_DCD * b, m, dt)
        t += dt

        if abs(t - burn_end) < dt:
            out.update(t_burnout=t, h_burnout=st["h"], v_burnout=st["vh"],
                       mach_burnout=mach)

    out.update(t_apogee=t, h_apogee=st["h"], x_apogee=st["x"],
               v_rail=st.get("v_rail", float("nan")))
    return out


def coast_to_apogee(h, vx, vh, m, dt=0.002):
    st = dict(x=0.0, h=h, vx=vx, vh=vh, s=0.0, on_rail=False)
    while st["vh"] > 0.0:
        v = math.hypot(st["vx"], st["vh"])
        _step(st, 0.0, cd_of_mach(v / a_sound(st["h"])), m, dt)
    return st["h"]


DROGUE_CDA_FT2 = 5.8
MAIN_CDA_FT2 = 81.6
MAIN_ALT_FT = 1500.0


def fly3(elev_deg=None, wind_n=0.0, wind_e=0.0, dt=0.002, descent=True):
    th = math.radians(ELEV_DEG if elev_deg is None else elev_deg)
    area = S_FT2 * 0.09290304
    m = (DRY_LB + PROP_LB) * LB2KG
    n, e, h = 0.0, 0.0, START_ALT_FT * FT2M
    vn = ve = vh = 0.0
    s_rail, on_rail, t = 0.0, True, 0.0
    burn_end = IGNITION_S + BURN_S
    apo = None
    cda = 0.0
    while True:
        f = thrust_n(t)
        m -= f / (ISP_S * G) * dt
        if on_rail:
            v = math.hypot(math.hypot(vn, ve), vh)
            drag = 0.5 * rho(h) * v * v * area * cd_of_mach(v / a_sound(h))
            v = max(0.0, v + ((f - drag) / m - G * math.sin(th)) * dt)
            s_rail += v * dt
            vn, ve, vh = v * math.cos(th), 0.0, v * math.sin(th)
            n, h = s_rail * math.cos(th), START_ALT_FT * FT2M + s_rail * math.sin(th)
            if s_rail >= RAIL_TRAVEL_FT * FT2M:
                on_rail = False
            t += dt
            continue
        rn, re_, rh = vn - wind_n, ve - wind_e, vh
        vr = math.sqrt(rn * rn + re_ * re_ + rh * rh)
        if vr > 1e-9:
            un, ue, uh = rn / vr, re_ / vr, rh / vr
        else:
            un, ue, uh = math.cos(th), 0.0, math.sin(th)
        if apo is not None:
            target = DROGUE_CDA_FT2 if h > MAIN_ALT_FT * FT2M else MAIN_CDA_FT2
            cda += (target - cda) * min(1.0, dt / 0.5)
        q = 0.5 * rho(h) * vr * vr
        drag = q * (area * cd_of_mach(vr / a_sound(h)) + cda * 0.09290304)
        thr = f if apo is None else 0.0
        an = (thr - drag) / m * un
        ae = (thr - drag) / m * ue
        ah = (thr - drag) / m * uh - G
        vn += an * dt; ve += ae * dt; vh += ah * dt
        n += vn * dt; e += ve * dt; h += vh * dt
        t += dt
        if apo is None and t > burn_end and vh <= 0.0:
            apo = (h / FT2M, n, e)
            if not descent:
                break
        if apo is not None and h <= 0.0:
            break
    out = dict(apogee_ft=apo[0], apo_n=apo[1], apo_e=apo[2])
    if descent:
        out.update(land_n=n, land_e=e, land_dist=math.hypot(n, e), t_land=t)
    return out


def cmd_verify(csv_out=None):
    r = fly()
    rows = [
        ("Burnout time (s)", r["t_burnout"], "{:.2f}"),
        ("Burnout altitude AGL (ft)", r["h_burnout"] / FT2M, "{:.0f}"),
        ("Burnout Mach number", r["mach_burnout"], "{:.2f}"),
        ("Rail exit speed (fps)", r["v_rail"] / FT2M, "{:.1f}"),
        ("Burnout vertical speed (fps)", r["v_burnout"] / FT2M, "{:.0f}"),
        ("Apogee time (s)", r["t_apogee"], "{:.2f}"),
        ("Apogee altitude AGL (ft)", r["h_apogee"] / FT2M, "{:.0f}"),
        ("Downrange at apogee (m)", r["x_apogee"], "{:.0f}"),
        ("Max dynamic pressure (kPa)", r["qmax"] / 1000.0, "{:.1f}"),
    ]
    print("3-DOF point-mass cross-check, brake stowed. These are the 3-DOF check")
    print("figures in Table 4. Compare against a JSBSim open-loop run:")
    print("  edit the gate altitude in myrocket.xml to a large number, or")
    print("  use  python sweep.py nominal --open\n")
    print(f"{'Quantity':32s}{'3-DOF':>12s}")
    print("-" * 44)
    for name, val, fmt in rows:
        print(f"{name:32s}{fmt.format(val):>12s}")
    if csv_out:
        with open(csv_out, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["quantity", "one_dof"])
            for name, val, fmt in rows:
                w.writerow([name, fmt.format(val)])
        print(f"\nwrote {csv_out}")
        print("Add a jsbsim column from your run, then plot the difference.")
    return rows


def cmd_authority(csv_out=None):
    base = fly()["h_apogee"] / FT2M
    windows = [(0, 2), (2, 4), (4, 6), (6, 8), (8, 10), (10, 15), (15, 20), (20, 22)]
    print(f"Un-braked apogee {base:.0f} ft. Each row holds the brake fully")
    print("deployed for that window only, all else stowed.\n")
    print(f"{'window after burnout (s)':>26s}{'apogee (ft)':>13s}"
          f"{'authority (ft)':>16s}{'per second':>12s}")
    print("-" * 67)
    rows, total = [], 0.0
    for t0, t1 in windows:
        ap = fly(brake_window=(t0, t1), brake_level=1.0)["h_apogee"] / FT2M
        auth = base - ap
        total += auth
        rows.append((f"{t0}-{t1}", ap, auth, auth / (t1 - t0)))
        print(f"{f'{t0}-{t1}':>26s}{ap:13.0f}{auth:16.0f}{auth/(t1-t0):12.1f}")
    full = fly(brake_window=(0, 25), brake_level=1.0)["h_apogee"] / FT2M
    first4 = base - fly(brake_window=(0, 4), brake_level=1.0)["h_apogee"] / FT2M
    print("-" * 67)
    print(f"{'full coast':>26s}{full:13.0f}{base-full:16.0f}")
    print(f"{'first 4 s (one run)':>26s}{base-first4:13.0f}{first4:16.0f}")
    print(f"\nFirst 4 s give {first4:.0f} ft of the {base-full:.0f} ft total, "
          f"{100*first4/(base-full):.0f}%. Note this is computed as a single")
    print("run, not by adding the 0-2 and 2-4 rows: authority is not additive.")
    if csv_out:
        with open(csv_out, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["window_s", "apogee_ft", "authority_ft", "authority_per_second_ft"])
            for r in rows:
                w.writerow([r[0], f"{r[1]:.1f}", f"{r[2]:.1f}", f"{r[3]:.2f}"])
        print(f"wrote {csv_out}")
    return rows


def cmd_estimator(csv_out=None):
    area = S_FT2 * 0.09290304
    dt = 0.001
    st, t, mm = _state0(), 0.0, (DRY_LB + PROP_LB) * LB2KG
    burn_end = IGNITION_S + BURN_S
    samples = []
    while st["vh"] >= 0.0 or t < IGNITION_S + 0.2:
        f = thrust_n(t)
        mm -= f / (ISP_S * G) * dt
        v = math.hypot(st["vx"], st["vh"])
        _step(st, f, cd_of_mach(v / a_sound(st["h"])), mm, dt)
        t += dt
        if t > burn_end and st["vh"] > 0 and round(t * 1000) % 500 == 0:
            samples.append((t - burn_end, st["h"], st["vx"], st["vh"], mm))

    print(f"{'t after burnout':>16s}{'truth':>9s}{'const 0.75':>12s}{'err':>8s}"
          f"{'drag form':>11s}{'err':>8s}   (ft)")
    print("-" * 66)
    e1s, e2s, rows = [], [], []
    for tc, hh, vxx, vv, mmm in samples:
        truth = coast_to_apogee(hh, vxx, vv, mmm) / FT2M
        e1 = (hh + 0.75 * vv * vv / (2 * G)) / FT2M
        vt2 = 2 * mmm * G / (rho(hh) * area * 0.45)
        e2 = (hh + (vt2 / (2 * G)) * math.log(1 + vv * vv / vt2)) / FT2M
        e1s.append(abs(e1 - truth)); e2s.append(abs(e2 - truth))
        rows.append((tc, truth, e1, e1 - truth, e2, e2 - truth))
        if tc <= 8:
            print(f"{tc:16.1f}{truth:9.0f}{e1:12.0f}{e1-truth:+8.0f}"
                  f"{e2:11.0f}{e2-truth:+8.0f}")
    print("-" * 66)
    print(f"mean |error| over the coast:  const {sum(e1s)/len(e1s):6.0f} ft"
          f"   drag form {sum(e2s)/len(e2s):6.0f} ft")
    if csv_out:
        with open(csv_out, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["t_after_burnout_s", "truth_ft", "const_ft", "const_err_ft",
                        "dragform_ft", "dragform_err_ft"])
            for r in rows:
                w.writerow([f"{r[0]:.2f}"] + [f"{x:.1f}" for x in r[1:]])
        print(f"wrote {csv_out}")
    return rows


def _closed_loop(estimator, kp=0.001, target_ft=10000.0, k_cd=1.0, k_T=1.0,
                 k_m=1.0, rate=2.0, dt=0.001):
    area = S_FT2 * 0.09290304
    m = (DRY_LB * k_m + PROP_LB) * LB2KG
    burn_end = IGNITION_S + BURN_S
    st, t, brake = _state0(), 0.0, 0.0
    while st["vh"] >= 0.0 or t < IGNITION_S + 0.2:
        powered = IGNITION_S <= t < burn_end
        f = thrust_n(t) * k_T
        m -= thrust_n(t) / (ISP_S * G) * dt
        h, v = st["h"], st["vh"]
        cmd = 0.0
        if not powered and t > burn_end and v > 0 and h > 1500 * FT2M:
            if estimator == "const":
                est = (h + 0.75 * v * v / (2 * G)) / FT2M
            else:
                vt2 = 2 * m * G / (rho(h) * area * 0.45)
                est = (h + (vt2 / (2 * G)) * math.log(1 + v * v / vt2)) / FT2M
            cmd = max(0.0, min(1.0, kp * (est - target_ft)))
        brake += max(-rate * dt, min(rate * dt, cmd - brake))
        vv = math.hypot(st["vx"], st["vh"])
        mach = vv / a_sound(h)
        _step(st, f, cd_of_mach(mach) * k_cd + BRAKE_DCD * brake, m, dt)
        t += dt
    return st["h"] / FT2M


def cmd_control(csv_out=None):
    cases = [("nominal", 1.0, 1.0, 1.0), ("Cd +10%", 1.1, 1.0, 1.0),
             ("Cd -10%", 0.9, 1.0, 1.0), ("thrust +5%", 1.0, 1.05, 1.0),
             ("thrust -5%", 1.0, 0.95, 1.0), ("dry mass +5%", 1.0, 1.0, 1.05),
             ("dry mass -5%", 1.0, 1.0, 0.95)]
    print("Closed-loop apogee error, target 10000 ft, in feet.")
    print("The estimator always assumes the nominal drag coefficient.\n")
    print(f"{'case':16s}{'const 0.75':>13s}{'drag form':>12s}")
    print("-" * 41)
    rows, e1, e2 = [], [], []
    for name, kc, kt, km in cases:
        a1 = _closed_loop("const", k_cd=kc, k_T=kt, k_m=km) - 10000.0
        a2 = _closed_loop("drag", k_cd=kc, k_T=kt, k_m=km) - 10000.0
        e1.append(abs(a1)); e2.append(abs(a2))
        rows.append((name, a1, a2))
        print(f"{name:16s}{a1:+13.0f}{a2:+12.0f}")
    print("-" * 41)
    print(f"{'worst case':16s}{max(e1):13.0f}{max(e2):12.0f}")
    if csv_out:
        with open(csv_out, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["case", "const_error_ft", "dragform_error_ft"])
            for r in rows:
                w.writerow([r[0], f"{r[1]:.1f}", f"{r[2]:.1f}"])
        print(f"wrote {csv_out}")
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=["verify", "authority", "estimator", "control", "all"])
    ap.add_argument("--csv", default=None, dest="csv_out")
    a = ap.parse_args()
    if a.mode in ("verify", "all"):
        cmd_verify(a.csv_out if a.mode == "verify" else None); print()
    if a.mode in ("authority", "all"):
        cmd_authority(a.csv_out if a.mode == "authority" else None); print()
    if a.mode in ("estimator", "all"):
        cmd_estimator(a.csv_out if a.mode == "estimator" else None)
        print()
    if a.mode in ("control", "all"):
        cmd_control(a.csv_out if a.mode == "control" else None)


if __name__ == "__main__":
    main()
