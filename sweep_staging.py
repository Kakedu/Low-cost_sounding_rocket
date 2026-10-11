# Booster dry-mass sweep to find where staging stops paying off.

import argparse
import csv
import math

import jsbsim

ROOT = "."
DT = 0.001
IGNITION_S = 0.0
NO_BRAKE_TARGET = 1.0e9
SEP_TIME = 0.2

PROP_PER_MOTOR = 6.837
SUSTAINER_DRY = 37.00
PAYLOAD = 1.80
ISP_S = 198.17
G0_FPS2 = 32.174

SERIES_APOGEE_FT = None
COAST_S = 0.0
K_STACK = 1.0
SERIES_DRY = 41.69


def fly_two_stage(stage1_dry_lb, t_max=200.0):
    fdm = jsbsim.FGFDMExec(ROOT)
    fdm.set_debug_level(0)
    fdm.load_model("myrocket2")
    fdm.set_dt(DT)
    fdm.load_ic("reset", True)
    fdm["ic/h-agl-ft"] = 11.5
    fdm.run_ic()

    fdm["fcs/target-apogee-ft"] = NO_BRAKE_TARGET
    fdm["fcs/sep-delay-sec"] = SEP_TIME
    fdm["fcs/coast-sec"] = COAST_S
    fdm["fcs/k-stack"] = K_STACK
    fdm["fcs/stage1-weight-lbs"] = stage1_dry_lb

    lit = False
    apogee = 0.0
    while fdm["simulation/sim-time-sec"] < t_max:
        if not lit and fdm["simulation/sim-time-sec"] >= IGNITION_S:
            fdm["propulsion/engine[0]/set-running"] = 1
            fdm["fcs/throttle-cmd-norm"] = 1.0
            lit = True
        fdm.run()
        h = fdm["position/h-agl-ft"]
        if h > apogee:
            apogee = h
        if fdm["simulation/sim-time-sec"] > 5.0 and h < 1.0:
            break
    return apogee


def ideal_dv_two_stage(stage1_dry_lb):
    mp = PROP_PER_MOTOR
    md1, md2 = stage1_dry_lb, SUSTAINER_DRY

    mi2 = md2 + mp + PAYLOAD
    mi1 = md1 + mp + mi2

    dv1 = ISP_S * G0_FPS2 * math.log(mi1 / (mi1 - mp))
    dv2 = ISP_S * G0_FPS2 * math.log(mi2 / (mi2 - mp))
    return dv1 + dv2


def ideal_dv_series():
    mp = 2.0 * PROP_PER_MOTOR
    mi = SERIES_DRY + mp + PAYLOAD
    return ISP_S * G0_FPS2 * math.log(mi / (mi - mp))


def delta_of(stage1_dry_lb):
    return stage1_dry_lb / (stage1_dry_lb + PROP_PER_MOTOR)


def crossing(xs, ys, level):
    for (x0, y0), (x1, y1) in zip(zip(xs, ys), zip(xs[1:], ys[1:])):
        if (y0 - level) * (y1 - level) <= 0.0 and y0 != y1:
            return x0 + (level - y0) * (x1 - x0) / (y1 - y0)
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--min", type=float, default=2.0, dest="lo")
    ap.add_argument("--max", type=float, default=40.0, dest="hi")
    ap.add_argument("--steps", type=int, default=16)
    ap.add_argument("--coast", type=float, default=0.0)
    ap.add_argument("--k-stack", type=float, default=1.0, dest="k_stack")
    ap.add_argument("--csv", default=None, dest="csv_out")
    a = ap.parse_args()

    dv_series = ideal_dv_series()

    print("=" * 78)
    print("Break-even stage-1 dry mass")
    print("=" * 78)
    global SERIES_APOGEE_FT, COAST_S, K_STACK
    COAST_S = a.coast
    K_STACK = a.k_stack
    from stage_compare import fly as _fly
    SERIES_APOGEE_FT = _fly("myrocket_series", None, closed_loop=False)[0]
    print(f"\nUnstaged control, series burn : {SERIES_APOGEE_FT:9.1f} ft simulated")
    print(f"                                {dv_series:9.1f} fps ideal delta-v")
    print(f"Stage-1 propellant fixed at     {PROP_PER_MOTOR:9.3f} lb, one M1670")
    print()
    print(f"{'stage 1 dry':>12s}{'delta':>8s}{'apogee':>11s}{'vs series':>11s}"
          f"{'ideal dv':>11s}{'vs series':>11s}")
    print(f"{'(lb)':>12s}{'':>8s}{'(ft)':>11s}{'(ft)':>11s}{'(fps)':>11s}{'(fps)':>11s}")
    print("-" * 78)

    masses, apogees, dvs = [], [], []
    for i in range(a.steps):
        m = a.lo + (a.hi - a.lo) * i / max(a.steps - 1, 1)
        try:
            ap_ft = fly_two_stage(m)
        except Exception as e:
            print(f"{m:12.2f}  ERROR {e}")
            continue
        dv = ideal_dv_two_stage(m)
        masses.append(m)
        apogees.append(ap_ft)
        dvs.append(dv)
        print(f"{m:12.2f}{delta_of(m):8.3f}{ap_ft:11.1f}"
              f"{ap_ft - SERIES_APOGEE_FT:+11.1f}{dv:11.1f}{dv - dv_series:+11.1f}")

    print("-" * 78)

    m_sim = crossing(masses, apogees, SERIES_APOGEE_FT)
    m_ideal = crossing(masses, dvs, dv_series)

    print("\nBreak-even")
    print("-" * 78)
    if m_sim is None:
        print("  Simulated: no crossing inside the swept range.")
        print("  Either widen the range, or staging never pays for this vehicle")
        print("  at any physically achievable stage-1 mass - which is itself the")
        print("  answer, and a stronger one than a single point comparison.")
    else:
        print(f"  Simulated (6-DOF, with drag and gravity loss):")
        print(f"      stage 1 dry mass {m_sim:6.2f} lb,  delta = {delta_of(m_sim):.3f}")

    if m_ideal is None:
        print("  Ideal delta-v: no crossing inside the swept range.")
    else:
        print(f"  Ideal delta-v (no drag, no gravity loss):")
        print(f"      stage 1 dry mass {m_ideal:6.2f} lb,  delta = {delta_of(m_ideal):.3f}")

    if m_sim is not None and m_ideal is not None:
        gap = m_sim - m_ideal
        print()
        print(f"  The simulated break-even is {abs(gap):.2f} lb "
              f"{'lighter' if gap < 0 else 'heavier'} than the ideal one.")
        if gap < 0:
            print("  Drag and gravity loss make staging HARDER to justify than the")
            print("  classical result suggests: the booster has to be lighter than")
            print("  ideal theory would allow before splitting the vehicle pays.")
        else:
            print("  Drag and gravity loss make staging EASIER to justify than the")
            print("  classical result suggests.")
        print("  That gap is the point of running a 6-DOF model rather than")
        print("  trusting the textbook form, and it is the result to report.")


    if a.csv_out:
        with open(a.csv_out, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["stage1_dry_lb", "delta", "apogee_ft",
                        "apogee_minus_series_ft", "ideal_dv_fps",
                        "ideal_dv_minus_series_fps"])
            for m, apg, dv in zip(masses, apogees, dvs):
                w.writerow([f"{m:.2f}", f"{delta_of(m):.4f}", f"{apg:.1f}",
                            f"{apg - SERIES_APOGEE_FT:.1f}", f"{dv:.1f}",
                            f"{dv - dv_series:.1f}"])
        print(f"\nwrote {a.csv_out}")




if __name__ == "__main__":
    main()
