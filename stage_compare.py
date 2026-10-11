# Same two motors flown three ways: cluster, series, two-stage.

import argparse
import csv

import jsbsim

ROOT = "."
DT = 0.001
IGNITION_S = 0.0
NO_BRAKE_TARGET = 1.0e9

IMPULSE_PER_MOTOR_NS = 6042.0


def fly(aircraft, sep_time=None, closed_loop=True, ignite_all=False, t_max=200.0):
    fdm = jsbsim.FGFDMExec(ROOT)
    fdm.set_debug_level(0)
    fdm.load_model(aircraft)
    fdm.set_dt(DT)
    fdm.load_ic("reset", True)
    fdm.run_ic()

    if not closed_loop:
        fdm["fcs/target-apogee-ft"] = NO_BRAKE_TARGET
    if sep_time is not None:
        try:
            fdm["fcs/sep-delay-sec"] = sep_time
        except Exception:
            pass

    lit = False
    apogee = 0.0
    maxq = 0.0
    vmax = 0.0
    t_last_thrust = 0.0
    while fdm["simulation/sim-time-sec"] < t_max:
        if not lit and fdm["simulation/sim-time-sec"] >= IGNITION_S:
            fdm["propulsion/engine[0]/set-running"] = 1
            fdm["fcs/throttle-cmd-norm"] = 1.0
            if ignite_all:
                fdm["propulsion/engine[1]/set-running"] = 1
                fdm["fcs/throttle-cmd-norm[1]"] = 1.0
            lit = True
        fdm.run()
        h = fdm["position/h-agl-ft"]
        if h > apogee:
            apogee = h
        q = fdm["aero/qbar-psf"]
        if q > maxq:
            maxq = q
        v = fdm["velocities/vt-fps"]
        if v > vmax:
            vmax = v
        thrust = fdm["propulsion/engine[0]/thrust-lbs"]
        try:
            thrust += fdm["propulsion/engine[1]/thrust-lbs"]
        except Exception:
            pass
        if thrust > 1.0:
            t_last_thrust = fdm["simulation/sim-time-sec"]
        if fdm["simulation/sim-time-sec"] > 5.0 and h < 1.0:
            break
    return apogee, maxq, vmax, t_last_thrust


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", default=None, dest="csv_out")
    a = ap.parse_args()

    cases = [
        ("single stage, open loop",      "myrocket",  None,  False, 1),
        ("single stage, closed loop",    "myrocket",  None,  True,  1),
        ("two stage, open loop",         "myrocket2", 0.2,   False, 2),
        ("two stage, closed loop",       "myrocket2", 0.2,   True,  2),
        ("series, open loop",            "myrocket_series",  None, False, 2),
        ("series, closed loop",          "myrocket_series",  None, True,  2),
        ("cluster, open loop",           "myrocket_cluster", None, False, 2),
        ("cluster, closed loop",         "myrocket_cluster", None, True,  2),
        ("two stage, staging inhibited", "myrocket2", 999.0, False, 1),
    ]

    print("=" * 78)
    print("Does staging pay for this vehicle?")
    print("=" * 78)
    print(f"\n{'configuration':30s}{'impulse':>9s}{'apogee':>10s}{'max q':>9s}"
          f"{'v max':>9s}{'burn end':>9s}")
    print(f"{'':30s}{'(N.s)':>9s}{'(ft)':>10s}{'(psf)':>9s}{'(fps)':>9s}{'(s)':>9s}")
    print("-" * 76)

    rows = []
    for label, ac, sep, closed, motors in cases:
        try:
            apogee, maxq, vmax, t_thrust = fly(
                ac, sep_time=sep, closed_loop=closed,
                ignite_all=ac.endswith("cluster"))
        except Exception as e:
            print(f"{label:32s}  ERROR: {e}")
            continue
        impulse = IMPULSE_PER_MOTOR_NS * motors
        eff = apogee / (impulse / 1000.0)
        rows.append((label, impulse, apogee, eff, maxq, vmax, t_thrust))
        print(f"{label:30s}{impulse:9.0f}{apogee:10.1f}{maxq:9.1f}"
              f"{vmax:9.1f}{t_thrust:9.2f}")

    print("-" * 78)

    def get(label):
        for r in rows:
            if r[0] == label:
                return r
        return None

    s_open = get("single stage, open loop")
    t_open = get("two stage, open loop")
    s_closed = get("single stage, closed loop")
    t_closed = get("two stage, closed loop")

    print("\nReading the table")
    print("-" * 76)

    c_open = get("cluster, open loop")
    c_closed = get("cluster, closed loop")
    r_open = get("series, open loop")
    r_closed = get("series, closed loop")

    print("DECOMPOSITION - three vehicles at 12084 N.s")
    print("  cluster   two motors together, one airframe, nothing dropped")
    print("  series    two motors in sequence, one airframe, nothing dropped")
    print("  two stage two motors in sequence, two airframes, one dropped")
    print()

    if c_open and r_open:
        d = c_open[2] - r_open[2]
        print(f"  thrust profile effect, cluster minus series : {d:+9.1f} ft")
    if r_open and t_open:
        d = t_open[2] - r_open[2]
        print(f"  PURE STAGING EFFECT,   two stage minus series: {d:+9.1f} ft"
              f"   ({100.0 * d / r_open[2]:+.1f}%)")
        print()
        if d > 0:
            print("  Staging pays. Dropping the spent airframe returns more than")
            print("  the second structure costs at this booster mass.")
        else:
            print("  Staging does not pay for this vehicle. Carrying a second")
            print("  airframe to burnout costs more than dropping it returns.")
            print("  Whether it pays depends on the booster dry mass:")
            print("  see sweep_staging.py for the break-even mass.")

    if c_open and r_open and t_open:
        print()
        print("  Check the mechanism against the columns rather than assuming it.")
        print("  If the cluster wins mainly on thrust profile, it should show a")
        print("  much shorter burn and a higher max q: high thrust-to-weight")
        print("  cuts the time spent climbing against gravity, at the cost of")
        print("  more drag. Gravity loss scales with time under thrust, so a")
        print(f"  burn ending at {c_open[6]:.1f} s against {r_open[6]:.1f} s is")
        print("  worth roughly g times that difference in lost velocity.")

    if s_closed and t_closed:
        print()
        print("CONTROLLER PERFORMANCE - target 10000 ft")
        for label in ("single stage, closed loop", "series, closed loop",
                      "cluster, closed loop", "two stage, closed loop"):
            r = get(label)
            if r is None:
                continue
            print(f"  {label:28s} apogee {r[2]:9.1f} ft, error {r[2] - 10000.0:+9.1f} ft")
        print("  The controller was tuned on the single-stage vehicle. Errors")
        print("  grow with the energy the brake has to remove, which is the")
        print("  gain-scheduling result again: loop gain has to track the")
        print("  vehicle, not just the brake.")

    if a.csv_out:
        with open(a.csv_out, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["configuration", "impulse_Ns", "apogee_ft", "ft_per_kNs"])
            for r in rows:
                w.writerow([r[0], f"{r[1]:.0f}", f"{r[2]:.1f}", f"{r[3]:.1f}"])
        print(f"\nwrote {a.csv_out}")




if __name__ == "__main__":
    main()
