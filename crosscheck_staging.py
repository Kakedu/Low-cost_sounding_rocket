# 3-DOF check of the staged / cluster / series vehicles and the break-even booster mass.

import math
import crosscheck as cc
from crosscheck import (LB2KG, FT2M, G, rho, a_sound, cd_of_mach, thrust_n,
                        ISP_S, PROP_LB, S_FT2, IGNITION_S, BURN_S)

ELEV_DEG = 89.0
PAD = {"two_stage": (10.73, 6.27), "series": (4.35, 12.65), "cluster": (4.35, 12.65)}
DCD_STACK = 0.247

PAYLOAD_LB = 1.8
STRUCT_LB = 32.30
HW_LB = 4.696
SEP_S = 4.0
AREA = S_FT2 * 0.09290304


def fly(config, stage1_dry_lb=9.0, dt=0.001):
    if config == "two_stage":
        m_dry = STRUCT_LB + HW_LB + PAYLOAD_LB + stage1_dry_lb
    else:
        m_dry = STRUCT_LB + 2 * HW_LB + PAYLOAD_LB
    m = (m_dry + 2 * PROP_LB) * LB2KG
    cc.ELEV_DEG = ELEV_DEG
    cc.START_ALT_FT, cc.RAIL_TRAVEL_FT = PAD[config]
    st, t = cc._state0(), 0.0
    t2 = SEP_S if config in ("two_stage", "series") else IGNITION_S
    dropped = False
    while st["vh"] >= 0.0 or t < IGNITION_S + 1.0:
        f1 = thrust_n(t)
        f2 = thrust_n(t - t2 + IGNITION_S)
        if config == "two_stage" and not dropped and t >= SEP_S:
            m -= stage1_dry_lb * LB2KG
            dropped = True
        f = f1 + f2
        m -= f / (ISP_S * G) * dt
        v = math.hypot(st["vx"], st["vh"])
        cd = cd_of_mach(v / a_sound(st["h"]))
        if config == "two_stage" and not dropped:
            cd += DCD_STACK
        cc._step(st, f, cd, m, dt)
        t += dt
    return st["h"] / FT2M


def main():
    jsb = {"series": 20549.7, "cluster": 20197.7, "two_stage": 20017.2}
    print(f"{'case':12s}{'3-DOF (ft)':>12s}{'JSBSim (ft)':>13s}{'diff':>9s}")
    res = {}
    for c in ("series", "cluster", "two_stage"):
        res[c] = fly(c)
        print(f"{c:12s}{res[c]:12.0f}{jsb[c]:13.0f}{100*(res[c]-jsb[c])/jsb[c]:+8.2f}%")
    print(f"\nthrust-profile effect (cluster - series): 3-DOF {res['cluster']-res['series']:+.0f} ft, "
          f"JSBSim {jsb['cluster']-jsb['series']:+.0f} ft")
    print(f"staging effect (two stage - series):      3-DOF {res['two_stage']-res['series']:+.0f} ft, "
          f"JSBSim {jsb['two_stage']-jsb['series']:+.0f} ft")
    lo, hi = 0.5, 37.0
    for _ in range(30):
        mid = 0.5 * (lo + hi)
        if fly("two_stage", mid) > res["series"]:
            lo = mid
        else:
            hi = mid
    print(f"\nbreak-even stage-1 dry mass: 3-DOF {0.5*(lo+hi):.2f} lb")


if __name__ == "__main__":
    main()
