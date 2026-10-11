# Launch angle x wind sweep, JSBSim against the 3-DOF model.

import argparse
import csv
import math
from pathlib import Path

import crosscheck as cc
import sweep

MPH = 1.46667
ELEVS = [70, 72, 74, 76, 78, 80, 82, 84, 86, 88, 89.9]
WINDS = [(0, 0, "none")] + [(s, d, n) for s in (10, 20)
                            for d, n in ((0, "head"), (180, "tail"), (270, "cross"))]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True)
    ap.add_argument("--csv", default="anglesweep_v6.csv")
    ap.add_argument("--rail-travel", type=float, default=None, dest="rail",
                    help="guided rail travel, ft, for both models (default: the model's 12.0 ft; "
                         "18.0 ft = a 23 ft rail, rail exit above 100 ft/s)")
    a = ap.parse_args()
    extra = {}
    if a.rail is not None:
        extra["rail_travel"] = a.rail
        cc.RAIL_TRAVEL_FT = a.rail
    for k, v in sweep._reroot(Path(a.root)).items():
        setattr(sweep, k, v)
    rows = []
    for spd, frm, name in WINDS:
        for el in ELEVS:
            w = (spd * MPH, frm) if spd else None
            sweep.build_case({"launch_elev": el, **extra}, open_loop=True, wind=w)
            r = sweep.run_case()
            to = math.radians(frm + 180.0)
            wn = spd * MPH * 0.3048 * math.cos(to)
            we = spd * MPH * 0.3048 * math.sin(to)
            p = cc.fly3(elev_deg=el, wind_n=wn, wind_e=we)
            row = dict(wind_mph=spd, wind=name, elev=el,
                       jsb_apogee_ft=r["apogee"], pm_apogee_ft=p["apogee_ft"],
                       jsb_apo_n=r["apo_n"], jsb_apo_e=r["apo_e"],
                       pm_apo_n=p["apo_n"], pm_apo_e=p["apo_e"],
                       jsb_land_n=r["land_n"], jsb_land_e=r["land_e"], jsb_land=r["land_dist"],
                       pm_land_n=p["land_n"], pm_land_e=p["land_e"], pm_land=p["land_dist"])
            rows.append(row)
            print(f"{name:5s} {spd:2d} mph  {el:5.1f} deg  apogee {r['apogee']:8.0f} / {p['apogee_ft']:8.0f} ft"
                  f"   landing {r['land_dist']:6.0f} / {p['land_dist']:6.0f} m", flush=True)
    with open(a.csv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print("wrote", a.csv)


if __name__ == "__main__":
    main()
