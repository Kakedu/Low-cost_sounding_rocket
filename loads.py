# Air-brake and fin loads for the cases the structure has to survive.

import argparse
import csv
import json
import math
from pathlib import Path

import barrowman as B
import sweep

_BJ = Path(__file__).resolve().parent / "barrowman_v6.json"
if _BJ.exists():
    B.FIN["S"] = json.load(open(_BJ))["geometry"]["fin"]["S"]

LBF = 4.4482216
PSF = 47.880259
FT2 = 0.09290304
N_FLAPS = 4


def fin_share():
    c0, _, fac = B.fin_set(B.FIN, B.L1 - B.FIN["Cr"])
    return c0, fac


def run(name, factors, wind=None, force_brake=False):
    sweep.build_case(factors, open_loop=force_brake, wind=wind)
    if force_brake:
        p = Path(sweep.CASE_AC) / f"{sweep.CASE_NAME}.xml"
        s = p.read_text()
        s = s.replace('<input>fcs/brake-cmd-src</input>',
                      '<input>fcs/force-brake</input>')
        s = s.replace('<property value="0.0">fcs/sil-mode</property>',
                      '<property value="0.0">fcs/sil-mode</property>\n  <property value="0.0">fcs/force-brake</property>')
        s = s.replace('   <kinematic name="fcs/airbrake-actuator">',
                      '''   <switch name="fcs/force-brake">
     <default value="0.0"/>
     <test logic="AND" value="1.0">
       propulsion/engine[0]/thrust-lbs lt 1.0
       simulation/sim-time-sec gt 1.0
       velocities/h-dot-fps gt 0.0
     </test>
   </switch>
   <kinematic name="fcs/airbrake-actuator">''')
        p.write_text(s)
    sweep.run_case()
    with open(sweep.CASE_CSV) as fh:
        r = list(csv.reader(fh))
    h = [x.strip() for x in r[0]]
    col = lambda k: h.index(next(c for c in h if c.startswith(k)))
    rows = [[float(v) for v in x] for x in r[1:]]
    i_q, i_a, i_b, i_m, i_t = (col("Dynamic pressure"), col("Alpha"), col("Brake drag"),
                               col("Mach"), col("Time"))
    i_be = col("Beta")
    i_h = col("Altitude AGL")
    apo_i = max(range(len(rows)), key=lambda j: rows[j][i_h])
    asc = rows[:apo_i]
    cna_f0, fac = fin_share()
    fin = []
    for x in asc:
        a_tot = math.radians(math.hypot(x[i_a], x[i_be]))
        a_tot = min(a_tot, math.radians(15.0))
        fin.append(x[i_q] * PSF * 0.0845 * FT2 * cna_f0 * fac(x[i_m]) * a_tot)
    jb = max(range(len(asc)), key=lambda j: asc[j][i_b])
    jf = max(range(len(asc)), key=lambda j: fin[j])
    f = B.FIN
    arm = f["S"] / 3 * (f["Cr"] + 2 * f["Ct"]) / (f["Cr"] + f["Ct"]) + B.R * 0
    out = dict(case=name,
               brake_total_N=asc[jb][i_b] * LBF, brake_per_flap_N=asc[jb][i_b] * LBF / N_FLAPS,
               brake_t_s=asc[jb][i_t], brake_mach=asc[jb][i_m], brake_q_kPa=asc[jb][i_q] * PSF / 1e3,
               fin_set_N=fin[jf], fin_per_fin_N=fin[jf] / 2, fin_root_moment_Nm=fin[jf] / 2 * arm,
               fin_t_s=asc[jf][i_t], fin_mach=asc[jf][i_m], fin_q_kPa=asc[jf][i_q] * PSF / 1e3,
               fin_alpha_deg=math.hypot(asc[jf][i_a], asc[jf][i_be]),
               max_q_kPa=max(x[i_q] for x in asc) * PSF / 1e3)
    print({k: (round(v, 2) if isinstance(v, float) else v) for k, v in out.items()}, flush=True)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default=".")
    ap.add_argument("--json", default="loads_v6.json")
    a = ap.parse_args()
    for k, v in sweep._reroot(Path(a.root)).items():
        setattr(sweep, k, v)
    cases = [run("nominal", {}),
             run("impulse +10%", {"prop_mass": 1.10}),
             run("impulse +20%", {"prop_mass": 1.20}),
             run("impulse +30%", {"prop_mass": 1.30}),
             run("full brake at burnout (fault)", {}, force_brake=True),
             run("gust 40 ft/s at 2.5 s", {"gust": 40.0}),
             run("crosswind 30 ft/s", {"wind": 30.0})]
    json.dump(dict(fin=B.FIN, n_flaps=N_FLAPS, cases=cases), open(a.json, "w"), indent=1)
    print("wrote", a.json)


if __name__ == "__main__":
    main()
