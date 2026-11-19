# Compares the Cavour model with the barometric log.

import csv, json, math, os
HERE = os.path.dirname(os.path.abspath(__file__))
PAD_M = 150.0
U_D_M = 5.0
U_NUM_M = float(os.environ.get("U_NUM_M", "1.0"))


from validation_atmosphere import baro_reading



if __name__ == "__main__":
    import numpy as np
    rows = [(float(r["ts"]), float(r["altitude[m]"])) for r in csv.DictReader(open(os.path.join(HERE, "data", "cavour", "altimeter_cavour.csv")))]
    tm, hm = max(rows, key=lambda r: r[1])
    pts = [r for r in rows if abs(r[0] - tm) <= 3.0]
    c = np.polyfit([p[0] - tm for p in pts], [p[1] for p in pts], 2)
    t0 = -c[1] / (2 * c[0]); D = float(np.polyval(c, t0))
    info = json.load(open(os.path.join(HERE, "build_info_cavour.json")))
    js = json.load(open(os.path.join(HERE, "jsb_result_cavour.json")))
    S_geo = js["apogee"]; S = baro_reading(S_geo, info["dT"])
    out = dict(record_max_m=hm, record_fit_m=D, t_fit_s=tm + t0, n_fit=len(pts),
               jsbsim_geometric_m=S_geo, jsbsim_t_apogee_s=js["t_apogee"], dT_K=info["dT"],
               jsbsim_as_baro_m=S, E_m=S - D, E_pct=100 * (S / D - 1),
               E_vs_official_pct=100 * (S / hm - 1))
    try:
        mc = json.load(open(os.path.join(HERE, "cavour_mc.json")))
        u_in = mc["baro_m"]["sd"]
        u_val = math.sqrt(u_in ** 2 + U_D_M ** 2 + U_NUM_M ** 2)
        out.update(u_input_m=u_in, u_D_m=U_D_M, u_val_m=u_val, E_over_u_val=(S - D) / u_val,
                   expanded_validation_halfwidth_m=2 * u_val, expanded_validation_halfwidth_pct=200 * u_val / D,
                   discrepancy_interval_low_m=S-D-2*u_val, discrepancy_interval_high_m=S-D+2*u_val,
                   interpretation="Conditional expanded uncertainty interval, not a pass/fail or detection threshold",
                   mc_baro_mean_m=mc["baro_m"]["mean"], mc_baro_p2_5=mc["baro_m"]["p2_5"],
                   mc_baro_p97_5=mc["baro_m"]["p97_5"],
                   flight_in_mc_95=int(mc["baro_m"]["p2_5"] <= D <= mc["baro_m"]["p97_5"]))
    except FileNotFoundError:
        pass
    json.dump(out, open(os.path.join(HERE, "cavour_check.json"), "w"), indent=1)
    for k, v in out.items():
        print(f"  {k:28s} {v:10.2f}" if isinstance(v, float) else f"  {k:28s} {v}")


def plot():
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"figure.facecolor": "#fcfcfb", "axes.facecolor": "#fcfcfb", "savefig.facecolor": "#fcfcfb",
        "axes.edgecolor": "#52514e", "xtick.color": "#52514e", "ytick.color": "#52514e", "font.size": 10,
        "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True, "grid.color": "#e4e3df",
        "legend.frameon": False})
    rows = [(float(r["ts"]), float(r["altitude[m]"])) for r in csv.DictReader(open(os.path.join(HERE, "data", "cavour", "altimeter_cavour.csv")))]
    info = json.load(open(os.path.join(HERE, "build_info_cavour.json")))
    js = json.load(open(os.path.join(HERE, "jsb_result_cavour.json")))
    tj = js["t"][::10]; hj = js["h"][::10]
    hb = [baro_reading(h, info["dT"], PAD_M) for h in hj]
    fig, ax = plt.subplots(figsize=(6.4, 4.0))
    ax.plot([r[0] for r in rows if r[0] <= 30], [r[1] for r in rows if r[0] <= 30], color="#0b0b0b", lw=1.6, label="Flight, barometric altimeter")
    ax.plot(tj, hj, color="#2a78d6", lw=1.2, ls="--", label="JSBSim, geometric")
    ax.plot(tj, hb, color="#2a78d6", lw=2, label="JSBSim, as a standard-atmosphere barometer reads it")
    ax.set_xlim(0, 30); ax.set_ylim(0, 3100)
    ax.set_xlabel("time from launch (s)"); ax.set_ylabel("height above pad (m)")
    ax.legend(loc="lower right", fontsize=8)
    fig.tight_layout(); fig.savefig(os.path.join(HERE, "fig_validation_cavour.png"), dpi=200)


if __name__ == "__main__":
    plot()
