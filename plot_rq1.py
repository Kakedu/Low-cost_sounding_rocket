# Point-mass error vs launch angle and wind.

import csv
import os
import sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from plot_results import INK2

C = {"head": "#2a78d6", "tail": "#eb6834", "cross": "#1baf7a"}
D = sys.argv[1] if len(sys.argv) > 1 else "results_v6"
SRC = sys.argv[2] if len(sys.argv) > 2 else "anglesweep_v6.csv"
OUT = sys.argv[3] if len(sys.argv) > 3 else "fig_rq1_angle_wind.png"
rows = list(csv.DictReader(open(os.path.join(D, SRC))))
fig, axs = plt.subplots(1, 2, figsize=(7.2, 3.9), sharey=True)
for ax, mph in zip(axs, (10, 20)):
    still = [r for r in rows if r["wind"] == "none"]
    ax.plot([float(r["elev"]) for r in still],
            [100 * (float(r["pm_apogee_ft"]) / float(r["jsb_apogee_ft"]) - 1) for r in still],
            color=INK2, lw=1.2, ls="--", label="still air")
    for w, mk in (("head", "o"), ("tail", "s"), ("cross", "^")):
        s = [r for r in rows if r["wind"] == w and int(r["wind_mph"]) == mph]
        ax.plot([float(r["elev"]) for r in s],
                [100 * (float(r["pm_apogee_ft"]) / float(r["jsb_apogee_ft"]) - 1) for r in s],
                color=C[w], marker=mk, ms=5, label=f"{w} wind")
    ax.axhline(0, color=INK2, lw=0.8)
    ax.set_title(f"{mph} mph wind", fontsize=10)
    ax.set_xlabel("launch elevation (deg)")
axs[0].set_ylabel("point-mass apogee error vs 6-DOF (%)")
h, l = axs[0].get_legend_handles_labels()
fig.legend(h, l, loc="lower center", ncol=4, fontsize=8.5)
fig.tight_layout(rect=(0, 0.08, 1, 1))
fig.savefig(os.path.join(D, OUT), dpi=200)
