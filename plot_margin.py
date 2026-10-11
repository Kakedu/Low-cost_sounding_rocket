# Hit rate against open-loop margin for the different flight computers.

import argparse
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import batch2_tables as bt
from plot_results import INK2

ap = argparse.ArgumentParser()
ap.add_argument("--dir", default="results_v6")
ap.add_argument("--open-nominal", type=float, default=10554.5, dest="open_nominal")
a = ap.parse_args()
bt.D = a.dir
T = (9600, 9800, 10000, 10200, 10400)
SETS = (("XML controller, true state, design Cd", "mc_closed_v6.csv", "mc_closed_t{t}_v6.csv", "#eb6834", "o"),
        ("3-state KF, ideal sensors", "mc_sil_kf_est_ideal_v6.csv", "mc_sil_kf_est_ideal_t{t}_v6.csv", "#2a78d6", "s"),
        ("4-state KF, sensor errors, temperature known", "mc_sil_kf4_est_notemp_v6.csv",
         "mc_sil_kf4_est_notemp_t{t}_v6.csv", "#1baf7a", "^"),
        ("4-state KF, all sensor errors", "mc_sil_kf4_est_v6.csv", "mc_sil_kf4_est_t{t}_v6.csv", "#7a5cc4", "D"))

fig, axs = plt.subplots(1, 2, figsize=(7.2, 3.6), sharex=True)
for lab, base, pat, col, mk in SETS:
    xs, pin, lin, hin, plo, llo, hlo = [], [], [], [], [], [], []
    for t in T:
        rows = bt.load(base if t == 10000 else pat.format(t=t))
        if not rows:
            continue
        s = bt.stats(rows, t)
        xs.append(100 * (a.open_nominal - t) / t)
        pin.append(100 * s["p_in"]); lin.append(100 * (s["p_in"] - s["ci_in"][0])); hin.append(100 * (s["ci_in"][1] - s["p_in"]))
        plo.append(100 * s["p_lo"]); llo.append(100 * (s["p_lo"] - s["ci_lo"][0])); hlo.append(100 * (s["ci_lo"][1] - s["p_lo"]))
    if xs:
        axs[0].errorbar(xs, pin, yerr=[lin, hin], color=col, marker=mk, ms=5, capsize=3, label=lab)
        axs[1].errorbar(xs, plo, yerr=[llo, hlo], color=col, marker=mk, ms=5, capsize=3, label=lab)
for y in (90, 95):
    axs[0].axhline(y, color=INK2, lw=0.6, ls=":")
axs[0].set_ylabel("flights within ±100 ft of target (%)")
axs[1].set_ylabel("flights more than 100 ft short (%)")
for ax in axs:
    ax.set_xlabel("open-loop margin above target (%)")
axs[0].set_ylim(0, 100); axs[1].set_ylim(0, 60)
h, l = axs[0].get_legend_handles_labels()
fig.legend(h, l, loc="lower center", ncol=2, fontsize=8)
fig.tight_layout(rect=(0, 0.13, 1, 1))
out = os.path.join(a.dir, "fig_margin_v6.png")
fig.savefig(out, dpi=200)
print("wrote", out)
