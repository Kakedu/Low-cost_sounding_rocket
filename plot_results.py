# Recovery range, design map and Monte Carlo figures.

import csv
import math
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

CLOSED = "#2a78d6"
OPEN = "#eb6834"
INK = "#0b0b0b"
INK2 = "#52514e"
GRID = "#e4e3df"
SURF = "#fcfcfb"
TARGET = 10000.0

plt.rcParams.update({
    "figure.facecolor": SURF, "axes.facecolor": SURF, "savefig.facecolor": SURF,
    "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2,
    "ytick.color": INK2, "text.color": INK, "font.size": 10,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.8,
    "legend.frameon": False, "lines.linewidth": 2.0,
})


import sys as _sys
DIR = _sys.argv[1] if __name__ == "__main__" and len(_sys.argv) > 1 else "results_v6"


def _p(name):
    return os.path.join(DIR, name)


def read(path):
    with open(path) as f:
        return [{k: (float(v) if v not in ("",) else math.nan) for k, v in r.items()}
                for r in csv.DictReader(f)]


def fig_recovery():
    r = read(_p("recovery_v6.csv"))
    f = [x["impulse_factor"] * 100 - 100 for x in r]
    fig, ax = plt.subplots(figsize=(6.4, 4.0))
    ax.axhline(TARGET, color=INK2, lw=1, ls="--")
    ax.text(f[0], TARGET + 120, "target 10,000 ft", color=INK2, fontsize=9)
    ax.plot(f, [x["open_apogee_ft"] for x in r], color=OPEN, marker="o", ms=5,
            label="open loop (no brake)")
    ax.plot(f, [x["closed_apogee_ft"] for x in r], color=CLOSED, marker="s", ms=5,
            label="closed loop (air brake)")
    ok = [fi for fi, x in zip(f, r) if abs(x["error_ft"]) <= 100]
    ax.axvspan(min(ok), max(ok), color=CLOSED, alpha=0.07, lw=0)
    ax.text((min(ok) + max(ok)) / 2, 8750, f"within 100 ft:\n{min(ok):+.0f}% to {max(ok):+.0f}%",
            ha="center", color=INK2, fontsize=9)
    ax.set_xlabel("total impulse error (%)")
    ax.set_ylabel("apogee (ft)")
    ax.legend(loc="upper left")
    fig.tight_layout()
    fig.savefig(_p("fig_recovery.png"), dpi=200)
    plt.close(fig)


def fig_dmap():
    r = read(_p("dmap_v6.csv"))
    dcds = sorted({x["dcd"] for x in r})
    kps = sorted({x["kp_factor"] for x in r})
    Z = np.array([[next(x["worst"] for x in r if x["dcd"] == d and x["kp_factor"] == k)
                   for k in kps] for d in dcds])
    fig, ax = plt.subplots(figsize=(6.0, 4.2))
    ax.grid(False)
    from matplotlib.colors import LinearSegmentedColormap, LogNorm
    cmap = LinearSegmentedColormap.from_list(
        "blue", ["#cde2fb", "#86b6ef", "#3987e5", "#1c5cab", "#0d366b"])
    im = ax.imshow(Z, cmap=cmap, norm=LogNorm(vmin=20, vmax=2000), origin="lower",
                   aspect="auto")
    for i, d in enumerate(dcds):
        for j, k in enumerate(kps):
            v = Z[i, j]
            ax.text(j, i, f"{v:.0f}", ha="center", va="center", fontsize=9,
                    color=INK if v < 150 else "#ffffff")
    i0, j0 = dcds.index(2.0), kps.index(1.0)
    ax.add_patch(plt.Rectangle((j0 - 0.5, i0 - 0.5), 1, 1, fill=False, ec=OPEN, lw=2))
    ax.text(j0, i0 - 0.32, "current design", color=INK, fontsize=7.5, ha="center")
    ax.set_xticks(range(len(kps)), [f"{k * 0.001:g}" for k in kps])
    ax.set_yticks(range(len(dcds)), [f"{d:g}" for d in dcds])
    ax.set_xlabel("proportional gain kp (per ft)")
    ax.set_ylabel("brake drag increment at full deployment, ΔCd")
    cb = fig.colorbar(im, ax=ax)
    cb.set_label("worst |apogee error| over 6 cases (ft)")
    fig.tight_layout()
    fig.savefig(_p("fig_design_map.png"), dpi=200)
    plt.close(fig)


def _ellipse(N, E, k):
    N, E = np.asarray(N), np.asarray(E)
    C = np.cov(np.vstack([E, N]))
    w, V = np.linalg.eigh(C)
    t = np.linspace(0, 2 * np.pi, 200)
    pts = V @ (k * np.sqrt(w)[:, None] * np.vstack([np.cos(t), np.sin(t)]))
    return pts[0] + E.mean(), pts[1] + N.mean()


def fig_mc():
    o, c = read(_p("mc_open_v6.csv")), read(_p("mc_closed_v6.csv"))
    fig, ax = plt.subplots(figsize=(6.4, 3.8))
    bins = np.arange(8800, 12000, 50)
    ax.hist([x["apogee"] for x in o], bins=bins, color=OPEN, alpha=0.85,
            label="open loop", edgecolor=SURF, linewidth=0.5)
    ax.hist([x["apogee"] for x in c], bins=bins, color=CLOSED, alpha=0.85,
            label="closed loop", edgecolor=SURF, linewidth=0.5)
    ax.axvline(TARGET, color=INK2, lw=1, ls="--")
    ax.set_xlabel("apogee (ft)")
    ax.set_ylabel("flights (of %d)" % len(o))
    ax.legend(loc="upper right")
    fig.tight_layout()
    fig.savefig(_p("fig_mc_apogee.png"), dpi=200)
    plt.close(fig)
    fig, ax = plt.subplots(figsize=(6.0, 6.0))
    for rows, col, lab, mk in ((o, OPEN, "open loop", "o"), (c, CLOSED, "closed loop", "s")):
        N = [x["land_n"] for x in rows]
        E = [x["land_e"] for x in rows]
        ax.scatter(E, N, s=9, color=col, alpha=0.45, marker=mk, linewidths=0, label=lab)
        for k, ls in ((1, "-"), (2, "--")):
            ex, ny = _ellipse(N, E, k)
            ax.plot(ex, ny, color=col, lw=1.6, ls=ls)
    ax.plot(0, 0, marker="^", color=INK, ms=9)
    ax.annotate("pad", (0, 0), xytext=(8, -12), textcoords="offset points", color=INK)
    ax.set_aspect("equal")
    ax.set_xlabel("east (m)")
    ax.set_ylabel("north (m)")
    lg = ax.legend(loc="upper right", title="ellipses: 1σ solid, 2σ dashed",
                   title_fontsize=8, markerscale=2.5)
    for h in lg.legend_handles:
        h.set_alpha(1)
    fig.tight_layout()
    fig.savefig(_p("fig_mc_landing.png"), dpi=200)
    plt.close(fig)
    fig, ax = plt.subplots(figsize=(6.4, 3.8))
    ax.scatter([x["wind_fps"] for x in o], [x["land_dist"] for x in o], s=10,
               color=OPEN, alpha=0.5, linewidths=0, label="open loop")
    ax.scatter([x["wind_fps"] for x in c], [x["land_dist"] for x in c], s=10,
               color=CLOSED, alpha=0.5, marker="s", linewidths=0, label="closed loop")
    ax.set_xlabel("wind speed (fps)")
    ax.set_ylabel("landing distance from pad (m)")
    ax.legend(loc="upper left")
    fig.tight_layout()
    fig.savefig(_p("fig_mc_wind.png"), dpi=200)
    plt.close(fig)


if __name__ == "__main__":
    if os.path.exists(_p("recovery_v6.csv")):
        fig_recovery()
    if os.path.exists(_p("dmap_v6.csv")):
        fig_dmap()
    if os.path.exists(_p("mc_open_v6.csv")) and os.path.exists(_p("mc_closed_v6.csv")):
        fig_mc()
    print("done")
