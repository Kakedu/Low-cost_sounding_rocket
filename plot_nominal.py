# Nominal flight plots, open and closed loop.

import csv
import os
import sys
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from plot_results import CLOSED, OPEN, INK2, TARGET

def load(fn):
    r = list(csv.reader(open(fn))); h = [x.strip() for x in r[0]]
    cols = {k: [float(x[j]) for x in r[1:]] for j, k in enumerate(h)}
    return cols

D = sys.argv[1] if len(sys.argv) > 1 else str(Path(__file__).resolve().parent/'results_revised')
c, o = load(os.path.join(D, "nominal_closed_v6.csv")), load(os.path.join(D, "nominal_open_v6.csv"))
H, T, B = "Altitude AGL (ft)", "Time", "Airbrake position"
fig, (a1, a2) = plt.subplots(2, 1, figsize=(6.4, 5.2), sharex=True,
                             gridspec_kw={"height_ratios": [3, 1.2]})
k = lambda d: [i for i, t in enumerate(d[T]) if t <= 30]
for d, col, lab in ((o, OPEN, "open loop (no brake)"), (c, CLOSED, "closed loop (air brake)")):
    idx = k(d)
    a1.plot([d[T][i] for i in idx], [d[H][i] for i in idx], color=col, label=lab)
a1.axhline(TARGET, color=INK2, lw=1, ls="--")
a1.text(5.0, TARGET + 150, "target 10,000 ft AGL", color=INK2, fontsize=10)
t_apo = lambda d: d[T][max(range(len(d[H])), key=d[H].__getitem__)]
thrust = o["Thrust (lbf)"]
last_powered = max(i for i, v in enumerate(thrust) if v > 1.0)
t_burnout = o[T][last_powered + 1]
a1.axvline(t_burnout, color=INK2, lw=0.8, ls=":")
a1.text(t_burnout + 0.17, 1200, f"burnout {t_burnout:.2f} s", color=INK2, fontsize=10)
a1.annotate(f"{max(o[H]):,.0f} ft", (t_apo(o), max(o[H])), xytext=(-70, 6), textcoords="offset points", color=INK2, fontsize=10)
a1.annotate(f"{max(c[H]):,.0f} ft", (t_apo(c), max(c[H])), xytext=(-70, -16), textcoords="offset points", color=INK2, fontsize=10)
a1.set_ylabel("altitude AGL (ft)"); a1.set_ylim(0, 11800); a1.legend(loc="lower right")
idx = k(c)
a2.plot([c[T][i] for i in idx], [c[B][i] for i in idx], color=CLOSED)
a2.set_ylabel("brake deployment"); a2.set_ylim(0, 1.0); a2.set_xlabel("time from ignition (s)")
fig.tight_layout(); fig.savefig(os.path.join(D, "fig_nominal.png"), dpi=200)
