# Compares the model with the Prometheus barometric and GPS apogee.

import csv, json, math, os
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({"figure.facecolor": "#fcfcfb", "axes.facecolor": "#fcfcfb", "savefig.facecolor": "#fcfcfb",
    "axes.edgecolor": "#52514e", "xtick.color": "#52514e", "ytick.color": "#52514e", "font.size": 10,
    "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True, "grid.color": "#e4e3df",
    "legend.frameon": False})
HERE = os.path.dirname(os.path.abspath(__file__))
D = os.path.join(HERE, "data")
LOGS = {"TeleMetrum": "2022-06-24-serial-5115-flight-0001-TeleMetrum.csv",
        "TeleMega": "2022-06-24-serial-6583-flight-0003-TeleMega.csv"}
PAD_M = 1401.0
U_GPS_M = 10.0
U_NUM_M = float(os.environ.get("U_NUM_M", "1.0"))


def read(fn):
    r = list(csv.reader(open(os.path.join(D, fn))))
    h = [c.strip() for c in r[0]]
    it, ih, ip = h.index("time"), h.index("height"), h.index("pressure")
    ig = [i for i, c in enumerate(h) if c == "altitude"][-1]
    raw = [(float(x[it]), float(x[ih]), float(x[ig]), float(x[ip])) for x in r[1:] if len(x) == len(r[0])]
    pad = raw[0][2]
    p_pad = sum(x[3] for x in raw[:5]) / 5
    fixes, prev = [], None
    for t, _, g, _ in raw:
        if g != prev:
            fixes.append((t, g - pad))
            prev = g
    rows = [(t, hb, g - pad) for t, hb, g, _ in raw if t <= 40]
    return rows, p_pad, fixes, [(t, p) for t, _, _, p in raw]


def parabola_peak(fixes, half=4.0):
    tm, _ = max((f for f in fixes if f[0] > 10), key=lambda f: f[1])
    pts = [f for f in fixes if abs(f[0] - tm) <= half]
    import numpy as np
    c = np.polyfit([p[0] - tm for p in pts], [p[1] for p in pts], 2)
    t0 = -c[1] / (2 * c[0])
    return float(np.polyval(c, t0)), tm + t0, len(pts)


def column_dT(p_pad, p_apo, dz, z_pad=PAD_M):
    from flight_temperature import column_temperature_K
    return column_temperature_K(p_pad,p_apo,dz,z_pad)


from validation_atmosphere import baro_reading



js = json.load(open(os.path.join(HERE, "jsb_result.json")))
mc = json.load(open(os.path.join(HERE, "prom_mc.json")))
out = {"jsbsim_apogee_m": js["apogee"], "jsbsim_dT_K": js.get("dT")}
fig, ax = plt.subplots(figsize=(6.4, 4.0))
style = {"TeleMetrum": ("#0b0b0b", "-"), "TeleMega": ("#eb6834", "--")}
dTs = []
for name, fn in LOGS.items():
    d, p_pad, fixes, pres = read(fn)
    gps_apo, t_apo, n_fit = parabola_peak(fixes)
    p_apo = min(p for t, p in pres if abs(t - t_apo) < 1.0)
    dT = column_dT(p_pad, p_apo, gps_apo)
    dTs.append(dT)
    out[f"{name}_baro_apogee_m"] = max(b for _, b, _ in d)
    out[f"{name}_gps_highest_fix_m"] = max(g for t, g in fixes if t > 10)
    out[f"{name}_gps_apogee_fit_m"] = gps_apo
    out[f"{name}_column_dT_K"] = dT
    c, ls = style[name]
    ax.plot([t for t, *_ in d], [b for _, b, _ in d], color=c, ls=ls, lw=1.6, label=f"{name}, barometric")
    ax.plot([t for t, *_ in d], [g for *_, g in d], color=c, ls=ls, lw=1.0, alpha=0.55, label=f"{name}, GPS")
gps = (out["TeleMetrum_gps_apogee_fit_m"] + out["TeleMega_gps_apogee_fit_m"]) / 2
baro = (out["TeleMetrum_baro_apogee_m"] + out["TeleMega_baro_apogee_m"]) / 2
dT_flight = sum(dTs) / len(dTs)
S, D_ = js["apogee"], gps
u_input = mc["sd"]
u_val = math.sqrt(u_input ** 2 + U_GPS_M ** 2 + U_NUM_M ** 2)
out.update(gps_mean_m=gps, baro_mean_m=baro, column_dT_K=dT_flight,
           jsbsim_as_baro_m=baro_reading(S, js["dT"], PAD_M),
           E_m=S - D_, E_pct=100 * (S / D_ - 1), u_input_m=u_input, u_D_m=U_GPS_M, u_num_m=U_NUM_M,
           u_val_m=u_val, E_over_u_val=(S - D_) / u_val,
           expanded_validation_halfwidth_m=2 * u_val, expanded_validation_halfwidth_pct=200 * u_val / D_,
           discrepancy_interval_low_m=S-D_-2*u_val, discrepancy_interval_high_m=S-D_+2*u_val,
           interpretation="Conditional expanded uncertainty interval, not a pass/fail or detection threshold",
           mc_mean_m=mc["mean"], mc_p2_5=mc.get("p2_5"), mc_p97_5=mc.get("p97_5"),
           flight_in_mc_95=int(mc.get("p2_5", -1e9) <= D_ <= mc.get("p97_5", 1e9)))
out["err_as_baro_vs_baro_pct"] = 100 * (out["jsbsim_as_baro_m"] / baro - 1)
from flight_temperature import infer_from_logs
out["flight_temperature_diagnostic"] = infer_from_logs(D)
tj = [t for t in js["t"] if t <= 40]
ax.plot(tj, js["h"][:len(tj)], color="#2a78d6", lw=2, label="JSBSim (geometric)")
ax.plot(tj, [baro_reading(h, js["dT"], PAD_M) for h in js["h"][:len(tj)]],
        color="#2a78d6", lw=1.3, ls="--", label="JSBSim, barometric model")
calibration_path = os.path.join(HERE, "..", "results_revised", "validation", "prometheus_2022_flight_calibrated", "result.json")
if os.path.exists(calibration_path):
    calibration = json.load(open(calibration_path))
    pairs = [(t, h) for t, h in zip(calibration["t"], calibration["h"]) if t <= 40]
    ax.plot([t for t, h in pairs], [h for t, h in pairs], color="#666666", lw=1.2, ls=":",
            label="Flight-informed temperature, same 2022 wind")
ax.axvline(29.6, color="#999999", lw=0.8, ls="--", label="Ejection pressure spike")
ax.set_xlim(0, 40); ax.set_ylim(0, 4500)
ax.set_xlabel("time from launch detection (s)"); ax.set_ylabel("height above pad (m)")
ax.legend(loc="lower right", fontsize=7.5, ncol=2)
fig.tight_layout(); fig.savefig(os.path.join(HERE, "fig_validation_gps.png"), dpi=200)
json.dump(out, open(os.path.join(HERE, "gps_check.json"), "w"), indent=1)
for k, v in out.items():
    print(f"  {k:30s} {v:10.2f}" if isinstance(v, float) else f"  {k:30s} {v}")
