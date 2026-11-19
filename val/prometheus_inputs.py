# Prometheus inputs: Barrowman aero, mass properties and the site weather profile.

import json, math, datetime, os
import numpy as np
import netCDF4

L = 2.229; R = 0.06985; D = 2 * R
Ln = 0.742
N, S, Cr, Ct, Xt = 3, 0.130, 0.268, 0.136, 0.066
fin_le_tail = 0.273
Aref = math.pi * R * R

cna_nose = 2.0
x_nose = 0.500 * Ln
Xb = L - fin_le_tail
Lf = math.sqrt(S ** 2 + (Xt + Ct / 2 - Cr / 2) ** 2)
K = 1 + R / (S + R)
cna_fin0 = K * 4 * N * (S / D) ** 2 / (1 + math.sqrt(1 + (2 * Lf / (Cr + Ct)) ** 2))
x_fin = Xb + Xt / 3 * (Cr + 2 * Ct) / (Cr + Ct) + (Cr + Ct - Cr * Ct / (Cr + Ct)) / 6
Af = 0.5 * (Cr + Ct) * S
AR = 2 * S ** 2 / Af
gam_c = math.atan((Xt + Ct / 2 - Cr / 2) / S)
def dieder(M):
    b = math.sqrt(max(1e-6, abs(1 - M * M)))
    return AR / (2 + math.sqrt(4 + (b * AR / math.cos(gam_c)) ** 2))
machs = [0.0, 0.3, 0.6, 0.8, 0.9, 1.0, 1.2, 1.5]
cna, cp = [], []
for M in machs:
    Me = min(M, 0.8)
    cf = cna_fin0 * dieder(Me) / dieder(0.0)
    tot = cna_nose + cf
    xcp = (cna_nose * x_nose + cf * x_fin) / tot
    cna.append((M, tot)); cp.append((M, L - xcp))

m_body, x_body, I_body = 13.93, 0.9549, 4.87
m_mdry, x_m = 2.981, 0.274
r_c, h_c = 0.064, 0.548
m_prop = 3.737
m_dry = m_body + m_mdry
x_cg = (m_body * x_body + m_mdry * x_m) / m_dry
I_case = m_mdry * (3 * r_c ** 2 + h_c ** 2) / 12
I_dry = I_body + m_body * (x_body - x_cg) ** 2 + I_case + m_mdry * (x_m - x_cg) ** 2

LAT, LON = 32.939377, -106.911986
d = netCDF4.Dataset(os.path.join(os.path.dirname(os.path.abspath(__file__)), "data", "spaceport_america_pressure_levels_2023_hourly.nc"))
lats, lons = d["latitude"][:], d["longitude"][:]
lon_q = LON % 360 if lons.max() > 180 else LON
times = netCDF4.num2date(d["time"][:], d["time"].units)
target = datetime.datetime(2023, 6, 24, 15, 0)
it = int(np.argmin([abs((datetime.datetime(t.year, t.month, t.day, t.hour) - target).total_seconds()) for t in times]))
def bilin(var):
    a = d[var][it]
    i = np.searchsorted(-lats, -LAT) if lats[0] > lats[-1] else np.searchsorted(lats, LAT)
    j = np.searchsorted(lons, lon_q)
    i0, i1, j0, j1 = i - 1, i, j - 1, j
    wy = (LAT - lats[i0]) / (lats[i1] - lats[i0]); wx = (lon_q - lons[j0]) / (lons[j1] - lons[j0])
    return ((1 - wy) * ((1 - wx) * a[:, i0, j0] + wx * a[:, i0, j1]) + wy * ((1 - wx) * a[:, i1, j0] + wx * a[:, i1, j1]))
z = bilin("z") / 9.80665; T = bilin("t"); u = bilin("u"); v = bilin("v")
p = d["level"][:] * 100.0
order = np.argsort(z)
z, T, u, v, p = z[order], T[order], u[order], v[order], p[order]
hs = np.arange(1401, 6001, 200)
prof = [(float(h), float(np.interp(h, z, u)), float(np.interp(h, z, v)), float(np.interp(h, z, T)),
         float(np.exp(np.interp(h, z, np.log(p)))), 0.0) for h in hs]
prof = [(h, uu, vv, TT, pp, pp / (287.05 * TT)) for h, uu, vv, TT, pp, _ in prof]

out = dict(source="computed without third-party simulators", weather_time_utc=str(times[it]),
           mass_dry=m_dry, mass0=m_dry + m_prop, com_dry=x_cg, I_dry=I_dry,
           cna_by_mach=cna, cp_by_mach=cp, cna_fin0=cna_fin0, x_fin_from_tip=x_fin,
           static_margin_cal=(x_cg - cp[0][1]) / D, profile=prof)
json.dump(out, open("prometheus_inputs.json", "w"), indent=1)
print(f"weather {times[it]}  dry {m_dry:.3f} kg  CG {x_cg:.4f} m from tail  Iyy {I_dry:.3f} kg m2")
print("CNa:", [f"{c:.2f}" for _, c in cna]); print("CP (m from tail):", [f"{c:.4f}" for _, c in cp])
print(f"static margin {out['static_margin_cal']:.2f} cal")
for r_ in prof[::4]: print(["%.2f" % x for x in r_])
