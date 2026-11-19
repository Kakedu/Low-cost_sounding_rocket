# Builds the JSBSim model of Cavour (PoliTo, EuRoC 2023), second validation case.

import csv, datetime, json, math, os, shutil
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from atmosphere_model import troposphere_density_kg_m3, EARTH_RADIUS_M
from motor_model import thrust_at, burned_mass, write_motor
import numpy as np
import netCDF4

CD_SCALE = float(os.environ.get("CD_SCALE", "1.0"))
DT_EXTRA = float(os.environ.get("DT_EXTRA", "0.0"))
DT_ABS = os.environ.get("DT_ABS", "reanalysis")
WIND_SCALE = float(os.environ.get("WIND_SCALE", "1.0"))
WIND_ROT_DEG = float(os.environ.get("WIND_ROT_DEG", "0.0"))
IMP_SCALE = float(os.environ.get("IMP_SCALE", "1.0"))
HERE = os.path.dirname(os.path.abspath(__file__))
DD = os.path.join(HERE, "data", "cavour")
ROOT = os.environ.get("JSB_ROOT", os.path.join(HERE, "jsb_cavour"))
INFO = os.environ.get("INFO", os.path.join(HERE, "build_info_cavour.json"))
M2IN, KG2LB, M2FT, LBF, G0 = 39.3701, 2.20462, 3.28084, 4.4482216, 9.80665
for d in ("aircraft/cavour", "engine", "scripts", "data_output"):
    os.makedirs(f"{ROOT}/{d}", exist_ok=True)
import jsbsim
shutil.copy(f"{os.path.dirname(jsbsim.__file__)}/engine/direct.xml", f"{ROOT}/engine/direct.xml")

L = 2.7224; R = 0.052; D = 2 * R
Ln = 0.52
N, S, Cr, Ct = 4, 0.10, 0.20, 0.07
Xt = Cr - Ct
fin_le_tail = 0.2104
st = lambda x_rp: (L - x_rp) * M2IN

cna_nose, x_nose = 2.0, 0.500 * Ln
Xb = L - fin_le_tail
Lf = math.sqrt(S ** 2 + (Xt + Ct / 2 - Cr / 2) ** 2)
K = 1 + R / (S + R)
cna_fin0 = K * 4 * N * (S / D) ** 2 / (1 + math.sqrt(1 + (2 * Lf / (Cr + Ct)) ** 2))
x_fin = Xb + Xt / 3 * (Cr + 2 * Ct) / (Cr + Ct) + (Cr + Ct - Cr * Ct / (Cr + Ct)) / 6
Af = 0.5 * (Cr + Ct) * S; AR = 2 * S ** 2 / Af
gam_c = math.atan((Xt + Ct / 2 - Cr / 2) / S)
def dieder(M):
    b = math.sqrt(max(1e-6, abs(1 - M * M)))
    return AR / (2 + math.sqrt(4 + (b * AR / math.cos(gam_c)) ** 2))
cna_tab, cp_tab = [], []
for M in [0.0, 0.3, 0.6, 0.8, 0.9, 1.0, 1.2, 1.5]:
    cf = cna_fin0 * dieder(min(M, 0.8)) / dieder(0.0)
    tot = cna_nose + cf
    cna_tab.append((M, tot)); cp_tab.append((M, L - (cna_nose * x_nose + cf * x_fin) / tot))
cp_sub = float(np.mean([c for m, c in cp_tab if m <= 0.9]))

m_dry, com_dry, I_dry_kg = 8.219, 1.1994, 4.449
x_prop = 0.19966
eng = [l.split() for l in open(os.path.join(DD, "Cesaroni_3618L995-P.eng")) if l.strip() and not l.startswith(";")]
m_prop_kg = float(eng[0][4])
tc = [(0.0, 0.0)] + [(float(a), float(b)) for a, b in eng[1:]]
t_arr = np.array([t for t, _ in tc]); f_arr = np.array([f for _, f in tc])
I_tot = float(np.trapezoid(f_arr, t_arr))
isp = I_tot / (m_prop_kg * G0) * IMP_SCALE
def thrust(t): return thrust_at(tc, t, IMP_SCALE)
def burned(t): return burned_mass(tc, m_prop_kg, t, IMP_SCALE)

cmq_l_fin = com_dry - (fin_le_tail - 0.5 * Cr)
cmq_l_nose = (L - x_nose) - com_dry
cmq = -2 * ((cna_tab[0][1] - cna_nose) * (cmq_l_fin / D) ** 2 + cna_nose * (cmq_l_nose / D) ** 2)

raw = [tuple(map(float, r)) for r in csv.reader(open(os.path.join(DD, "drag_coefficient_power_off.csv")))]
raw = sorted({round(m, 3): c for m, c in raw}.items())
grid = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.85, 0.9, 0.95, 1.0, 1.1, 1.2, 1.5]
ms, cs = np.array([m for m, _ in raw]), np.array([c for _, c in raw])
cd_rows = [(g, float(np.interp(g, ms, cs))) for g in grid if g <= ms.max()]
def cd(m): return float(np.interp(m, ms, cs))

LAT, LON, PAD = 39.388692, -8.287814, 150.0
w = netCDF4.Dataset(os.path.join(DD, "euroc_2023_all_windows.nc"))
lats, lons = w["latitude"][:], w["longitude"][:]
times = netCDF4.num2date(w["time"][:], w["time"].units)
def col(it):
    i = int(np.argmin(abs(lats - LAT))); j = int(np.argmin(abs(lons - LON)))
    i0, i1 = sorted((i, i + (1 if lats[i] > LAT else -1)))
    j0, j1 = sorted((j, j + (1 if lons[j] < LON else -1)))
    wy = (LAT - lats[i0]) / (lats[i1] - lats[i0]); wx = (LON - lons[j0]) / (lons[j1] - lons[j0])
    def bil(v):
        a = w[v][it]
        return (1 - wy) * ((1 - wx) * a[:, i0, j0] + wx * a[:, i0, j1]) + wy * ((1 - wx) * a[:, i1, j0] + wx * a[:, i1, j1])
    return bil("z") / G0, bil("t"), bil("u"), bil("v")
t_launch = datetime.datetime(2023, 10, 13, 13, 30)
idx = [k for k, t in enumerate(times) if datetime.datetime(t.year, t.month, t.day, t.hour) in
       (datetime.datetime(2023, 10, 13, 13), datetime.datetime(2023, 10, 13, 14))]
cols = [col(k) for k in idx]
z, T, u, v = [np.mean([c[q] for c in cols], axis=0) for q in range(4)]
p = w["level"][:] * 100.0
z = EARTH_RADIUS_M * z / (EARTH_RADIUS_M - z)
rotation = math.radians(WIND_ROT_DEG)
u, v = u * math.cos(rotation) - v * math.sin(rotation), u * math.sin(rotation) + v * math.cos(rotation)
o = np.argsort(z); z, T, u, v, p = z[o], T[o], u[o], v[o], p[o]
hs = np.arange(PAD, 4501, 100.0)
prof = np.array([(h, np.interp(h, z, u), np.interp(h, z, v), np.interp(h, z, T),
                  math.exp(np.interp(h, z, np.log(p)))) for h in hs])
rho = prof[:, 4] / (287.05 * prof[:, 3])
def isa_rho(h, dT):
    return troposphere_density_kg_m3(h, dT)
sel = prof[:, 0] <= PAD + 3200
dTs = np.arange(-10, 30, 0.05)
dT_re = float(dTs[int(np.argmin([np.sum((isa_rho(prof[sel, 0], d) - rho[sel]) ** 2) for d in dTs]))])
dT = (dT_re if DT_ABS == "reanalysis" else float(DT_ABS)) + DT_EXTRA
T_pad_re = float(np.interp(PAD, z, T)) - 273.15

RAIL_LEN, INCL, HEAD, BTN = 12.0, 84.0, 133.0, 1.0954
RAIL = RAIL_LEN - BTN
t = s = vv = 0.0; dt = 1e-4; m = m_dry + m_prop_kg; rho0 = isa_rho(PAD, dT)
while s < RAIL:
    if t >= tc[-1][0] * IMP_SCALE: raise RuntimeError("Motor burned out before rail exit")
    F = thrust(t)
    a = (F - 0.5 * rho0 * vv * vv * math.pi * R * R * cd(vv / 340.0) * CD_SCALE) / m - G0 * math.sin(math.radians(INCL))
    vv = max(0.0, vv + a * dt); s += vv * dt; m -= F / (isp * G0) * dt; t += dt
t_exit, v_exit = t, vv
p_burn = burned(t_exit)

motor_info = write_motor(tc, m_prop_kg, f"{ROOT}/engine/L995.xml", name="L995",
                         start_time=t_exit, isp_scale=IMP_SCALE,
                         sample_dt=float(os.environ.get("MOTOR_SAMPLE_DT", ".002")))
prop_left_lb = motor_info["remaining_prop_lb"]

S_ft2 = math.pi * R * R * M2FT ** 2; d_ft = D * M2FT
tbl = lambda rr: "\n".join(f"              {a:6.2f}  {b:8.4f}" for a, b in rr)
wind_e = "\n".join(f"              {h*M2FT:9.1f}  {WIND_SCALE*uu*M2FT:8.3f}" for h, uu in zip(prof[:, 0], prof[:, 1]))
wind_n = "\n".join(f"              {h*M2FT:9.1f}  {WIND_SCALE*vn*M2FT:8.3f}" for h, vn in zip(prof[:, 0], prof[:, 2]))
gain = '<property>fcs/body-aero-gain</property>'
ac = f'''<?xml version="1.0"?>
<!-- Cavour, PoliTo Rocket Team, EuRoC 2023. Built by build_cavour.py. -->
<fdm_config name="cavour" version="2.0" release="BETA">
 <fileheader><author>generated</author><description>validation case 2</description></fileheader>
 <metrics>
   <wingarea unit="FT2"> {S_ft2:.5f} </wingarea><wingspan unit="FT"> {d_ft:.5f} </wingspan><chord unit="FT"> {d_ft:.5f} </chord>
   <htailarea unit="FT2">0</htailarea><htailarm unit="FT">0</htailarm><vtailarea unit="FT2">0</vtailarea><vtailarm unit="FT">0</vtailarm>
   <location name="AERORP" unit="IN"><x> {st(cp_sub):.3f} </x><y>0</y><z>0</z></location>
   <location name="EYEPOINT" unit="IN"><x>0</x><y>0</y><z>0</z></location>
   <location name="VRP" unit="IN"><x>0</x><y>0</y><z>0</z></location>
 </metrics>
 <mass_balance>
   <ixx unit="SLUG*FT2"> {0.014634*0.737562:.5f} </ixx>
   <iyy unit="SLUG*FT2"> {I_dry_kg*0.737562:.4f} </iyy><izz unit="SLUG*FT2"> {I_dry_kg*0.737562:.4f} </izz>
   <emptywt unit="LBS"> {m_dry*KG2LB:.4f} </emptywt>
   <location name="CG" unit="IN"><x> {st(com_dry):.3f} </x><y>0</y><z>0</z></location>
 </mass_balance>
 <ground_reactions>
  <contact type="STRUCTURE" name="NOSE"><location unit="IN"><x>0</x><y>0</y><z>0</z></location>
   <static_friction>0.8</static_friction><dynamic_friction>0.5</dynamic_friction>
   <spring_coeff unit="LBS/FT">1000</spring_coeff><damping_coeff unit="LBS/FT/SEC">100</damping_coeff></contact>
  <contact type="STRUCTURE" name="TAIL"><location unit="IN"><x>{st(0):.2f}</x><y>0</y><z>0</z></location>
   <static_friction>0.8</static_friction><dynamic_friction>0.5</dynamic_friction>
   <spring_coeff unit="LBS/FT">1000</spring_coeff><damping_coeff unit="LBS/FT/SEC">100</damping_coeff></contact>
 </ground_reactions>
 <propulsion>
  <engine file="L995"><feed>0</feed>
   <thruster file="direct"><location unit="IN"><x> {st(0):.2f} </x><y>0</y><z>0</z></location>
    <orient unit="DEG"><pitch>0</pitch><roll>0</roll><yaw>0</yaw></orient></thruster>
  </engine>
  <tank type="FUEL" number="0"><location unit="IN"><x> {st(x_prop):.2f} </x><y>0</y><z>0</z></location>
   <capacity unit="LBS"> {prop_left_lb:.15g} </capacity><contents unit="LBS"> {prop_left_lb:.15g} </contents></tank>
 </propulsion>
 <flight_control name="FCS">
  <property value="0.0">fcs/drogue-cmd</property>
  <channel name="Wind">
   <fcs_function name="fcs/wind-east-fps"><function><table><independentVar lookup="row">position/h-sl-ft</independentVar>
     <tableData>
{wind_e}
     </tableData></table></function><output>atmosphere/wind-east-fps</output></fcs_function>
   <fcs_function name="fcs/wind-north-fps"><function><table><independentVar lookup="row">position/h-sl-ft</independentVar>
     <tableData>
{wind_n}
     </tableData></table></function><output>atmosphere/wind-north-fps</output></fcs_function>
  </channel>
  <channel name="Recovery">
   <fcs_function name="fcs/body-aero-gain"><function><difference><value>1</value><property>fcs/drogue-cmd</property></difference></function></fcs_function>
  </channel>
 </flight_control>
 <aerodynamics>
  <axis name="DRAG">
   <function name="aero/force/Drag_mach"><description>RASAero II power-off Cd</description>
    <product><property>aero/qbar-psf</property><property>metrics/Sw-sqft</property><value> {CD_SCALE} </value>
     <table><independentVar lookup="row">velocities/mach</independentVar><tableData>
{tbl(cd_rows)}
     </tableData></table></product></function>
  </axis>
  <axis name="LIFT">
   <function name="aero/force/Lift_alpha"><product>{gain}<property>aero/qbar-psf</property><property>metrics/Sw-sqft</property>
     <sin><property>aero/alpha-rad</property></sin>
     <table><independentVar lookup="row">velocities/mach</independentVar><tableData>
{tbl(cna_tab)}
     </tableData></table></product></function>
  </axis>
  <axis name="SIDE">
   <function name="aero/force/Side_beta"><product>{gain}<property>aero/qbar-psf</property><property>metrics/Sw-sqft</property>
     <sin><property>aero/beta-rad</property></sin><value>-1</value>
     <table><independentVar lookup="row">velocities/mach</independentVar><tableData>
{tbl(cna_tab)}
     </tableData></table></product></function>
  </axis>
  <axis name="PITCH">
   <function name="aero/moment/Pitch_damp"><product>{gain}<property>aero/qbar-psf</property><property>metrics/Sw-sqft</property>
     <property>metrics/cbarw-ft</property><property>aero/ci2vel</property><property>velocities/q-aero-rad_sec</property><value> {cmq:.1f} </value></product></function>
  </axis>
  <axis name="YAW">
   <function name="aero/moment/Yaw_damp"><product>{gain}<property>aero/qbar-psf</property><property>metrics/Sw-sqft</property>
     <property>metrics/bw-ft</property><property>aero/bi2vel</property><property>velocities/r-aero-rad_sec</property><value> {cmq:.1f} </value></product></function>
  </axis>
  <axis name="ROLL">
   <function name="aero/moment/Roll_damp"><product><property>aero/qbar-psf</property><property>metrics/Sw-sqft</property>
     <property>metrics/bw-ft</property><property>aero/bi2vel</property><property>velocities/p-aero-rad_sec</property><value> -0.8 </value></product></function>
  </axis>
 </aerodynamics>
 <output name="data_output/cavour.csv" type="CSV" rate="50">
  <property caption="h_agl_ft"> position/h-agl-ft </property>
  <property caption="hdot_fps"> velocities/h-dot-fps </property>
  <property caption="mach"> velocities/mach </property>
  <property caption="theta_deg"> attitude/theta-deg </property>
  <property caption="thrust_lbf"> propulsion/engine[0]/thrust-lbs </property>
 </output>
</fdm_config>
'''
open(f"{ROOT}/aircraft/cavour/cavour.xml", "w").write(ac)
h0 = RAIL * math.sin(math.radians(INCL))
open(f"{ROOT}/aircraft/cavour/rail_exit.xml", "w").write(f'''<?xml version="1.0"?>
<initialize name="rail_exit">
  <ubody unit="FT/SEC"> {v_exit*M2FT:.3f} </ubody><vbody unit="FT/SEC"> 0.0 </vbody><wbody unit="FT/SEC"> 0.0 </wbody>
  <phi unit="DEG"> 0.0 </phi><theta unit="DEG"> {INCL} </theta><psi unit="DEG"> {HEAD} </psi>
  <latitude unit="DEG" type="geod"> {LAT} </latitude><longitude unit="DEG"> {LON} </longitude>
  <elevation unit="FT"> {PAD*M2FT:.2f} </elevation><altitudeAGL unit="FT"> {h0*M2FT:.2f} </altitudeAGL>
</initialize>
''')
open(f"{ROOT}/scripts/cavour.xml", "w").write(f'''<?xml version="1.0"?>
<runscript name="cavour validation">
  <use aircraft="cavour" initialize="rail_exit"/>
  <run start="0.0" end="120" dt="{os.environ.get("SIM_DT", "0.001")}">
    <event name="Ignition"><condition> simulation/sim-time-sec ge 0.0 </condition>
      <set name="atmosphere/delta-T" value="{dT*1.8:.3f}"/>
      <set name="propulsion/engine[0]/set-running" value="1"/>
      <set name="fcs/throttle-cmd-norm" value="1.0"/></event>
    <event name="Apogee"><condition> simulation/sim-time-sec gt 5.0
                  velocities/h-dot-fps lt 0.0 </condition>
      <set name="simulation/terminate" value="1"/></event>
  </run>
</runscript>
''')
info = dict(motor=motor_info, t_exit=t_exit, v_exit=v_exit, dT=dT, dT_reanalysis=dT_re, T_pad_reanalysis_C=T_pad_re,
            isp=isp, impulse=I_tot, m_prop=m_prop_kg, cmq=cmq, cna0=cna_tab[0][1],
            cp_sub_from_tail=cp_sub, static_margin_cal=(com_dry - cp_sub) / D,
            static_margin_launch_cal=((m_dry * com_dry + m_prop_kg * x_prop) / (m_dry + m_prop_kg) - cp_sub) / D,
            weather=[str(times[k]) for k in idx], pad_m=PAD)
json.dump(info, open(INFO, "w"), indent=1)
print(json.dumps(info, indent=1))
