# Builds the JSBSim model of Prometheus (Western Engineering, 2022) for validation.

import json, math, os, shutil
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from atmosphere_model import troposphere_density_kg_m3, EARTH_RADIUS_M
from motor_model import thrust_at, burned_mass, write_motor
import numpy as np
CD_SCALE = float(os.environ.get("CD_SCALE", "1.0"))
DT_EXTRA = float(os.environ.get("DT_EXTRA", "0.0"))
DT_ABS = os.environ.get("DT_ABS", "reanalysis")
WIND_SCALE = float(os.environ.get("WIND_SCALE", "1.0"))
WIND_ROT_DEG = float(os.environ.get("WIND_ROT_DEG", "0.0"))
IMP_SCALE = float(os.environ.get("IMP_SCALE", "1.0"))
HERE = os.path.dirname(os.path.abspath(__file__))
os.chdir(HERE)

M2IN, KG2LB, M2FT = 39.3701, 2.20462, 3.28084
LBF = 4.4482216
G0 = 9.80665
input_file = os.environ.get("PROM_INPUTS", "prometheus_inputs_2022.json")
o = json.load(open(input_file))
if "WEATHER_WEIGHT" in os.environ:
    hs = np.array(o["profile"])[:, 0]
    columns = []
    for stamp in ("2022062412", "2022062500"):
        a = np.loadtxt(os.path.join(HERE, "data", f"prometheus_sounding_{stamp}.csv"), delimiter=",", skiprows=1)
        columns.append(np.column_stack([hs] + [np.interp(hs, a[:,0], a[:,j]) for j in range(1,5)]))
    weight = float(os.environ["WEATHER_WEIGHT"])
    a = (1-weight)*columns[0]+weight*columns[1]
    o["profile"] = np.column_stack((a, a[:,4]/(287.05287*a[:,3]))).tolist()
rotation = math.radians(WIND_ROT_DEG)
for row in o["profile"]:
    u, v = row[1], row[2]
    row[1] = u*math.cos(rotation)-v*math.sin(rotation)
    row[2] = u*math.sin(rotation)+v*math.cos(rotation)
ROOT = os.environ.get("JSB_ROOT", os.path.join(HERE, "jsb"))
INFO = os.environ.get("INFO", os.path.join(HERE, "build_info.json"))
for d in ("aircraft/prometheus", "engine", "scripts", "data_output"):
    os.makedirs(f"{ROOT}/{d}", exist_ok=True)
import jsbsim
P = os.path.dirname(jsbsim.__file__)
shutil.copy(f"{P}/engine/direct.xml", f"{ROOT}/engine/direct.xml")

L = 2.229
st = lambda x_rp: (L - x_rp) * M2IN
R = 0.06985; D = 2 * R
S_ft2 = math.pi * R * R * M2FT ** 2
d_ft = D * M2FT
cg_dry = st(o["com_dry"])
x_motor = st(0.274)
x_noz = st(0.0)
cp_tab = o["cp_by_mach"]
cna_tab = o["cna_by_mach"]
cp_sub = np.mean([c for m, c in cp_tab if m <= 0.9])
x_cp = st(cp_sub)
I_dry = o["I_dry"] * 0.737562
m_dry_lb = o["mass_dry"] * KG2LB
m_prop_kg = 3.737

cna_nose, cna_tot = 2.0, cna_tab[0][1]
cna_fin = cna_tot - cna_nose
x_fin_rp = 0.273 - 0.5 * 0.268 + 0.10
x_nose_rp = 2.229 - 0.5 * 0.742
l_fin = o["com_dry"] - x_fin_rp
l_nose = x_nose_rp - o["com_dry"]
cmq = -2 * (cna_fin * (l_fin / D) ** 2 + cna_nose * (l_nose / D) ** 2)

eng = [l.split() for l in open(os.path.join(HERE, "data", "Cesaroni_7579M1520-P.eng"))
       if l.strip() and not l.startswith(";")][1:]
tc = [(0.0, 0.0)] + [(float(a), float(b)) for a, b in eng]
t_arr = np.array([t for t, _ in tc]); f_arr = np.array([f for _, f in tc])
I_tot = np.trapezoid(f_arr, t_arr)
isp = I_tot / (m_prop_kg * G0) * IMP_SCALE
def thrust(t): return thrust_at(tc, t, IMP_SCALE)
def burned(t): return burned_mass(tc, m_prop_kg, t, IMP_SCALE)

pts = [(0.15,0.422),(0.45,0.38),(0.77,0.32),(0.82,0.3),(0.88,0.3),(0.94,0.32),(0.99,0.37),
       (1.04,0.44),(1.24,0.43),(1.33,0.42),(1.49,0.39)]
def cd(m):
    if m <= 0.15: return 0.422
    for (m0,c0),(m1,c1) in zip(pts, pts[1:]):
        if m <= m1: return c0 + (m-m0)*(c1-c0)/(m1-m0)
    return 0.39
prof = np.array(o["profile"])
def isa_rho(h, dT):
    return troposphere_density_kg_m3(h, dT)
sel = prof[:, 0] <= 1401 + 4300
dTs = np.arange(-30, 60.001, 0.05)
err = [np.sum((isa_rho(prof[sel, 0], d) - prof[sel, 5]) ** 2) for d in dTs]
dT_reanalysis = float(dTs[int(np.argmin(err))])
dT = (float(DT_ABS) if DT_ABS != "reanalysis" else dT_reanalysis) + DT_EXTRA
rel = (isa_rho(prof[sel, 0], dT) / prof[sel, 5] - 1) * 100
print(f"delta-T {dT:.2f} K (profile density fit {dT_reanalysis:.2f} K), density error vs reanalysis {rel.min():+.1f}% to {rel.max():+.1f}%")

RAIL_LEN, INCL, HEAD = 5.18, 80.0, 75.0
RAIL = RAIL_LEN - 0.69
t = s = v = 0.0; dt = 1e-4; m = o["mass0"]; rho = isa_rho(1401., dT)
while s < RAIL:
    if t >= tc[-1][0] * IMP_SCALE: raise RuntimeError("Motor burned out before rail exit")
    F = thrust(t); mdot = F / (isp * G0)
    a = (F - 0.5 * rho * v * v * math.pi * R * R * cd(v / 340.) * CD_SCALE * 1.02) / m - G0 * math.sin(math.radians(INCL))
    v = max(0.0, v + a * dt); s += v * dt; m -= mdot * dt; t += dt
t_exit, v_exit = t, v
p_burn = burned(t_exit)
print(f"rail exit t {t_exit:.3f} s  v {v_exit:.2f} m/s")

motor_info = write_motor(tc, m_prop_kg, f"{ROOT}/engine/M1520.xml", name="M1520",
                         start_time=t_exit, isp_scale=IMP_SCALE,
                         sample_dt=float(os.environ.get("MOTOR_SAMPLE_DT", ".002")))
prop_left_lb = motor_info["remaining_prop_lb"]


wind_e = "\n".join(f"              {h*M2FT:9.1f}  {WIND_SCALE*u*M2FT:8.3f}" for h, u in zip(prof[:, 0], prof[:, 1]))
wind_n = "\n".join(f"              {h*M2FT:9.1f}  {WIND_SCALE*w*M2FT:8.3f}" for h, w in zip(prof[:, 0], prof[:, 2]))

def tbl(rows_, fmt="{:6.2f}  {:8.4f}"):
    return "\n".join("              " + fmt.format(a, b) for a, b in rows_)
cd_rows = [(0.0, 0.422)] + pts + [(3.0, 0.39)]
cna_rows = [(m_, c_) for m_, c_ in cna_tab]

ac = f'''<?xml version="1.0"?>
<!-- Prometheus, Western Engineering Rocketry Team, Spaceport America Cup 2022.
     Built by build_prometheus.py from the team's published definition. -->
<fdm_config name="prometheus" version="2.0" release="BETA">
 <fileheader><author>generated</author><description>validation case</description></fileheader>
 <metrics>
   <wingarea unit="FT2"> {S_ft2:.5f} </wingarea>
   <wingspan unit="FT"> {d_ft:.5f} </wingspan>
   <chord unit="FT"> {d_ft:.5f} </chord>
   <htailarea unit="FT2">0</htailarea><htailarm unit="FT">0</htailarm>
   <vtailarea unit="FT2">0</vtailarea><vtailarm unit="FT">0</vtailarm>
   <location name="AERORP" unit="IN"><x> {x_cp:.3f} </x><y>0</y><z>0</z></location>
   <location name="EYEPOINT" unit="IN"><x>0</x><y>0</y><z>0</z></location>
   <location name="VRP" unit="IN"><x>0</x><y>0</y><z>0</z></location>
 </metrics>
 <mass_balance>
   <ixx unit="SLUG*FT2"> 0.044 </ixx>
   <iyy unit="SLUG*FT2"> {I_dry:.4f} </iyy>
   <izz unit="SLUG*FT2"> {I_dry:.4f} </izz>
   <emptywt unit="LBS"> {m_dry_lb:.4f} </emptywt>
   <location name="CG" unit="IN"><x> {cg_dry:.3f} </x><y>0</y><z>0</z></location>
 </mass_balance>
 <ground_reactions>
  <contact type="STRUCTURE" name="NOSE"><location unit="IN"><x>0</x><y>0</y><z>0</z></location>
   <static_friction>0.8</static_friction><dynamic_friction>0.5</dynamic_friction>
   <spring_coeff unit="LBS/FT">1000</spring_coeff><damping_coeff unit="LBS/FT/SEC">100</damping_coeff></contact>
  <contact type="STRUCTURE" name="TAIL"><location unit="IN"><x>{x_noz:.2f}</x><y>0</y><z>0</z></location>
   <static_friction>0.8</static_friction><dynamic_friction>0.5</dynamic_friction>
   <spring_coeff unit="LBS/FT">1000</spring_coeff><damping_coeff unit="LBS/FT/SEC">100</damping_coeff></contact>
 </ground_reactions>
 <propulsion>
  <engine file="M1520">
   <feed>0</feed>
   <thruster file="direct">
    <location unit="IN"><x> {x_noz:.2f} </x><y>0</y><z>0</z></location>
    <orient unit="DEG"><pitch>0</pitch><roll>0</roll><yaw>0</yaw></orient>
   </thruster>
  </engine>
  <tank type="FUEL" number="0">
   <location unit="IN"><x> {x_motor:.2f} </x><y>0</y><z>0</z></location>
   <capacity unit="LBS"> {prop_left_lb:.15g} </capacity>
   <contents unit="LBS"> {prop_left_lb:.15g} </contents>
  </tank>
 </propulsion>
 <flight_control name="FCS">
  <property value="0.0">fcs/drogue-cmd</property>
  <property value="0.0">fcs/main-cmd</property>
  <channel name="Wind">
   <!-- ECMWF reanalysis profile at the site (2023-06-24 09:00 MDT) -->
   <fcs_function name="fcs/wind-east-fps">
    <function><table><independentVar lookup="row">position/h-sl-ft</independentVar>
     <tableData>
{wind_e}
     </tableData></table></function>
    <output>atmosphere/wind-east-fps</output>
   </fcs_function>
   <fcs_function name="fcs/wind-north-fps">
    <function><table><independentVar lookup="row">position/h-sl-ft</independentVar>
     <tableData>
{wind_n}
     </tableData></table></function>
    <output>atmosphere/wind-north-fps</output>
   </fcs_function>
  </channel>
  <channel name="Recovery">
   <!-- drogue 24 in Cd 1.6, main 72 in Cd 2.2 (team data) -->
   <fcs_function name="fcs/chute-cda-sqft">
    <function><sum>
     <product><property>fcs/drogue-cmd</property>
      <difference><value>1</value><property>fcs/main-cmd</property></difference>
      <value>{1.6*math.pi*0.3048**2*M2FT**2:.4f}</value></product>
     <product><property>fcs/main-cmd</property><value>{2.2*math.pi*0.9144**2*M2FT**2:.4f}</value></product>
    </sum></function>
   </fcs_function>
   <fcs_function name="fcs/body-aero-gain">
    <function><difference><value>1</value><property>fcs/drogue-cmd</property></difference></function>
   </fcs_function>
  </channel>
 </flight_control>
 <aerodynamics>
  <axis name="DRAG">
   <function name="aero/force/Drag_mach">
    <description>Team's power-off Cd curve</description>
    <product><property>aero/qbar-psf</property><property>metrics/Sw-sqft</property>
     <value> {CD_SCALE} </value>
     <table><independentVar lookup="row">velocities/mach</independentVar><tableData>
{tbl(cd_rows)}
     </tableData></table></product>
   </function>
   <function name="aero/force/Drag_chute">
    <product><property>aero/qbar-psf</property><property>fcs/chute-cda-sqft</property></product>
   </function>
  </axis>
  <axis name="LIFT">
   <function name="aero/force/Lift_alpha">
    <description>CN_alpha(M) sin(alpha), Barrowman (prometheus_inputs.py)</description>
    <product><property>fcs/body-aero-gain</property><property>aero/qbar-psf</property><property>metrics/Sw-sqft</property>
     <sin><property>aero/alpha-rad</property></sin>
     <table><independentVar lookup="row">velocities/mach</independentVar><tableData>
{tbl(cna_rows)}
     </tableData></table></product>
   </function>
  </axis>
  <axis name="SIDE">
   <function name="aero/force/Side_beta">
    <product><property>fcs/body-aero-gain</property><property>aero/qbar-psf</property><property>metrics/Sw-sqft</property>
     <sin><property>aero/beta-rad</property></sin>
     <value>-1</value>
     <table><independentVar lookup="row">velocities/mach</independentVar><tableData>
{tbl(cna_rows)}
     </tableData></table></product>
   </function>
  </axis>
  <axis name="PITCH">
   <function name="aero/moment/Pitch_damp">
    <product><property>fcs/body-aero-gain</property><property>aero/qbar-psf</property><property>metrics/Sw-sqft</property>
     <property>metrics/cbarw-ft</property><property>aero/ci2vel</property>
     <property>velocities/q-aero-rad_sec</property><value> {cmq:.1f} </value></product>
   </function>
  </axis>
  <axis name="YAW">
   <function name="aero/moment/Yaw_damp">
    <product><property>fcs/body-aero-gain</property><property>aero/qbar-psf</property><property>metrics/Sw-sqft</property>
     <property>metrics/bw-ft</property><property>aero/bi2vel</property>
     <property>velocities/r-aero-rad_sec</property><value> {cmq:.1f} </value></product>
   </function>
  </axis>
  <axis name="ROLL">
   <function name="aero/moment/Roll_damp">
    <product><property>aero/qbar-psf</property><property>metrics/Sw-sqft</property>
     <property>metrics/bw-ft</property><property>aero/bi2vel</property>
     <property>velocities/p-aero-rad_sec</property><value> -0.8 </value></product>
   </function>
  </axis>
 </aerodynamics>
 <output name="data_output/prometheus.csv" type="CSV" rate="50">
  <property caption="h_agl_ft"> position/h-agl-ft </property>
  <property caption="h_sl_ft"> position/h-sl-ft </property>
  <property caption="vt_fps"> velocities/vt-fps </property>
  <property caption="hdot_fps"> velocities/h-dot-fps </property>
  <property caption="mach"> velocities/mach </property>
  <property caption="theta_deg"> attitude/theta-deg </property>
  <property caption="alpha_deg"> aero/alpha-deg </property>
  <property caption="thrust_lbf"> propulsion/engine[0]/thrust-lbs </property>
  <property caption="mass_lb"> inertia/weight-lbs </property>
  <property caption="lat_deg"> position/lat-geod-deg </property>
  <property caption="lon_deg"> position/long-gc-deg </property>
 </output>
</fdm_config>
'''
open(f"{ROOT}/aircraft/prometheus/prometheus.xml", "w").write(ac)

h0_m = 1401 + RAIL * math.sin(math.radians(INCL))
v_ft = v_exit * M2FT
open(f"{ROOT}/aircraft/prometheus/rail_exit.xml", "w").write(f'''<?xml version="1.0"?>
<initialize name="rail_exit">
  <!-- state at rail exit from the 1-D rail integration; terrain at 1401 m -->
  <ubody unit="FT/SEC"> {v_ft:.3f} </ubody>
  <vbody unit="FT/SEC"> 0.0 </vbody>
  <wbody unit="FT/SEC"> 0.0 </wbody>
  <phi unit="DEG"> 0.0 </phi>
  <theta unit="DEG"> {INCL} </theta>
  <psi unit="DEG"> {HEAD} </psi>
  <latitude unit="DEG" type="geod"> 32.939377 </latitude>
  <longitude unit="DEG"> -106.911986 </longitude>
  <elevation unit="FT"> {1401*M2FT:.2f} </elevation>
  <altitudeAGL unit="FT"> {(h0_m-1401)*M2FT:.2f} </altitudeAGL>
</initialize>
''')
open(f"{ROOT}/scripts/prometheus.xml", "w").write(f'''<?xml version="1.0"?>
<runscript name="prometheus validation">
  <use aircraft="prometheus" initialize="rail_exit"/>
  <run start="0.0" end="600" dt="{os.environ.get("SIM_DT", "0.001")}">
    <event name="Ignition">
      <condition> simulation/sim-time-sec ge 0.0 </condition>
      <set name="atmosphere/delta-T" value="{1.8 * dT:.3f}"/>  <!-- JSBSim delta-T is in Rankine -->
      <set name="propulsion/engine[0]/set-running" value="1"/>
      <set name="fcs/throttle-cmd-norm" value="1.0"/>
    </event>
    <event name="Apogee">
      <condition> simulation/sim-time-sec gt 5.0
                  velocities/h-dot-fps lt 0.0 </condition>
      <set name="fcs/drogue-cmd" value="1"/>
      <notify><property>position/h-agl-ft</property></notify>
    </event>
    <event name="Main">
      <condition> fcs/drogue-cmd ge 1.0
                  position/h-agl-ft lt 1500.0 </condition>
      <set name="fcs/main-cmd" value="1"/>
    </event>
    <event name="Ground">
      <condition> simulation/sim-time-sec gt 5.0
                  position/h-agl-ft lt 3.0 </condition>
      <set name="simulation/terminate" value="1"/>
    </event>
  </run>
</runscript>
''')
json.dump(dict(t_exit=t_exit, v_exit=v_exit, dT=dT, isp=isp, impulse=I_tot, cmq=cmq,
               motor=motor_info, input_file=input_file, weather=o.get("weather_metadata", {"time_utc":o.get("weather_time_utc"),"role":"calibration proxy"}),
               cg_dry_in=cg_dry, x_cp_in=x_cp, static_margin_cal=(o["com_dry"] - cp_sub) / D),
          open(INFO, "w"), indent=1)
print(json.load(open(INFO)))
