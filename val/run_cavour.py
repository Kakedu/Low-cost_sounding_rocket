"""Run cavour; require the intended termination event and complete propulsion."""
import csv,json,os,math
from pathlib import Path
import jsbsim
HERE=Path(__file__).resolve().parent
ROOT=Path(os.environ.get("JSB_ROOT",str(HERE/"jsb_cavour")))
INFO=Path(os.environ.get("INFO",str(HERE/"build_info_cavour.json")))
OUT=Path(os.environ.get("OUT",str(HERE/"jsb_result_cavour.json")))
os.chdir(ROOT)
csvpath=ROOT/"data_output/cavour.csv"
if csvpath.exists():csvpath.unlink()
info=json.loads(INFO.read_text())
f=jsbsim.FGFDMExec(str(ROOT));f.set_debug_level(0)
if not f.load_script("scripts/cavour.xml"):raise RuntimeError("JSBSim load_script failed")
if not f.run_ic():raise RuntimeError("JSBSim run_ic failed")
impulse=0.;burnout=None;last_time=0.
while f.run():
    thrust=f['propulsion/engine[0]/thrust-lbs'];now=f.get_sim_time()
    impulse+=thrust*(now-last_time)*4.4482216152605;last_time=now
    if burnout is None and f['propulsion/tank[0]/contents-lbs']<1e-9:burnout=now+info['t_exit']
terminal={'t_s':f.get_sim_time()+info['t_exit'],'h_agl_ft':f['position/h-agl-ft'],
          'v_vertical_fps':f['velocities/h-dot-fps'],'terminated':f['simulation/terminate'],
          'fuel_lb':f['propulsion/tank[0]/contents-lbs']}
if terminal['terminated']!=1:raise RuntimeError('Run ended without expected event '+str(terminal))
if 'cavour'=='prometheus' and terminal['h_agl_ft']>=3.1:raise RuntimeError('Landing incomplete '+str(terminal))
if 'cavour'=='cavour' and terminal['v_vertical_fps']>=0:raise RuntimeError('Apogee event incomplete '+str(terminal))
if burnout is None:raise RuntimeError('Motor did not reach burnout')
del f
rows=list(csv.DictReader(csvpath.open()))
t=[float(r['Time'])+info['t_exit'] for r in rows];h=[float(r['h_agl_ft'])*.3048 for r in rows]
if not all(map(math.isfinite,t+h)) or len(h)<2:raise RuntimeError('Invalid flight output')
i=max(range(len(h)),key=h.__getitem__)
result=dict(t=t,h=h,apogee=h[i],t_apogee=t[i],dT=info['dT'],
            max_mach=max(float(r['mach']) for r in rows),terminal=terminal,burnout_s=burnout,
            impulse_after_rail_Ns=impulse,motor=info.get('motor'),weather=info.get('weather'),
            simulated_dt_s=float(os.environ.get('SIM_DT','0.001')))
OUT.write_text(json.dumps(result),encoding='utf8')
print('Completed cavour: apogee',h[i],'m; burnout',burnout,'s; terminal',terminal)
