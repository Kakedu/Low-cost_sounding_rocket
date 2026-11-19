"""Read frozen templates and run four isolated deterministic rail traces."""
from pathlib import Path
import sys,json,math,csv,shutil,hashlib
V=Path(__file__).resolve().parent;P=V
sys.path.insert(0,str(P))
import sweep,jsbsim
R=V/'.reproduction/rail_checks';R.mkdir(parents=True,exist_ok=True)
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
sources=[P/'aircraft/myrocket/myrocket.xml',P/'aircraft/myrocket/reset.xml',P/'engine/M1670.xml',P/'engine/direct.xml',P/'scripts/launch.xml',P/'sweep.py']
bound={str(p.relative_to(P)):sha(p) for p in sources}
out=[]
for angle,wind in [(89,0),(85,0),(89,29.333333),(85,29.333333)]:
 root=R/f'tilt{90-angle}_wind{wind:g}';root.mkdir(exist_ok=True)
 for name in ['aircraft','engine','scripts','data_output']:(root/name).mkdir(exist_ok=True)
 shutil.copyfile(P/'engine/direct.xml',root/'engine/direct.xml')
 sweep.JSBSIM_ROOT=root;sweep.CASE_AC=root/'aircraft/_sweep';sweep.CASE_ENG=root/'engine/_sweep_engine.xml';sweep.CASE_SCRIPT=root/'scripts/_sweep.xml';sweep.CASE_CSV=root/'data_output/_sweep.csv'
 sweep.build_case({'launch_elev':angle},open_loop=True,wind=(wind,270))
 f=jsbsim.FGFDMExec(str(root));f.set_debug_level(0)
 assert f.load_script('scripts/_sweep.xml')
 assert f.run_ic();f.disable_output()
 def axis():
  th=f['attitude/theta-rad'];ps=f['attitude/psi-rad']
  return (math.cos(th)*math.cos(ps),math.cos(th)*math.sin(ps),-math.sin(th))
 initial=axis();pos=[0.,0.,0.];maximum=0.;maxlat=0.;rows=[]
 while f['simulation/sim-time-sec']<2.01:
  assert f.run();t=f['simulation/sim-time-sec'];a=axis()
  departure=math.degrees(math.acos(max(-1,min(1,sum(x*y for x,y in zip(initial,a))))))
  for i,name in enumerate(['north','east','down']):pos[i]+=f[f'velocities/v-{name}-fps']*.001
  along=sum(x*y for x,y in zip(initial,pos));lateral=math.sqrt(max(0,sum(x*x for x in pos)-along*along))
  row=dict(t_s=t,off_rail=f['fcs/off-rail'],margin_ft=f['fcs/rail-margin-ft'],axis_departure_deg=departure,lateral_displacement_ft=lateral,theta_deg=f['attitude/theta-deg'],p_rad_s=f['velocities/p-rad_sec'],q_rad_s=f['velocities/q-rad_sec'],r_rad_s=f['velocities/r-rad_sec'])
  rows.append(row)
  if row['off_rail']>.5:break
  maximum=max(maximum,departure);maxlat=max(maxlat,lateral)
 with (root/'trace.csv').open('w',newline='',encoding='utf-8') as stream:
  writer=csv.DictWriter(stream,fieldnames=rows[0]);writer.writeheader();writer.writerows(rows)
 out.append(dict(launch_elevation_deg=angle,wind_fps=wind,wind_from_deg=270,release=rows[-1],max_guided_axis_departure_deg=maximum,max_guided_lateral_displacement_ft=maxlat,fallback_used=rows[-1]['margin_ft']<=0,trace_path=str((root/'trace.csv').relative_to(V)),trace_sha256=sha(root/'trace.csv')))
 assert all(sha(P/name)==h for name,h in bound.items())
(V/'.reproduction/rail_check_results.json').write_text(json.dumps(dict(status='PASS',source_sha256=bound,new_deterministic_runs=4,new_Monte_Carlo_flights=0,cases=out),indent=2),encoding='utf-8')
print(json.dumps(out,indent=2))
