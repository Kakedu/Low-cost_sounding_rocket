"""Mass-table and integration convergence against the input RASP time curve."""
from pathlib import Path
import argparse, csv, json, math, shutil, time, xml.etree.ElementTree as ET
import numpy as np
import jsbsim
from make_motor import convert_rasp, read_rasp, thrust_at, LBF


def check(root, scratch, sample_dt, sim_dt):
    dest=scratch/f's{sample_dt:g}_dt{sim_dt:g}'
    for folder in ('aircraft','engine'):
        shutil.copytree(root/folder,dest/folder,dirs_exist_ok=True)
    (dest/'data_output').mkdir(exist_ok=True)
    info=convert_rasp(root/'M1670.eng',dest/'engine/M1670.xml',sample_dt)
    ac=dest/'aircraft/myrocket/myrocket.xml'
    tree=ET.parse(ac);tanks=tree.getroot().findall('./propulsion/tank')
    if len(tanks)!=1: raise RuntimeError('Single-stage convergence model must have exactly one tank')
    for key in ('capacity','contents'):
        element=tanks[0].find(key)
        if element is None or element.get('unit')!='LBS': raise RuntimeError('Expected tank mass in LBS')
        element.text=format(info['prop_lb'],'.17g')
    tree.write(ac,encoding='utf-8',xml_declaration=True)
    f=jsbsim.FGFDMExec(str(dest));f.set_debug_level(0)
    if not f.load_model('myrocket'): raise RuntimeError('Model load failed')
    f.set_dt(sim_dt)
    if not f.load_ic('reset',True): raise RuntimeError('Initial conditions failed')
    if not f.run_ic(): raise RuntimeError('Initial condition initialization failed')
    f.disable_output()
    fuel_initial=f['propulsion/tank[0]/contents-lbs']
    if not math.isclose(fuel_initial,info['prop_lb'],abs_tol=1e-12,rel_tol=0):
        raise RuntimeError('Initial tank mass does not match converter mass')
    f['fcs/target-apogee-ft']=1e9
    f['propulsion/engine[0]/set-running']=1;f['fcs/throttle-cmd-norm']=1
    points=read_rasp(root/'M1670.eng')[1]
    trace=[];last=0.;impulse=0.;fuel0=info['prop_lb'];fuel_end=None
    start=time.perf_counter();apogee=0.
    while f['simulation/sim-time-sec']<40:
        if not f.run(): raise RuntimeError('Integration failed')
        t=f['simulation/sim-time-sec'];F=f['propulsion/engine[0]/thrust-lbs']*LBF
        if not all(math.isfinite(f[k]) for k in ('position/h-agl-ft','velocities/h-dot-fps','aero/qbar-psf','velocities/mach','propulsion/tank[0]/contents-lbs')) or not math.isfinite(F):
            raise RuntimeError('Non-finite state in motor verification')
        fuel=f['propulsion/tank[0]/contents-lbs'];impulse+=F*sim_dt
        if F>1: last=t
        if fuel_end is None and fuel<1e-10: fuel_end=t
        if t<4.2:
            # propulsion is evaluated at the preceding integration instant.
            raw=thrust_at(points,max(0,t-sim_dt))
            trace.append((t,F,raw,fuel))
        apogee=max(apogee,f['position/h-agl-ft'])
        if t>5 and f['velocities/h-dot-fps']<0: break
    else:
        raise RuntimeError('Motor verification did not reach apogee by 40 seconds')
    errors=[x[1]-x[2] for x in trace]
    out={**info,'sim_dt_s':sim_dt,'observed_impulse_Ns':impulse,'impulse_error_pct':100*(impulse/info['impulse_Ns']-1),
         'last_thrust_over_1N_s':last,'fuel_depletion_s':fuel_end,'remaining_prop_lb':fuel,'initial_prop_lb':fuel_initial,
         'thrust_rms_error_N':float(np.sqrt(np.mean(np.square(errors)))),
         'thrust_max_error_N':max(map(abs,errors)),'apogee_ft':apogee,'runtime_s':time.perf_counter()-start}
    with open(dest/'motor_trace.csv','w',newline='') as fh:
        writer=csv.writer(fh);writer.writerow(['time_s','jsb_thrust_N','rasp_thrust_at_step_start_N','prop_lb']);writer.writerows(trace)
    return out


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--root',default=str(Path(__file__).resolve().parent))
    parser.add_argument('--scratch',required=True);parser.add_argument('--json',required=True)
    parser.add_argument('--selected',action='store_true',help='Check only the standard 2 ms table / 1 ms simulation pair')
    args=parser.parse_args();root=Path(args.root);scratch=Path(args.scratch);results=[]
    pairs=((.002,.001),) if args.selected else ((.01,.001),(.005,.001),(.002,.001),(.001,.001),(.0005,.001),(.002,.0005),(.002,.00025))
    for sample_dt,sim_dt in pairs:
        r=check(root,scratch,sample_dt,sim_dt);results.append(r);print(json.dumps(r),flush=True)
        Path(args.json).write_text(json.dumps(results,indent=2),encoding='utf-8')
    selected=[r for r in results if r['sample_dt_s']==.002 and r['sim_dt_s']==.001][0]
    sampling_error=0 if args.selected else abs(selected['apogee_ft']-next(r for r in results if r['sample_dt_s']==.0005)['apogee_ft'])
    if abs(selected['impulse_error_pct'])>.05 or sampling_error>1 or selected['fuel_depletion_s'] is None or abs(selected['fuel_depletion_s']-3.9)>.01:
        raise RuntimeError('Default motor table failed acceptance thresholds')


if __name__=='__main__': main()
