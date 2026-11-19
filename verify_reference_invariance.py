"""Verify compiled body moments when the same wind force is moved to its CP.

This builds private, isolated force-family fixtures from the current XML.
One fixture retains the current reference-point couple; the other removes
that couple and moves AERORP by the same structural x offset. Kinematics
are held with JSBSim's zero-dt integration suspension. No flight output or
source model is edited. The comparison tests a coordinate identity, not
the empirical validity of the aerodynamic force law at large incidence.
"""
from pathlib import Path
import argparse,hashlib,json,math,shutil,xml.etree.ElementTree as ET
import jsbsim

FAMILIES=(
 ('stack_normal','myrocket2',0,2.4738,('Lift_alpha','Side_beta'),'Pitch_stack_cp','Yaw_stack_cp'),
 ('stack_crossflow','myrocket2',0,4.0356,('Lift_crossflow','Side_crossflow','Drag_crossflow','Drag_beta'),'Pitch_crossflow_cp','Yaw_crossflow_cp'),
 ('sust_crossflow','myrocket2',1,-.0675,('Lift_crossflow_sust','Side_crossflow_sust','Drag_crossflow_sust','Drag_beta_sust'),'Pitch_crossflow_cp_sust','Yaw_crossflow_cp_sust'),
 ('single_crossflow','myrocket',0,-.0675,('Lift_crossflow','Side_crossflow','Drag_crossflow','Drag_beta'),'Pitch_crossflow_cp','Yaw_crossflow_cp'),
 ('cluster_crossflow','myrocket_cluster',0,-.7196,('Lift_crossflow','Side_crossflow','Drag_crossflow','Drag_beta'),'Pitch_crossflow_cp','Yaw_crossflow_cp'),
 ('series_crossflow','myrocket_series',0,-.7196,('Lift_crossflow','Side_crossflow','Drag_crossflow','Drag_beta'),'Pitch_crossflow_cp','Yaw_crossflow_cp'))
FORCE_KEYS=('forces/fbx-aero-lbs','forces/fby-aero-lbs','forces/fbz-aero-lbs')
MOMENT_KEYS=('moments/l-aero-lbsft','moments/m-aero-lbsft','moments/n-aero-lbsft')
OFFICIAL_SOURCES={
 'version':'JSBSim v1.3.1',
 'wind_to_body':'https://raw.githubusercontent.com/JSBSim-Team/jsbsim/v1.3.1/src/models/FGAuxiliary.cpp',
 'force_and_moment_frames':'https://raw.githubusercontent.com/JSBSim-Team/jsbsim/v1.3.1/src/models/FGAerodynamics.cpp',
 'implementation_lines':'FGAuxiliary.cpp 296-328; FGAerodynamics.cpp 192-201, 246-267, 391-407, 468-470 (zero-based online raw lines)'}
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()

def build(root,scratch,family,mode):
    title,aircraft,phase,dx,kept,pitch,yaw=family
    directory=scratch/title/mode
    tree=ET.parse(root/'aircraft'/aircraft/(aircraft+'.xml'));doc=tree.getroot()
    rp=doc.find("metrics/location[@name='AERORP']/x")
    if mode=='reference':rp.text=f'{float(rp.text)+dx*12:.17g}'
    aero=doc.find('aerodynamics')
    original_pitch=next(f for f in aero.findall('.//function') if f.get('name')=='aero/moment/'+pitch)
    original_yaw=next(f for f in aero.findall('.//function') if f.get('name')=='aero/moment/'+yaw)
    if not math.isclose(float(original_pitch.find('product/value').text),dx,abs_tol=1e-12):
        raise RuntimeError(f'{title}: pitch is not the current body-axis couple convention')
    if not math.isclose(float(original_yaw.find('product/value').text),-dx,abs_tol=1e-12):
        raise RuntimeError(f'{title}: yaw offset mismatch')
    count=0
    for axis in aero.findall('axis'):
        for fn in list(axis):
            name=(fn.get('name') or '').split('/')[-1]
            retain=name in kept or (mode=='current' and name in (pitch,yaw))
            if not retain:axis.remove(fn);continue
            if name in kept:count+=1
            # The selected family is tested at its active gate; all other
            # forces and moments are absent. Source force expressions remain.
            for parent in fn.iter():
                for child in list(parent):
                    if child.tag=='property' and child.text in ('fcs/off-rail','fcs/body-aero-gain','fcs/staged'):
                        value=phase if child.text=='fcs/staged' else 1
                        index=list(parent).index(child);parent.remove(child)
                        new=ET.Element('value');new.text=str(value);parent.insert(index,new)
    if count!=len(kept):raise RuntimeError(f'{title}: missing force component')
    (directory/'aircraft'/aircraft).mkdir(parents=True,exist_ok=True)
    tree.write(directory/'aircraft'/aircraft/(aircraft+'.xml'),encoding='utf-8',xml_declaration=True)
    shutil.copy2(root/'aircraft'/aircraft/'reset.xml',directory/'aircraft'/aircraft/'reset.xml')
    shutil.copytree(root/'engine',directory/'engine',dirs_exist_ok=True)
    (directory/'data_output').mkdir(exist_ok=True)
    return directory

def run(directory,aircraft,alpha,beta):
    f=jsbsim.FGFDMExec(str(directory));f.set_debug_level(0)
    if not f.load_model(aircraft) or not f.load_ic('reset',True):raise RuntimeError('Fixture model/IC load failed')
    f['ic/h-agl-ft']=1000;f['ic/u-fps']=300*math.cos(alpha)*math.cos(beta)
    f['ic/v-fps']=300*math.sin(beta);f['ic/w-fps']=300*math.sin(alpha)*math.cos(beta)
    if not f.run_ic():raise RuntimeError('Fixture run_ic failed')
    f.disable_output();f['fcs/rail-travel-ft']=-1;f['fcs/target-apogee-ft']=1e9
    f.suspend_integration()
    for _ in range(3):
        if not f.run():raise RuntimeError('Compiled fixture run failed')
    keys=('aero/alpha-rad','aero/beta-rad','aero/qbar-psf','metrics/Sw-sqft','simulation/sim-time-sec',*FORCE_KEYS,*MOMENT_KEYS)
    record={k:f[k] for k in keys}
    if not all(math.isfinite(x) for x in record.values()):raise RuntimeError('Nonfinite fixture force/moment/state')
    if abs(record['simulation/sim-time-sec'])>1e-12:raise RuntimeError('Fixture integration was not suspended')
    if abs(record['aero/alpha-rad']-alpha)>1e-12 or abs(record['aero/beta-rad']-beta)>1e-12:
        raise RuntimeError('Fixture kinematics changed')
    if record['aero/qbar-psf']*record['metrics/Sw-sqft']<=0:raise RuntimeError('Nonpositive qS')
    return record

def main():
    root=Path(__file__).resolve().parent
    ap=argparse.ArgumentParser();ap.add_argument('--scratch',type=Path,default=root/'results_revised/verification/reference_sandbox')
    ap.add_argument('--json',type=Path,default=root/'results_revised/verification/reference_invariance_current.json')
    ap.add_argument('--tolerance',type=float,default=1e-9);args=ap.parse_args()
    scratch=args.scratch.resolve()
    if scratch==root or any(scratch==root/d for d in ('aircraft','engine','scripts')):raise ValueError('Scratch must be separate from source models')
    if args.tolerance<=0:raise ValueError('Positive tolerance required')
    values=[]
    for family in FAMILIES:
        fixtures={mode:build(root,scratch,family,mode) for mode in ('current','reference')}
        for deg in (-45,-15,-1,1,15,45):
            for axis,ad,bd in [('alpha',deg,0),('beta',0,deg),('combined',deg,deg)]:
                raw={mode:run(directory,family[1],math.radians(ad),math.radians(bd)) for mode,directory in fixtures.items()}
                errors={k:raw['current'][k]-raw['reference'][k] for k in (*FORCE_KEYS,*MOMENT_KEYS)}
                if max(abs(x) for x in errors.values())>args.tolerance:raise RuntimeError(f'{family[0]} {axis} {deg}: reference invariance failed: {errors}')
                values.append(dict(family=family[0],axis=axis,alpha_deg=ad,beta_deg=bd,dx_structural_ft=family[3],raw=raw,errors=errors))
    if len(values)!=108:raise RuntimeError('Expected 108 reference-invariance cases')
    evidence=dict(status='PASS',jsbsim_version=jsbsim.__version__,cases=len(values),compiled_models=2*len(values),max_force_error_lbf=max(abs(r['errors'][k]) for r in values for k in FORCE_KEYS),max_moment_error_lbf_ft=max(abs(r['errors'][k]) for r in values for k in MOMENT_KEYS),tolerance_lbf_or_lbf_ft=args.tolerance,
      fixture='300 fps,1000ft; zero-dt fixed kinematics; isolated active force family; compare current reference-point couple to force at shifted AERORP',
      limitation='Large-angle cases test coordinate consistency only, not empirical aerodynamic accuracy.',
      xml_sha256={ac:sha(root/'aircraft'/ac/(ac+'.xml')) for ac in ('myrocket','myrocket2','myrocket_cluster','myrocket_series')},engine_sha256=sha(root/'engine/M1670.xml'),source_sha256=sha(Path(__file__)),official_sources=OFFICIAL_SOURCES,values=values)
    args.json.parent.mkdir(parents=True,exist_ok=True);args.json.write_text(json.dumps(evidence,indent=2),encoding='utf-8')
    print(f'PASS: {len(values)} compiled reference-invariance cases, maximum moment difference {evidence["max_moment_error_lbf_ft"]:.3g} lbf ft')

if __name__=='__main__':main()
