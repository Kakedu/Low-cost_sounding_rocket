"""Check compiled JSBSim alpha/beta forces, normalized by the live qS.

The portable fixture copies inputs into a private sandbox and releases the
launch rail there. It never modifies the source model. Run verify_aero.py
first to produce the optional matching XML-hash evidence.
"""
from pathlib import Path
import argparse,hashlib,json,math,shutil,xml.etree.ElementTree as ET
import jsbsim

KINDS=('Lift_alpha','Side_beta','Lift_crossflow','Side_crossflow','Drag_crossflow','Drag_beta')
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def interpolate(rows,x):
    if x<=rows[0][0]:return rows[0][1]
    if x>=rows[-1][0]:return rows[-1][1]
    for (x0,y0),(x1,y1) in zip(rows,rows[1:]):
        if x0<=x<=x1:return y0+(y1-y0)*(x-x0)/(x1-x0)
    raise RuntimeError('Mach lookup interval not found')

def main():
    root=Path(__file__).resolve().parent
    p=argparse.ArgumentParser()
    p.add_argument('--scratch',type=Path,default=root/'.reproduction/compiled_sandbox')
    p.add_argument('--json',type=Path,default=root/'.reproduction/compiled_aero_runtime.json')
    p.add_argument('--xml-manifest',type=Path,default=root/'results_revised/aero_check.json')
    p.add_argument('--tolerance',type=float,default=1e-10)
    args=p.parse_args();scratch=args.scratch.resolve()
    if scratch==root or any(scratch==root/f for f in ('aircraft','engine','scripts')):
        raise ValueError('Scratch must be separate from source input directories')
    if args.tolerance<=0:raise ValueError('Positive comparison tolerance required')
    xml_hashes={ac:sha(root/'aircraft'/ac/(ac+'.xml')) for ac in ('myrocket','myrocket2','myrocket_cluster','myrocket_series')}
    if args.xml_manifest.exists():
        expected={r['aircraft']:r['sha256'] for r in json.loads(args.xml_manifest.read_text(encoding='utf-8'))}
        if xml_hashes!=expected:raise RuntimeError('Current XML hashes differ from structural-force verification evidence')
    for d in ('aircraft','engine'):shutil.copytree(root/d,scratch/d,dirs_exist_ok=True)
    (scratch/'data_output').mkdir(exist_ok=True)
    model=ET.parse(root/'aircraft/myrocket/myrocket.xml').getroot()
    fn=next(x for x in model.findall('.//function') if x.get('name')=='aero/force/Lift_alpha')
    mach_rows=[tuple(map(float,line.split())) for line in fn.find('.//tableData').text.strip().splitlines()]
    coefficient=json.loads((root/'barrowman_v6.json').read_text())['crossflow']['single']['coef']
    values=[];errors=[]
    for angle in ('alpha','beta'):
        for deg in (-15,-1,1,15):
            f=jsbsim.FGFDMExec(str(scratch));f.set_debug_level(0)
            if not f.load_model('myrocket'):raise RuntimeError('Fixture model load failed')
            f.set_dt(.001)
            if not f.load_ic('reset',True):raise RuntimeError('Fixture IC load failed')
            rad=math.radians(deg)
            f['ic/h-agl-ft']=1000;f['ic/u-fps']=300*math.cos(rad)
            f['ic/v-fps']=300*math.sin(rad) if angle=='beta' else 0
            f['ic/w-fps']=300*math.sin(rad) if angle=='alpha' else 0
            if not f.run_ic():raise RuntimeError('Fixture run_ic failed')
            f.disable_output();f['fcs/target-apogee-ft']=1e9;f['fcs/rail-travel-ft']=-1
            for _ in range(3):
                if not f.run():raise RuntimeError('Compiled-force fixture run failed')
            a=f['aero/alpha-rad'];b=f['aero/beta-rad'];q=f['aero/qbar-psf'];area=f['metrics/Sw-sqft'];mach=f['velocities/mach'];gain=f['fcs/body-aero-gain']
            raw={k:f['aero/force/'+k] for k in KINDS}
            if not all(math.isfinite(v) for v in (a,b,q,area,mach,gain,*raw.values())) or q*area<=0:
                raise RuntimeError('Nonfinite force fixture or nonpositive qS')
            if gain!=1:raise RuntimeError('Fixture did not release rail/body aero suppression')
            cn=interpolate(mach_rows,mach);c=coefficient
            expected={'Lift_alpha':cn*.5*math.sin(2*a),'Side_beta':-cn*.5*math.sin(2*b),'Lift_crossflow':c*math.sin(a)*abs(math.sin(a))*math.cos(a),'Side_crossflow':-c*math.sin(b)*abs(math.sin(b))*math.cos(b),'Drag_crossflow':c*abs(math.sin(a))**3,'Drag_beta':c*abs(math.sin(b))**3}
            normalized={k:raw[k]/(q*area) for k in KINDS}
            err={k:abs(normalized[k]-expected[k]) for k in KINDS}
            if max(err.values())>args.tolerance:raise RuntimeError(f'Compiled-force mismatch {angle} {deg}: {err}')
            if raw['Drag_crossflow']<0 or raw['Drag_beta']<0:raise RuntimeError('Crossflow drag must be nonnegative')
            if angle=='alpha' and raw['Lift_alpha']*deg<=0:raise RuntimeError('Lift sign mismatch')
            if angle=='beta' and raw['Side_beta']*deg>=0:raise RuntimeError('Side-force sign mismatch')
            values.append(dict(angle=angle,input_deg=deg,alpha_rad=a,beta_rad=b,mach=mach,qS_lbf=q*area,raw_force_lbf=raw,normalized=normalized,expected=expected,errors=err))
            errors.extend(err.values())
    if len(errors)!=48:raise RuntimeError('Expected eight cases and six compiled force components')
    evidence=dict(status='PASS',cases=8,compiled_force_components=48,max_coefficient_error=max(errors),tolerance=args.tolerance,xml_sha256=xml_hashes,engine_sha256=sha(root/'engine/M1670.xml'),source_sha256=sha(Path(__file__)),xml_hash_evidence_checked=args.xml_manifest.exists(),fixture='1000 ft, 300 fps; alpha or beta ±1°/±15°; private rail release; three 1 ms steps',values=values)
    args.json.parent.mkdir(parents=True,exist_ok=True);args.json.write_text(json.dumps(evidence,indent=2),encoding='utf-8')
    print(f'PASS: {len(errors)} compiled-force components, max coefficient error {max(errors):.3g}')

if __name__=='__main__':main()
