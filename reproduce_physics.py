"""Reproduce RQ1, RQ3, point-mass checks, structural loads and FEA.

Each command runs in a separate snapshot: generated models and data_output
are isolated. Failure stops publication of that command's result. Seed 1
is shared by the paired 500-flight staging datasets; post-processing seed 3.
"""
import argparse
import concurrent.futures
import json
import hashlib
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--workers',type=int,default=5)
    parser.add_argument('--sandbox',required=True);args=parser.parse_args()
    root=Path(__file__).resolve().parent;out=root/'results_revised';out.mkdir(exist_ok=True)
    sand=Path(args.sandbox).resolve();sand.mkdir(parents=True,exist_ok=True)
    jobs=[]
    def add(name,argv): jobs.append((name,argv))
    for rail,suffix in [(12,'anglesweep_v6.csv'),(18,'anglesweep_rail18_v6.csv')]:
        add(suffix,['anglesweep.py','--root','.','--rail-travel',str(rail),'--csv',str(out/suffix)])
    add('stage_compare_v6.csv',['stage_compare.py','--csv',str(out/'stage_compare_v6.csv')])
    for tilt in (1,5):
        for mode in ('coast','wind','tipoff','mc'):
            name=f'staging_{mode}_tilt{tilt}_v6.csv'
            argv=['staging_seq.py',mode,'--root','.','--tilt',str(tilt),'--csv',str(out/name)]
            if mode=='mc': argv+=['--n','500','--coasts','0,5,10,14','--seed','1']
            add(name,argv)
    for coast in (0,5):
        for k in (1.,1.3,1.5):
            name=f'sweep_staging_c{coast}_k{k:.1f}_v6.csv'
            add(name,['sweep_staging.py','--coast',str(coast),'--k-stack',str(k),'--csv',str(out/name)])
    for mode in ('verify','authority','estimator','control'):
        add(f'crosscheck_{mode}.csv',['crosscheck.py',mode,'--csv',str(out/f'crosscheck_{mode}.csv')])
    add('loads_v6.json',['loads.py','--root','.','--json',str(out/'loads_v6.json')])
    records=[]
    def run(job):
        name,argv=job;dest=sand/name.replace('.','_');dest.mkdir(exist_ok=True)
        for directory in ('aircraft','engine','scripts'):
            shutil.copytree(root/directory,dest/directory,dirs_exist_ok=True)
        for p in root.glob('*.py'): shutil.copy2(p,dest/p.name)
        for name2 in ('M1670.eng','barrowman_v6.json','motor_build_info.json'): shutil.copy2(root/name2,dest/name2)
        (dest/'data_output').mkdir(exist_ok=True)
        command=[sys.executable,*argv];start=time.perf_counter()
        with open(dest/'run.log','w',encoding='utf-8') as log:
            process=subprocess.run(command,cwd=dest,stdout=log,stderr=subprocess.STDOUT,env={**os.environ,'JSBSIM_DEBUG':'0'})
        record=dict(result=name,command=command,cwd=str(dest),runtime_s=time.perf_counter()-start,
                    seed=1 if 'staging_mc_' in name else None,exit_code=process.returncode,
                    failed_runs=0 if process.returncode==0 else 1,log=str(dest/'run.log'))
        inputs=[dest/p for p in ('engine/M1670.xml','aircraft/myrocket/myrocket.xml','aircraft/myrocket2/myrocket2.xml','aircraft/myrocket_cluster/myrocket_cluster.xml','aircraft/myrocket_series/myrocket_series.xml',argv[0])]
        record['inputs_sha256']={str(p.relative_to(dest)):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs}
        output=out/name
        if output.exists():record['output_sha256']=hashlib.sha256(output.read_bytes()).hexdigest()
        print(json.dumps(record),flush=True)
        return record
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as pool:
        for r in pool.map(run,jobs):
            records.append(r);(out/'physics_run_manifest.json').write_text(json.dumps(records,indent=2))
    if any(r['exit_code'] for r in records): raise RuntimeError('Physics run failed; inspect run manifests')
    # All dependencies are now present. Reuse independent snapshots for post-processing.
    followups=[]
    for tilt in (1,5):
        name=f'staging_post_tilt{tilt}_v6.txt';dest=sand/f'post_tilt{tilt}';dest.mkdir(exist_ok=True)
        command=[sys.executable,str(root/'staging_post.py'),str(out/f'staging_mc_tilt{tilt}_v6.csv'),'--seed','3','--draws','200']
        start=time.perf_counter()
        with open(out/name,'w') as fh: subprocess.run(command,cwd=dest,stdout=fh,stderr=subprocess.STDOUT,check=True)
        followups.append(dict(result=name,command=command,runtime_s=time.perf_counter()-start,seed=3,failed_runs=0))
    for name,argv in [('flap_fea_v6.json',['flap_fea.py','--loads',str(out/'loads_v6.json'),'--json',str(out/'flap_fea_v6.json')]),
                      ('crosscheck_staging.json',['crosscheck_staging.py','--jsb-csv',str(out/'stage_compare_v6.csv'),'--json',str(out/'crosscheck_staging.json')])]:
        followups.append(run((name,argv)))
    for src,fig in [('anglesweep_v6.csv','fig_rq1_angle_wind.png'),('anglesweep_rail18_v6.csv','fig_rq1_angle_wind_rail18.png')]:
        subprocess.run([sys.executable,str(root/'plot_rq1.py'),str(out),src,fig],cwd=root,check=True)
    (out/'physics_run_manifest.json').write_text(json.dumps(records+followups,indent=2))


if __name__=='__main__': main()
