"""Chunked paired SIL runs. A fresh Python process loads each frozen source job.

Chunking changes orchestration only: each sample consumes the original seed1
vehicle draw sequence and its original independent per-flight sensor RNG seed.
Completed chunks are validated before resuming. Failures are never dropped.
"""
import argparse,csv,hashlib,json,os
from pathlib import Path
import random,shutil,subprocess,sys,time,traceback
from concurrent.futures import ThreadPoolExecutor,as_completed
HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
SOURCE=HERE/'repair_source_v3'
OUT=HERE/'repaired_runs_v3'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def loaded(home,expected):
    result={}
    for name in ['sil','sweep','kf_estimator','kf.attitude_sensor','atmosphere_model','motor_model']:
        module=sys.modules.get(name)
        if module is None:continue
        path=Path(module.__file__).resolve();rel=path.relative_to(home).as_posix();checksum=sha(path)
        if checksum!=expected[rel]:raise RuntimeError('Loaded source hash mismatch: '+name)
        result[name]=dict(file=str(path),relative_path=rel,sha256=checksum)
    return result
def worker(job,start,stop):
    destination=OUT/'chunks'/job['name'].removesuffix('.csv')/f'{start:04}_{stop:04}'
    home=destination/'source';shutil.copytree(SOURCE,home);(home/'data_output').mkdir(exist_ok=True)
    os.chdir(home);sys.path.insert(0,str(home))
    import sil,sweep
    for key,value in sweep._reroot(home).items():setattr(sweep,key,value)
    sil.TARGET_FT=job.get('target',10000.);sil.KF_NAIVE=not job.get('qaware',False)
    sil.USE_PAD_TEMP=not job.get('isa',False);sil.PAD_TEMP_SD_K=job.get('pad_temp_sd',3.)
    sil.CTRL_DELAY_S=0.
    expected=json.loads((SOURCE/'source_manifest.json').read_text())['source_hashes']
    rec=dict(job=job,start_sample=start,stop_sample=stop,status='RUNNING',completed=0,
        source_snapshot=home.relative_to(ROOT).as_posix(),source_hashes=expected,loaded_sources=loaded(home,expected),
        runner_sha256=sha(__file__),started_unix=time.time())
    metadata=destination/'run.json';metadata.write_text(json.dumps(rec,indent=2))
    rng=random.Random(job['seed']);sig=dict(sweep.MC_SIGMA)
    if 'cd_sigma' in job:sig['cd0']=job['cd_sigma']
    rows=[]
    try:
        for sample in range(stop):
            factors,wind=sweep.draw_inputs(rng,sig)
            if sample<start:continue
            if job.get('fix_brake'):factors['airbrake_dcd']=1.
            noise_seed=job['seed']*1000003+sample
            metrics=sil.run_case(factors,wind,job['est'],job['cd'],random.Random(noise_seed),
                ideal=job.get('ideal',False),errors=job.get('errors',sil.ALL_ERRORS))
            rows.append(dict(sample=sample,**{k:factors[k] for k in sorted(factors)},wind_fps=wind[0],wind_from=wind[1],**metrics))
            rec.update(completed=len(rows),elapsed=time.time()-rec['started_unix'],loaded_sources=loaded(home,expected))
            metadata.write_text(json.dumps(rec,indent=2))
            if len(rows)%25==0:print(job['name'],sample,'completed',flush=True)
        output=destination/'metrics.csv'
        with output.open('w',newline='') as f:
            w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
        rec.update(status='SUCCESS',output=output.relative_to(ROOT).as_posix(),output_sha256=sha(output),
            elapsed=time.time()-rec['started_unix'])
    except BaseException as error:
        rec.update(status='FAILED',exception=repr(error),traceback=traceback.format_exc(),failed_sample=sample,
            factors=factors,wind=wind,sensor_rng_seed=job['seed']*1000003+sample)
        metadata.write_text(json.dumps(rec,indent=2));raise
    metadata.write_text(json.dumps(rec,indent=2))
def dispatch(job,start,stop,resume):
    folder=OUT/'chunks'/job['name'].removesuffix('.csv')/f'{start:04}_{stop:04}'
    if resume and (folder/'run.json').exists():
        record=json.loads((folder/'run.json').read_text());output=ROOT/record.get('output','missing')
        if record['status']=='SUCCESS' and record['job']==job and output.is_file() and sha(output)==record['output_sha256']:
            expected=json.loads((SOURCE/'source_manifest.json').read_text())['source_hashes']
            if record['source_hashes']!=expected:raise RuntimeError('Resume source mismatch')
            return record
        raise RuntimeError('Existing incomplete/failed chunk requires a new recorded attempt: '+str(folder))
    folder.mkdir(parents=True,exist_ok=False)
    args=[sys.executable,str(Path(__file__).resolve()),'--worker','--job',json.dumps(job),'--start',str(start),'--stop',str(stop)]
    with (folder/'run.log').open('w') as f:p=subprocess.run(args,stdout=f,stderr=subprocess.STDOUT)
    record=json.loads((folder/'run.json').read_text()) if (folder/'run.json').exists() else dict(job=job,status='FAILED',start_sample=start,stop_sample=stop,returncode=p.returncode)
    if p.returncode!=0 and record.get('status')=='SUCCESS':raise RuntimeError('Nonzero subprocess contradicts SUCCESS')
    return record
def merge(job,records):
    ordered=sorted(records,key=lambda r:r['start_sample'])
    if any(r['status']!='SUCCESS' for r in ordered):return dict(job=job,status='FAILED',chunks=ordered)
    rows=[]
    for r in ordered:rows.extend(csv.DictReader((ROOT/r['output']).open()))
    if [int(r['sample']) for r in rows]!=list(range(job['n'])):raise RuntimeError('Incomplete/duplicate sample IDs')
    output=OUT/job['name']
    with output.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    return dict(job=job,status='SUCCESS',completed=len(rows),sha256=sha(output),output=output.relative_to(ROOT).as_posix(),chunks=ordered)
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--worker',action='store_true');ap.add_argument('--job');ap.add_argument('--start',type=int);ap.add_argument('--stop',type=int)
    ap.add_argument('--workers',type=int,default=8);ap.add_argument('--chunk-size',type=int,default=100)
    ap.add_argument('--quality',action='store_true');ap.add_argument('--resume',action='store_true');ap.add_argument('--manifest-name',default='repair_run_manifest.json');a=ap.parse_args()
    if a.worker:worker(json.loads(a.job),a.start,a.stop);return
    protocol=json.loads((HERE/'repair_protocol_v3.json').read_text());jobs=protocol['jobs']
    if a.quality:jobs=[j for j in jobs if j['name'] in protocol['quality_jobs']]
    OUT.mkdir(parents=True,exist_ok=True)
    manifest=dict(status='RUNNING',protocol_sha256=sha(HERE/'repair_protocol_v3.json'),source_manifest_sha256=sha(SOURCE/'source_manifest.json'),runner_sha256=sha(__file__),jobs=jobs,started_unix=time.time(),command=sys.argv)
    mf=OUT/a.manifest_name;mf.write_text(json.dumps(manifest,indent=2));records=[]
    with ThreadPoolExecutor(max_workers=a.workers) as pool:
        futures=[pool.submit(dispatch,j,start,min(start+a.chunk_size,j['n']),a.resume) for j in jobs for start in range(0,j['n'],a.chunk_size)]
        for f in as_completed(futures):
            r=f.result();records.append(r);print(r['job']['name'],r['start_sample'],r['stop_sample'],r['status'],flush=True)
    groups=[merge(j,[r for r in records if r['job']==j]) for j in jobs]
    manifest.update(status='SUCCESS' if all(r['status']=='SUCCESS' for r in groups) else 'FAILED',records=groups,elapsed=time.time()-manifest['started_unix'])
    mf.write_text(json.dumps(manifest,indent=2));print(manifest['status'],len(groups),'studies',sum(r.get('completed',0) for r in groups),'flights')
    sys.exit(0 if manifest['status']=='SUCCESS' else 1)
if __name__=='__main__':main()
