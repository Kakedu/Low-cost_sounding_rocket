"""Reproduce both external flight comparisons with fail-fast paired manifests.

Prometheus uses independent 2022 NOAA regional soundings, not its apogee-fitted
temperature. Both cases include amplitude and direction weather envelopes.
All uncertainty distributions are conditional engineering assumptions.
"""
from pathlib import Path
import json,math,random,subprocess,sys,os,tempfile,concurrent.futures,time,hashlib
import numpy as np
from validation_atmosphere import baro_reading
HERE=Path(__file__).resolve().parent

def fly(name,variables,persist=None):
    with tempfile.TemporaryDirectory(prefix=f'validation_{name}_') as td:
        root=Path(td);info=root/'info.json';out=root/'result.json'
        env=dict(os.environ,JSB_ROOT=str(root/'jsb'),INFO=str(info),OUT=str(out))
        env.update({k:str(v) for k,v in variables.items()})
        logs=[]
        for driver in ('build_'+name+'.py','run_'+name+'.py'):
            result=subprocess.run([sys.executable,str(HERE/driver)],env=env,capture_output=True,text=True)
            logs.append(dict(driver=driver,exit=result.returncode,stdout=result.stdout,stderr=result.stderr))
            if result.returncode:
                raise RuntimeError(json.dumps(dict(name=name,variables=variables,logs=logs)))
        res=json.loads(out.read_text());build=json.loads(info.read_text())
        if persist:
            persist=Path(persist);persist.mkdir(parents=True,exist_ok=True)
            (persist/'result.json').write_text(json.dumps(res),encoding='utf8')
            (persist/'build_info.json').write_text(json.dumps(build,indent=2),encoding='utf8')
            (persist/'execution.json').write_text(json.dumps(dict(variables=variables,logs=logs),indent=2),encoding='utf8')
        return res,build

def summarize(v):
    a=np.asarray(v,dtype=float)
    return dict(mean=float(a.mean()),sd=float(a.std(ddof=1)),p2_5=float(np.percentile(a,2.5)),
                p5=float(np.percentile(a,5)),p50=float(np.percentile(a,50)),p95=float(np.percentile(a,95)),
                p97_5=float(np.percentile(a,97.5)))

def monte_carlo(name,n=200,workers=3):
    seed=7 if name=='prometheus' else 11;rng=random.Random(seed);sl=math.sqrt(math.log(1+.2**2))
    inputs=[]
    for i in range(n):
        cd=math.exp(rng.gauss(-.5*sl*sl,sl));imp=rng.gauss(1,.025)
        if name=='prometheus':imp*=1+rng.uniform(0,.03)
        delta=rng.gauss(0,1.5);ws=rng.uniform(0,1.5) if name=='prometheus' else rng.uniform(.5,1.5)
        rot=rng.uniform(-30,30)
        variables=dict(CD_SCALE=cd,IMP_SCALE=imp,DT_ABS='reanalysis',DT_EXTRA=delta,WIND_SCALE=ws,WIND_ROT_DEG=rot)
        if name=='prometheus':variables['WEATHER_WEIGHT']=rng.uniform(0,1)
        inputs.append(dict(sample=i,cd=cd,imp=imp,delta_T_extra_K=delta,wind_scale=ws,wind_rotation_deg=rot,environment=variables))
    prefix='prom' if name=='prometheus' else 'cavour'
    start=time.time();out=[]
    (HERE/(prefix+'_mc_inputs.json')).write_text(json.dumps(dict(seed=seed,n=n,inputs=inputs),indent=2),encoding='utf8')
    def one(item):
        result,build=fly(name,item['environment'])
        return dict(item,apogee_m=result['apogee'],baro_m=baro_reading(result['apogee'],result['dT'],1401 if name=='prometheus' else 150),
                    dT=result['dT'],burnout_s=result['burnout_s'],terminated=result['terminal']['terminated'],
                    terminal_h_agl_ft=result['terminal']['h_agl_ft'])
    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as pool:
        futures={pool.submit(one,item):item for item in inputs}
        for fut in concurrent.futures.as_completed(futures):
            item=futures[fut]
            try:out.append(fut.result())
            except Exception as ex:
                (HERE/(prefix+'_mc_failure.json')).write_text(json.dumps(dict(input=item,error=str(ex)),indent=2),encoding='utf8')
                for other in futures:other.cancel()
                raise
            if len(out)%20==0:print(name,len(out),'/',n,'elapsed',round(time.time()-start,1),flush=True)
    out.sort(key=lambda x:x['sample'])
    assert len(out)==n and [x['sample'] for x in out]==list(range(n))
    common=dict(n=n,seed=seed,failures=0,completed=n,quantile='numpy linear',samples=out,
                assumptions={'cd_lognormal_cv':.2,'isp_sd_fraction':.025,'additional_column_fit_sd_K':1.5,
                             'wind_rotation_uniform_deg':[-30,30],'weather':'independent 2022 regional bracketing soundings' if name=='prometheus' else 'launch-hour 2023 ERA5',
                             'scope':'conditional input uncertainty; weather envelope does not quantify every spatial or profile-shape error'})
    if name=='prometheus':common.update(summarize([x['apogee_m'] for x in out]))
    else:
        common['apogee_m']=summarize([x['apogee_m'] for x in out]);common['baro_m']=summarize([x['baro_m'] for x in out])
    (HERE/(prefix+'_mc.json')).write_text(json.dumps(common,indent=2),encoding='utf8')
    return common

def main():
    resultdir=HERE.parent/'results_revised/validation';resultdir.mkdir(parents=True,exist_ok=True)
    checks={}
    for name in ('prometheus','cavour'):
        for dt in (.001,.0005,.00025):
            result,build=fly(name,{'SIM_DT':dt},resultdir/f'{name}_dt_{dt:g}')
            checks[f'{name}_{dt}']=dict(apogee_m=result['apogee'],t_apogee_s=result['t_apogee'],burnout_s=result['burnout_s'],
                impulse_Ns_after_rail=result['impulse_after_rail_Ns'],dt_s=dt,delta_T_K=result['dT'])
            print('convergence',name,dt,checks[f'{name}_{dt}'],flush=True)
            if dt==.001:
                (HERE/('jsb_result.json' if name=='prometheus' else 'jsb_result_cavour.json')).write_text(json.dumps(result),encoding='utf8')
                (HERE/('build_info.json' if name=='prometheus' else 'build_info_cavour.json')).write_text(json.dumps(build,indent=2),encoding='utf8')
        result,build=fly(name,{'IMP_SCALE':1.1},resultdir/f'{name}_isp_1p1')
        source=build['motor']['source_impulse_Ns']
        prop_used=build['motor']['burned_kg']
        before=source*1.1*(prop_used/build['motor']['prop_kg'])
        integrated=result['impulse_after_rail_Ns']+before
        checks[name+'_isp_scale_test']=dict(expected_Ns=source*1.1,integrated_Ns=integrated,
                                            rel_error=(integrated/(source*1.1)-1),burnout_s=result['burnout_s'])
        assert abs(integrated/(source*1.1)-1)<.001
    from flight_temperature import infer_from_logs
    flight_temperature=infer_from_logs()
    result,build=fly('prometheus',{'PROM_INPUTS':'prometheus_inputs_2022.json','DT_ABS':flight_temperature['column_dT_K']},resultdir/'prometheus_2022_flight_calibrated')
    checks['prometheus_calibration_proxy']=dict(apogee_m=result['apogee'],dT_K=result['dT'],
                                               role='apogee/pressure-informed calibration with the same 2022 regional wind',flight_temperature=flight_temperature)
    (resultdir/'convergence_and_scaling.json').write_text(json.dumps(checks,indent=2),encoding='utf8')
    for name in ('prometheus','cavour'):monte_carlo(name,200,3)
    u_num=max(abs(checks['prometheus_0.001']['apogee_m']-checks['prometheus_0.00025']['apogee_m']),
              abs(checks['cavour_0.001']['apogee_m']-checks['cavour_0.00025']['apogee_m']),1.)
    env=dict(os.environ,U_NUM_M=str(u_num))
    for driver in ('gps_check.py','cavour_check.py'):
        run=subprocess.run([sys.executable,str(HERE/driver)],env=env,capture_output=True,text=True)
        (resultdir/(driver+'.log')).write_text(run.stdout+'\n'+run.stderr,encoding='utf8')
        if run.returncode:raise RuntimeError(run.stderr)
    for name in ('prom_mc.json','cavour_mc.json','gps_check.json','cavour_check.json','build_info.json','build_info_cavour.json'):
        (resultdir/name).write_bytes((HERE/name).read_bytes())
    for name in ('fig_validation_gps.png','fig_validation_cavour.png'):
        (HERE.parent/'results_revised'/name).write_bytes((HERE/name).read_bytes())
    print('VALIDATION SUITE COMPLETED; numerical uncertainty scale',u_num,flush=True)

if __name__=='__main__':main()
