"""Recompute control statistics from complete corrected studies.

Quantiles use NumPy's linear percentile definition throughout. Bootstrap
resamples flights or paired flight differences, never Gaussian substitutes.
"""
import argparse, csv, hashlib, json, math, pathlib, statistics as st
import numpy as np
from scipy.stats import binomtest

ROOT=pathlib.Path(__file__).resolve().parent
TARGET=10000.0

def rows(path,expected_n=None):
    with path.open(encoding='utf-8') as f: data=list(csv.DictReader(f))
    if not data or (expected_n is not None and len(data)!=expected_n):
        raise ValueError(f'Incomplete study {path}: expected {expected_n}, found {len(data)}')
    if 'sample' in (data[0] if data else {}):
        ids=[int(r['sample']) for r in data]
        if sorted(ids)!=list(range(len(data))): raise ValueError(f'Invalid sample IDs: {path}')
        data.sort(key=lambda r:int(r['sample']))
    return data

def wilson(k,n,z=1.96):
    p=k/n; d=1+z*z/n; c=(p+z*z/(2*n))/d
    h=z*math.sqrt(p*(1-p)/n+z*z/(4*n*n))/d
    return [c-h,c+h]

def apogees(data,metric='geometric'):
    a=np.array([float(r['apogee']) for r in data])
    if metric=='barometric':
        h=a*.3048; H=h*6356766/(6356766+h)
        dt=np.array([float(r.get('dT_K',0) or 0) for r in data])
        a=H/.3048*288.15/(288.15+dt)
    return a

def bootstrap_ci(values,statistic='median',n_boot=2000,seed=7):
    a=np.asarray(values); rng=np.random.default_rng(seed)
    sampled=a[rng.integers(0,len(a),size=(n_boot,len(a)))]
    fn={'median':np.median,'mean':np.mean,'sd':lambda x,axis:np.std(x,axis=axis,ddof=1)}[statistic]
    return np.percentile(fn(sampled,axis=1),[2.5,97.5],method='linear').tolist()

def stats(data,target=TARGET,metric='geometric',study_name=''):
    a=apogees(data,metric); err=a-target; n=len(a)
    if not np.isfinite(a).all(): raise ValueError(f'Nonfinite apogee in {study_name}')
    counts=dict(hit=int(np.sum(np.abs(err)<=100)),low=int(np.sum(err<-100)),high=int(np.sum(err>100)))
    result=dict(n=n,mean=float(a.mean()),sd=float(a.std(ddof=1)),min=float(a.min()),max=float(a.max()),
                quantiles=dict(zip(['p5','p50','p95'],np.percentile(a,[5,50,95],method='linear').tolist())),
                counts=counts,rates={k:v/n for k,v in counts.items()},wilson={k:wilson(v,n) for k,v in counts.items()},
                sd_bootstrap95=bootstrap_ci(a,'sd'),tolerances={str(t):int(np.sum(np.abs(err)<=t))/n for t in [50,75,100]})
    if 'cd_projected_updates' in data[0]:
        projections=np.array([int(r['cd_projected_updates']) for r in data])
        if (projections<0).any():raise ValueError('Negative projection counter')
        result['cd_projection']=dict(active_flights=int((projections>0).sum()),total_updates=int(projections.sum()),max_updates=int(projections.max()))
        if 'h_rms_m' in data[0]:
            heights=np.array([float(r['h_rms_m']) for r in data])
            groups={}
            for label,mask in [('no_projection',projections==0),('projection_active',projections>0)]:
                count=int(mask.sum());values=heights[mask]
                groups[label]=dict(n=count,hits=int((np.abs(err[mask])<=100).sum()),
                    h_rms_median=float(np.median(values)) if count else None,
                    h_rms_p95=float(np.percentile(values,95,method='linear')) if count else None,
                    h_rms_over100m=int((values>100).sum()))
            result['cd_projection'].update(descriptive_groups=groups,
                interpretation='Projection activity is an endogenous filter outcome; subgroup associations are descriptive, not a randomized treatment effect or proof of cause.')
    for key in ['h_rms_m','v_rms_ms','h_gate_rms_m','v_gate_rms_ms','cd_rms','cd_low_mach_rms','attitude_rms_deg','h_max_abs_m','v_max_abs_ms']:
        if key in data[0]:
            values=np.array([float(r[key]) for r in data]); valid=np.isfinite(values)
            for i in np.flatnonzero(~valid):
                no_gate=int(data[i].get('gate_sensor_samples',-1))==0
                no_cd_estimator='_truth_' in study_name and key in ['cd_rms','cd_low_mach_rms']
                gate_metric=key in ['h_gate_rms_m','v_gate_rms_ms','cd_rms','cd_low_mach_rms']
                if not (no_cd_estimator or (no_gate and gate_metric)):
                    raise ValueError(f'Nonfinite applicable metric {key}: {study_name}, sample {data[i].get("sample",i)}')
            values=values[valid]
            result[key]=dict(mean=float(values.mean()),median=float(np.median(values)),p95=float(np.percentile(values,95,method='linear')),p99=float(np.percentile(values,99,method='linear')),max=float(values.max()),
                             mean_bootstrap95=bootstrap_ci(values,'mean'),defined_n=int(valid.sum()),undefined_n=int((~valid).sum())) if len(values) else None
    if 'baro_rejections' in data[0]:
        result['observation_rejections']={}
        for label in ['baro','accel']:
            updates=np.array([int(r[label+'_updates']) for r in data]);rejected=np.array([int(r[label+'_rejections']) for r in data])
            if (updates<rejected).any() or (rejected<0).any():raise ValueError('Invalid observation rejection counter')
            fractions=np.divide(rejected,updates,out=np.zeros(n),where=updates>0)
            result['observation_rejections'][label]=dict(attempted=int(updates.sum()),rejected=int(rejected.sum()),pooled_fraction=float(rejected.sum()/max(1,updates.sum())),flights_with_rejection=int((rejected>0).sum()),per_flight_fraction_median=float(np.median(fractions)),per_flight_fraction_p95=float(np.percentile(fractions,95)),per_flight_fraction_max=float(fractions.max()))
    return result

def validate_pair(first,second,allow_brake_difference=False):
    if len(first) not in [300,500] or len(second) not in [300,500]:
        raise ValueError('Unplanned/incomplete paired study count')
    for data in [first,second]:
        if [int(r['sample']) for r in data]!=list(range(len(data))): raise ValueError('Missing paired sample ID')
    n=min(len(first),len(second))
    inputs=['isp','prop_mass','cd0','cn_alpha','cp_x','dry_mass','drogue_cda','main_cda','launch_elev','launch_az','wind_fps','dT_K']
    if not allow_brake_difference: inputs.append('airbrake_dcd')
    for i in range(n):
        for name in inputs:
            if name in first[i] and name in second[i] and not math.isclose(float(first[i][name]),float(second[i][name]),rel_tol=1e-11,abs_tol=1e-11):
                raise ValueError(f'Unpaired vehicle input {name}, sample {i}')
        def direction(row): return float(row.get('wind_from_deg',row.get('wind_from',0)))
        if not math.isclose(direction(first[i]),direction(second[i]),rel_tol=1e-11,abs_tol=1e-11): raise ValueError(f'Unpaired wind direction, sample {i}')
    return n

def paired(first,second,target=TARGET,metric='geometric',allow_brake_difference=False):
    n=validate_pair(first,second,allow_brake_difference); a=apogees(first[:n],metric); b=apogees(second[:n],metric)
    ha=np.abs(a-target)<=100; hb=np.abs(b-target)<=100
    only_first=int(np.sum(ha&~hb)); only_second=int(np.sum(hb&~ha)); d=np.abs(a-target)-np.abs(b-target)
    p=float(binomtest(only_first,only_first+only_second,.5).pvalue) if only_first+only_second else 1.
    return dict(n=n,comparison='planned first-300 prefix' if len(first)!=len(second) else 'same sample IDs',only_first=only_first,only_second=only_second,mcnemar_exact_p=p,
                median_abs_error_difference=float(np.median(d)),median_difference_bootstrap95=bootstrap_ci(d),
                mean_abs_error_difference=float(d.mean()),mean_difference_bootstrap95=bootstrap_ci(d,'mean'))

def budget(open_data,controlled,target=TARGET):
    n=validate_pair(open_data,controlled,allow_brake_difference=True); o=apogees(open_data[:n]); c=apogees(controlled[:n]); short=o<target-100
    reach=~short; hit=np.abs(c-target)<=100; residual=np.abs(c[reach]-target)
    return dict(n=n,energy_short=int(short.sum()),energy_short_and_missed=int(np.sum(short&~hit)),
                other_low=int(np.sum((c<target-100)&reach)),high=int(np.sum(c>target+100)),
                reachable=int(reach.sum()),reachable_hit=int(np.sum(reach&hit)),reachable_hit_rate=float(hit[reach].mean()),
                residual_p50_p95=np.percentile(residual,[50,95],method='linear').tolist())

def regression(data):
    names=['isp','cd0','cn_alpha','cp_x','dry_mass','launch_elev','wind_fps','airbrake_dcd','dT_K']
    if 'isp' not in data[0]: names[0]='prop_mass'
    X=np.array([[float(r[k]) for k in names] for r in data]); y=apogees(data)
    keep=X.std(0)>0; X=X[:,keep]; names=[k for k,take in zip(names,keep) if take]
    Z=(X-X.mean(0))/X.std(0); yz=(y-y.mean())/y.std(); A=np.column_stack([Z,np.ones(len(Z))])
    beta=np.linalg.lstsq(A,yz,rcond=None)[0]; r2=1-float(np.sum((yz-A@beta)**2)/np.sum(yz**2))
    return dict(r_squared=r2,standardized_coefficients={k:float(b) for k,b in zip(names,beta[:-1])})

def landing(data,n_boot=1000):
    pts=np.array([[float(r['land_e']),float(r['land_n'])] for r in data]); distance=np.linalg.norm(pts,axis=1)
    if not np.isfinite(pts).all(): raise ValueError('Nonfinite landing coordinate')
    def geometry(p):
        C=np.cov(p.T); w,V=np.linalg.eigh(C); centered=p-p.mean(0)
        md=np.sqrt(np.einsum('ij,jk,ik->i',centered,np.linalg.pinv(C),centered))
        radii=np.percentile(md,[68,95],method='linear')
        return w,md,radii
    w,md,radii=geometry(pts); rng=np.random.default_rng(11); boots=[]
    for _ in range(n_boot):
        wb,_,rb=geometry(pts[rng.integers(0,len(pts),len(pts))]); boots.append([np.sqrt(wb[-1])*r for r in rb])
    boots=np.array(boots)
    return dict(n=len(pts),center_east_north=pts.mean(0).tolist(),distance_quantiles=dict(zip(['p5','p50','p95','max'],np.r_[np.percentile(distance,[5,50,95],method='linear'),distance.max()].tolist())),
                ellipses={str(c):dict(radius=float(r),major_semiaxis=float(np.sqrt(w[-1])*r),minor_semiaxis=float(np.sqrt(w[0])*r),
                                      sample_coverage=float(np.mean(md<=r)),major_semiaxis_bootstrap95=np.percentile(boots[:,i],[2.5,97.5],method='linear').tolist())
                          for i,(c,r) in enumerate(zip([68,95],radii))},
                original_sigma_coverage={str(k):float(np.mean(md<=k)) for k in [1,2]},
                wind_distance_correlation=float(np.corrcoef([float(r['wind_fps']) for r in data],distance)[0,1]))

def nominal(path):
    data=rows(path); data=[{k.strip():float(v) for k,v in r.items()} for r in data]
    apo=max(data,key=lambda r:r['Altitude AGL (ft)']); asc=[r for r in data if r['Time']<=apo['Time']]
    powered=[i for i,r in enumerate(data) if r['Thrust (lbf)']>1]
    bo=data[min(powered[-1]+1,len(data)-1)]; peakb=max(asc,key=lambda r:r['Airbrake position']); peakq=max(asc,key=lambda r:r['Dynamic pressure(psf)'])
    return dict(apogee=apo['Altitude AGL (ft)'],t_apogee=apo['Time'],burnout=bo,max_brake=peakb['Airbrake position'],t_max_brake=peakb['Time'],
                max_q_kpa=peakq['Dynamic pressure(psf)']*47.880259/1000,t_max_q=peakq['Time'],landing_distance=data[-1]['Ground track (m)'],t_land=data[-1]['Time'])

def build(project,allow_partial=False):
    dest=project/'results_revised'
    protocol=json.loads((project/'study_protocol.json').read_text(encoding='utf-8'))['jobs']
    expected={j['name']:j['n'] for j in protocol if j['kind'] in ['xml','sil']}
    present={p.name:p for p in dest.glob('mc_*.csv')}
    if set(present)-set(expected): raise ValueError('Unexpected MC result file')
    missing=set(expected)-set(present)
    if missing and not allow_partial: raise ValueError(f'Incomplete Monte Carlo suite: missing {sorted(missing)}')
    mc={name:rows(path,expected[name]) for name,path in present.items()}
    old={p.name:rows(p) for p in (project/'results_v6').glob('mc_*.csv')}
    summary=dict(schema=1,quantile_method='numpy.percentile linear',bootstrap='nonparametric flight resampling; paired differences where applicable; seed7,2000replicates',
                 mc={name:stats(data, next((j for j in [9600,9800,10200,10400] if f'_t{j}_' in name),10000),study_name=name) for name,data in mc.items()},
                 old_mc={name:stats(data, next((j for j in [9600,9800,10200,10400] if f'_t{j}_' in name),10000),study_name=name) for name,data in old.items()},
                 paired={},budget={},landing={},regression={},barometric={})
    pairs=[('mc_closed_v6.csv','mc_oracle_v6.csv'),('mc_sil_truth_nom_v6.csv','mc_closed_v6.csv'),
           ('mc_sil_kf4_est_v6.csv','mc_sil_kf_est_v6.csv'),('mc_sil_kf4_est_v6.csv','mc_sil_kf4_nom_v6.csv'),
           ('mc_sil_kf4_orc_v6.csv','mc_sil_kf4_est_v6.csv'),('mc_sil_kf4_est_notemp_v6.csv','mc_sil_kf4_est_v6.csv'),
           ('mc_sil_kf_est_ideal_v6.csv','mc_oracle_v6.csv'),('mc_sil_kf4_est_ideal_v6.csv','mc_sil_kf_est_ideal_v6.csv'),
           ('mc_sil_kf_est_ideal_fixbrake_v6.csv','mc_sil_kf_est_ideal_v6.csv'),
           ('mc_sil_kf4_est_ideal_fixbrake_v6.csv','mc_sil_kf_est_ideal_fixbrake_v6.csv')]
    pairs.append(('mc_sil_kf4_est_ideal_fixbrake_v6.csv','mc_sil_kf4_est_ideal_v6.csv'))
    for e in ['port','scale','quant','lat','temp','att']: pairs.append((f'mc_sil_kf_est_only_{e}_v6.csv','mc_sil_kf_est_basenoise_v6.csv'))
    for sd in [1,2,5]: pairs.append((f'mc_sil_kf4_est_T{sd}_v6.csv',f'mc_sil_kf4_nom_T{sd}_v6.csv'))
    for a,b in pairs:
        if a in mc and b in mc: summary['paired'][a+' | '+b]=paired(mc[a],mc[b],allow_brake_difference='_fixbrake_' in a or '_fixbrake_' in b)
    if 'mc_open_v6.csv' in mc:
        for name,data in mc.items():
            suffix=next((s for s in ['s10','s30','i01','i05'] if f'_{s}_' in name),None)
            open_name=f'mc_open_{suffix}_v6.csv' if suffix else 'mc_open_v6.csv'
            target=next((t for t in [9600,9800,10200,10400] if f'_t{t}_' in name),10000)
            if open_name in mc: summary['budget'][name]=budget(mc[open_name],data,target)
        for name in ['mc_open_v6.csv','mc_closed_v6.csv']:
            if name in mc: summary['landing'][name]=landing(mc[name])
    if all(name in mc for name in ['mc_open_v6.csv','mc_closed_v6.csv']):
        o,c=mc['mc_open_v6.csv'],mc['mc_closed_v6.csv'];validate_pair(o,c)
        oa,ca=apogees(o),apogees(c);reach=oa>=9900
        def mean_factor(mask,key): return float(np.mean([float(r[key]) for r,take in zip(c,mask) if take])) if mask.any() else None
        low=reach&(ca<9900);high=reach&(ca>10100)
        summary['design_miss_mechanisms']=dict(other_low=int(low.sum()),other_high=int(high.sum()),
              mean_cd_factor_other_low=mean_factor(low,'cd0'),mean_cd_factor_other_high=mean_factor(high,'cd0'))
    summary['temperature_prefix']={}
    for sd,name in [(0,'mc_sil_kf4_est_notemp_v6.csv'),(1,'mc_sil_kf4_est_T1_v6.csv'),(2,'mc_sil_kf4_est_T2_v6.csv'),(3,'mc_sil_kf4_est_v6.csv'),(5,'mc_sil_kf4_est_T5_v6.csv')]:
        if name in mc: summary['temperature_prefix'][str(sd)]=stats(mc[name][:300],study_name=name)
    for name in ['mc_open_v6.csv','mc_closed_v6.csv','mc_sil_kf4_est_v6.csv']:
        if name in mc: summary['regression'][name]=regression(mc[name])
    for name in ['mc_sil_kf4_est_isa_v6.csv','mc_sil_kf4_est_v6.csv']:
        if name in mc: summary['barometric'][name]=stats(mc[name],metric='barometric',study_name=name)
    summary['nominal']={kind:nominal(dest/f'nominal_{kind}_v6.csv') for kind in ['open','closed'] if (dest/f'nominal_{kind}_v6.csv').exists()}
    summary['oat']={p.name:rows(p) for p in dest.glob('oat_*.csv')}
    summary['morris']={p.name:rows(p) for p in dest.glob('morris_*_summary.csv')}
    summary['recovery']=rows(dest/'recovery_v6.csv') if (dest/'recovery_v6.csv').exists() else []
    summary['design_map']=rows(dest/'dmap_v6.csv') if (dest/'dmap_v6.csv').exists() else []
    summary['file_sha256']={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in dest.glob('*.csv') if p.name in mc or p.name in summary['oat'] or p.name in summary['morris']}
    return summary

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--project',type=pathlib.Path,default=ROOT); ap.add_argument('--output',type=pathlib.Path,default=ROOT.parent/'revision_control_summary.json'); ap.add_argument('--allow-partial',action='store_true'); a=ap.parse_args()
    result=build(a.project,a.allow_partial); a.output.write_text(json.dumps(result,indent=2,allow_nan=False),encoding='utf-8')
    lines=['Corrected control statistics','',f"Complete Monte Carlo CSVs analysed: {len(result['mc'])}",'','file | n | mean ft | SD ft | hit% | low% | high%','--- | --- | --- | --- | --- | --- | ---']
    for name,s in result['mc'].items(): lines.append(f"{name} | {s['n']} | {s['mean']:.1f} | {s['sd']:.1f} | {100*s['rates']['hit']:.1f} | {100*s['rates']['low']:.1f} | {100*s['rates']['high']:.1f}")
    lines+=['','RMS definitions: height and speed over every 50 Hz ascent sample; gate-only RMS is a separate metric. Cd RMS is against the instantaneous Mach-dependent equivalent drag coefficient at the assumed filter mass.','Initial attitude SD0.1deg, gyro bias SD0.02deg/s and per-sample gyro noise SD0.1deg/s are engineering assumptions. Ideal ablation disables measurement noise/bias and these attitude errors.','Main impulse perturbation uses fixed propellant mass and Isp scaling; it changes burn duration with the same mass-indexed thrust curve.','All reported ellipse coverages are empirical in-sample fractions. Bootstrap ellipse-axis intervals resample actual flights and do not establish out-of-sample predictive coverage.']
    a.output.with_suffix('.md').write_text('\n'.join(lines),encoding='utf-8'); print('analysed',len(result['mc']),'MC sets')

if __name__=='__main__': main()
