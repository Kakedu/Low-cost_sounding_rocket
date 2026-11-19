"""Explicit source fields for the main-report RQ2 numerical template."""
import csv,hashlib,re
from pathlib import Path
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def rows(path):return list(csv.DictReader(path.open()))
def add(summary,root):
    root=Path(root);d=root/'results_revised';log=d/'logs/physics/crosscheck_authority_csv.log'
    text=log.read_text(encoding='utf-8');m=re.search(r'First 4 s give (\d+) ft of the (\d+) ft total, (\d+)%',text)
    if not m:raise ValueError('Missing source-bound single-run first4 authority')
    estimator=rows(d/'crosscheck_estimator.csv');control=rows(d/'crosscheck_control.csv');recovery=rows(d/'recovery_v6.csv')
    valid=[float(r['impulse_factor']) for r in recovery if abs(float(r['error_ft']))<=100]
    reachable=[r for r in control if r['case']!='Isp -5%']
    nominal=summary['nominal']['open']['apogee']
    result=dict(selected_filter='kf4',selected_filter_definition='Predefined mainRQ2 four-state estimator with Cp state; selected before outcomes, not by posthoc highest hit rate.',
        brake_first4_gain_ft=float(m[1]),brake_full_coast_gain_ft=float(m[2]),brake_first4_authority_percent=float(m[3]),
        authority_scope='Single-run point-mass internal cross-check; values rounded in captured originalstdout. Window effects are not additive.',
        predictor_fixed_mae_ft=sum(abs(float(r['const_err_ft'])) for r in estimator)/len(estimator),
        predictor_closed_mae_ft=sum(abs(float(r['dragform_err_ft'])) for r in estimator)/len(estimator),
        predictor_fixed_worst_ft=max(abs(float(r['const_error_ft'])) for r in control),
        predictor_closed_worst_ft=max(abs(float(r['dragform_error_ft'])) for r in control),
        predictor_reachable_fixed_worst_ft=max(abs(float(r['const_error_ft'])) for r in reachable),
        predictor_reachable_closed_worst_ft=max(abs(float(r['dragform_error_ft'])) for r in reachable),
        predictor_scope='Worst all7 cases includes the energy-short Isp-5% case; reachable worst excludes that single unbraked-below9900ft case, retaining drymass+5%.',
        isp_recovery_min_factor=min(valid),isp_recovery_max_factor=max(valid),isp_recovery_factors=valid,
        recovery_scope='Discrete tested Isp factors at fixed fuel, XML controller; not a continuous guaranteed range.',
        margin_min_percent=100*(nominal-10400)/10400,margin_max_percent=100*(nominal-9600)/9600,
        source_sha256={p.relative_to(root).as_posix():sha(p) for p in [log,d/'crosscheck_estimator.csv',d/'crosscheck_control.csv',d/'recovery_v6.csv',d/'nominal_open_v6.csv']})
    summary['document_metrics']=result;return summary
