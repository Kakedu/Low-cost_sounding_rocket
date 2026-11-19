"""Reproducible hinge-frequency screen and concrete assembly measurements.

Uses accepted flap geometry/FE frequencies and loads; adds no fitted
coefficient and makes no claim of assembly certification or flutter.
"""
from pathlib import Path
import argparse,hashlib,json,math
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    root=Path(__file__).resolve().parent;ap=argparse.ArgumentParser()
    ap.add_argument('--fea',type=Path,default=root/'results_revised/flap_fea_v6.json')
    ap.add_argument('--loads',type=Path,default=root/'results_revised/loads_v6.json')
    ap.add_argument('--json',type=Path,default=root/'results_revised/structure_hinge_screen.json')
    ap.add_argument('--markdown',type=Path,default=root/'STRUCTURE_VALIDATION.md');args=ap.parse_args()
    fe=json.loads(args.fea.read_text());loads=json.loads(args.loads.read_text());dy=fe['dynamics']
    inertia=dy['I_hinge_kgm2'];f1=fe['meshes'][-1]['f_Hz'][0]
    shedding=[value for velocities in dy['shedding_Hz'].values() for value in velocities.values()]
    low,high=min(shedding),max(shedding)
    threshold=lambda f:inertia*(2*math.pi*f)**2
    kf1=threshold(f1);k2=threshold(2*f1)
    if not math.isclose(k2,dy['k_for_2x_f1_Nm_per_rad'],rel_tol=1e-12):raise RuntimeError('FEA hinge threshold formula disagrees')
    stiffnesses=[('assumed k=1',1.),('assumed k=10',10.),('lower shedding-envelope threshold',threshold(low)),('assumed k=100',100.),('first plate-mode threshold',kf1),('assumed k=1000',1000.),('twice first plate-mode threshold',k2),('upper shedding-envelope threshold',threshold(high))]
    table=[]
    for label,k in stiffnesses:
        fh=math.sqrt(k/inertia)/(2*math.pi)
        table.append(dict(label=label,hinge_stiffness_Nm_per_rad=k,rigid_hinge_frequency_Hz=fh,frequency_over_plate_f1=fh/f1,frequency_over_twice_plate_f1=fh/(2*f1),within_shedding_envelope=low-1e-9<=fh<=high+1e-9))
    # Match the frozen sampler: B=max(0.5,X), X~N(1,0.20^2).
    # This is lower clipping (with a point mass at the floor), not resampling.
    nominal_mean=1.;sd=.2;floor=.5;factor=1.4
    cdf=lambda z:.5*math.erfc(-z/math.sqrt(2))
    phi=lambda z:math.exp(-.5*z*z)/math.sqrt(2*math.pi)
    a=(floor-nominal_mean)/sd;floor_probability=cdf(a)
    clipped_mean=nominal_mean*(1-floor_probability)+sd*phi(a)+floor*floor_probability
    second_moment=(nominal_mean**2+sd**2)*(1-floor_probability)+sd*(nominal_mean+floor)*phi(a)+floor**2*floor_probability
    clipped_sd=math.sqrt(second_moment-clipped_mean**2)
    coefficient_coverage=cdf((factor*nominal_mean-nominal_mean)/sd)
    evidence=root/'provenance/brake_coefficient_distribution_clarification.json'
    clarification=json.loads(evidence.read_text())
    if clarification['status']!='CONFIRMED_FROZEN_SAMPLER' or not clarification['all_draw_inputs_AST_identical']:
        raise RuntimeError('Missing frozen brake-coefficient sampler evidence')
    for value,key in [(floor_probability,'floor_probability'),(clipped_mean,'actual_clipped_mean'),
                      (coefficient_coverage,'probability_B_le_1_4_nominal')]:
        if not math.isclose(value,clarification[key],rel_tol=1e-12,abs_tol=1e-15):
            raise RuntimeError('Frozen sampler evidence and clipped-normal calculation disagree: '+key)
    values=dict(status='SCREEN_ONLY',hinge_inertia_kg_m2=inertia,first_plate_frequency_Hz=f1,shedding_envelope_Hz=[low,high],formula='f_h=sqrt(k_theta/I_h)/(2*pi); k_theta(target)=I_h*(2*pi*f_target)^2',stiffness_table=table,
     load_factor_choice=dict(coefficient_factor=factor,additional_load_factor=1.5,combined_factor=2.1,
      coefficient_cv_assumed=sd/nominal_mean,coefficient_distribution='B=max(0.5,X), X~Normal(mean=1,SD=0.20)',
      underlying_normal_mean=nominal_mean,underlying_normal_sd=sd,lower_clip=floor,
      lower_clip_probability=floor_probability,actual_clipped_mean=clipped_mean,actual_clipped_sd=clipped_sd,
      actual_clipped_cv=clipped_sd/clipped_mean,
      marginal_coefficient_coverage_at_1p4_nominal=coefficient_coverage,
      coverage_formula='P(B<=1.4)=Phi((1.4-1)/0.20)=Phi(2), because 1.4 is above the 0.5 lower clip.',
      scope='Engineering choices of this study. 1.4 is nominal mean plus 2 assigned SD and has 97.725% marginal coefficient coverage under the chosen clipped-normal model. Lower clipping slightly changes its actual mean and CV. This coefficient probability does not cover complete flight loads, temporal maxima or all uncertain inputs, and the additional 1.5 multiplier and 2.1 product are not calibrated global reliability or capacity bounds.'),
     source_sha256={str(args.fea):sha(args.fea),str(args.loads):sha(args.loads),'structure_screen.py':sha(Path(__file__)),
      'provenance/brake_coefficient_distribution_clarification.json':sha(evidence)},
     load_allocation='Four flaps share modeled brake force equally. Four-fin set divided by two is a conservative maximum-single-fin envelope for two fins in an incidence plane; combined alpha/beta loads do not resolve each physical fin. Fin incidence is capped at 15deg, and actual/used angles and exceedance flags are separate.',
     limitations=['Frequency overlap is an envelope screen; trajectory-synchronous speed, deployment height and forcing spectrum are not jointly evaluated.','Rigid hinge frequency assumes a single flap inertia and linear rotational spring. Coupled plate, hinge, servo, linkage and body modes can differ.','The factor two relative to plate frequency is a stated screening separation criterion, not an established stability or certification threshold.','Uniform pressure and a fully fixed edge do not establish hinge-pin, lug, fastener, contact, actuator, fatigue or aeroelastic capacity.'],
     verification_plan=[
      dict(step=1,measure='Actual geometry and mass',method='Record flap thickness/width/length, hinge axis and pin diameter, lug spacing, fastener locations, servo arm, linkage angles, and deployed projected height. Weigh each flap and measure or compute its hinge-axis inertia from the actual assembly.',output='Dimensioned assembly, mass/inertia uncertainty, and load-path coordinates.'),
      dict(step=2,measure='Effective rotational stiffness and backlash',method='Measure applied hinge torque versus angle through the relevant deployment range with the actuator powered and unpowered. Repeat loading/unloading to identify backlash, hysteresis and compliance; calculate the local derivative dM/dtheta in N m/rad.',output='Measured k_theta(theta), backlash/hysteresis curves, and an updated coupled-mode model.'),
      dict(step=3,measure='Assembly modes and damping',method='Measure a low-amplitude impulse or frequency-response test of the assembled flap/hinge/linkage/body in representative deployed states. Fit observed natural frequencies and damping; compare to the plate-only and rigid-hinge predictions.',output='Measured coupled-mode frequencies, damping and mode shapes; no stability claim from frequency separation alone.'),
      dict(step=4,measure='Static load path and actuator torque',method='Use the accepted modeled load as a calculation input and compare distributed-pressure and representative linkage/fastener load paths in an assembly model. Instrument a bench specimen to compare deflection, strain and actuator torque with those predictions under a separately approved hardware test procedure.',output='Load/deflection and torque data with geometry, fixture and force distribution documented. The existing 2.1 factor remains a screen, not certified capacity.'),
      dict(step=5,measure='Flow forcing at matched deployment states',method='Use synchronized flight-model velocity and projected deployed height to replace the rectangular shedding envelope. Representative flow testing should measure pressure/force spectra and response before any flutter claim.',output='Matched velocity/deployment/forcing spectra and a coupled aeroelastic assessment.')])
    args.json.parent.mkdir(parents=True,exist_ok=True);args.json.write_text(json.dumps(values,indent=2),encoding='utf-8')
    rows='\n'.join(f'| {r["label"]} | {r["hinge_stiffness_Nm_per_rad"]:.3f} | {r["rigid_hinge_frequency_Hz"]:.2f} | {r["frequency_over_plate_f1"]:.3f} | {"yes" if r["within_shedding_envelope"] else "no"} |' for r in table)
    text=f'''# Flap, hinge and actuator verification screen

The plate-only first frequency is {f1:.2f} Hz. The assumed flap hinge inertia is {inertia:.9g} kg m². The rigid-hinge relation is f_h = sqrt(k_theta/I_h)/(2π). The assumed shedding relation f=0.146 V/h over V=100–270 m/s and h=10–55 mm gives {low:.2f}–{high:.2f} Hz. These rectangular envelopes are not matched flight/deployment histories.

| stiffness case | k_theta (N m/rad) | rigid hinge f (Hz) | f/f_plate1 | in shedding envelope |
|---|---:|---:|---:|---|
{rows}

The assumed 1 and 10 N m/rad hinge cases lie below this shedding envelope; 100 and 1000 N m/rad lie within it. A rigid hinge frequency equal to the first plate mode corresponds to {kf1:.2f} N m/rad, and twice that mode corresponds to {k2:.2f} N m/rad. The factor-two separation is a screening convention, not a validated stability condition. Changing the actual assembly inertia, compliance or damping changes the result.

Four modeled flaps share total brake force equally. Dividing the four-fin set load by two is a maximum-single-fin screening envelope for two fins in one incidence plane; it does not resolve each fin under combined incidence. Fin incidence is capped at 15°, with actual/used angles and exceedance flags recorded separately. Uniform-pressure plate FEA excludes the complete hinge, contact, mounting, actuator and linkage load paths.

The design force is the maximum among seven selected simulation scenarios, not a physical maximum over all uncertainty. The study chooses a 1.4 brake-coefficient factor (nominal mean plus two assigned standard deviations) and a 1.5 additional multiplier, giving 2.1. The frozen sampler uses B=max(0.5,X), X~Normal(1,0.20²), with lower clipping rather than resampling. Because 1.4 is above the floor, its marginal coefficient coverage is Phi((1.4-1)/0.20)=Phi(2)={100*coefficient_coverage:.3f}%. The floor has probability {100*floor_probability:.5f}%, and clipping gives actual mean {clipped_mean:.9f} and SD {clipped_sd:.9f}; the assigned underlying CV is 20%. This probability concerns that single coefficient, not complete flight loads, temporal maxima or all uncertain inputs. Neither the additional 1.5 engineering multiplier nor the 2.1 product guarantees global load coverage or certifies capacity.

## Measurements needed before an assembly claim

'''+ '\n\n'.join(f'{i+1}. **{v["measure"]}.** {v["method"]} Required evidence: {v["output"]}' for i,v in enumerate(values['verification_plan']))+'''

Run `python structure_screen.py` after accepted `loads_v6.json`, `flap_fea_v6.json` and `provenance/brake_coefficient_distribution_clarification.json` are present. Output `results_revised/structure_hinge_screen.json` binds the accepted load/FE inputs, frozen-sampler evidence and script SHA256. No new fitted factors, FEA solve, physical packaging proposal or assembly certification is introduced.
'''
    args.markdown.write_text(text,encoding='utf-8');print(f'Wrote {len(table)} hinge stiffness cases and five hardware measurement steps')
if __name__=='__main__':main()
