# RASP time curve to JSBSim burned-mass table

The source RASP thrust points are interpreted as a piecewise-linear function of time, in newtons. Its exact trapezoidal integral is `J`. With source propellant mass `mp`, constant nominal specific impulse is `Isp = J/(mp*g0)`, where `g0 = 9.80665 m/s²`. In the JSBSim pound units this is equivalently `J/(LBF*prop_lb)`. The table coordinate is cumulative expended propellant, not remaining propellant:

`p(t) = mp * integral_0^t F(s) ds / J`.

The input timestamps and a uniform sampling grid are both retained. This preserves each source-curve corner. A model that begins at rail exit starts at `p=0` with the original curve evaluated at that physical elapsed time, and loads only the propellant that remains. The fixed-mass Isp multiplier stretches time while preserving thrust amplitude: `F_scaled(t)=F_source(t/scale)`. Its complete impulse is `scale*J`; propellant loading does not change.

## Endpoint regularization

A zero initial thrust in a mass-indexed table cannot advance the mass coordinate. Near the source endpoint, suppose a time-linear ramp goes from zero to `F1` over `dt`. Its expelled mass in kilograms is `dp=F1*dt/(2*Isp*g0)`. Replace the zero endpoint thrust by `r*F1` in the mass-linear interpolation. Integrating `dt=Isp*g0*dp/F(p)` gives an interval duration of `dt*ln(1/r)/(2*(1-r))`. Matching the original interval duration therefore requires

`ln(1/r)=2*(1-r)`.

The nontrivial root is `r=0.20318786997997995`, recorded as `END_FRACTION` in `motor_model.py`. The same construction regularizes the terminal ramp; the empty tank then shuts off thrust. This constant follows from the endpoint integration equation. It was not fitted to apogee, a flight target or reported experimental altitude. The regularized endpoint force tends to zero with grid spacing, so the approximation must be assessed through sampling and time-step convergence.

The accepted standard uses a 2 ms source-table grid and 1 ms JSBSim step. `python verify_motor.py --scratch SEPARATE_DIRECTORY --json results_revised/verification/motor_convergence_cp.json` checks five table resolutions and two further time-step refinements against the time-domain source. The JSON reports actual initial/final tank mass, measured impulse, depletion time, thrust differences and apogee. These trials end after observed apogee, not at landing. Equal impulse alone does not establish equal acceleration history or altitude, and table-resolution convergence does not replace integrator refinement.

`start_time` is physical elapsed motor time including the Isp stretch. The helper rejects negative force, non-increasing times and zero-flow plateaus; an arbitrary reigniting curve with a zero-thrust interval requires another engine representation.
