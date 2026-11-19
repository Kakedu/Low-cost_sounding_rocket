# Flight computer in the loop: JSBSim flies the rocket, this script reads
# simulated sensors every 20 ms, runs the Kalman filter + apogee predictor,
# and writes the brake command back.

import argparse
import csv
import json
import math
import random
from pathlib import Path

import jsbsim

import sweep
from kf.attitude_sensor import AttitudeSensor, named_rng, quaternion_from_euler, up_projection

FT = 0.3048
G0 = 9.80665
SENSORS = dict(baro_sigma_m=0.5, acc_sigma_g=0.04, acc_bias_sigma_g=0.02)
# A negative specific force is possible while drag exceeds residual thrust.
# Require the known nominal motor duration plus three assigned Isp SDs before
# counting the three negative-force packets. This is a prior-based protection,
# not knowledge of the sampled motor's actual depletion time. The unbounded
# Normal Isp scenario still has a nonzero probability beyond this time.
_NOMINAL_BURN_S = json.loads((Path(__file__).resolve().parent / "motor_build_info.json").read_text())["burn_end_s"]
COAST_GUARD_S = _NOMINAL_BURN_S * (1.0 + 3.0 * sweep.MC_SIGMA["isp"])
REAL = dict(port_cp_mean=0.0, port_cp_sd=0.03, acc_scale_sd=0.03, acc_lsb_g=0.049,
            latency_cycles=1)
PSF_TO_PA = 47.880259
FC_RATE_HZ = 50
CD_NOMINAL = 0.45
DCD_BRAKE = 2.0
S_REF_M2 = 0.0845 * 0.09290304
W_BURNOUT_LB = 38.8
KP_PER_FT = 0.001
TARGET_FT = 10000.0
GATE_FT = 1500.0
CD_MIN, CD_MAX = 0.2, 1.0  # clamp the filter's Cd before the predictor uses it


from atmosphere_model import troposphere_density_kg_m3

def isa_alt_m(p_pa, dT=0.0):
    T0 = 288.15 + dT
    H = T0 / 0.0065 * (1.0 - (p_pa / 101325.0) ** (287.05287 * 0.0065 / G0))
    return H * 6356766.0 / (6356766.0 - H)

def isa_rho(h_m, dT=0.0):
    return troposphere_density_kg_m3(h_m, dT)


class TruthEstimator:
    def __init__(self, dt):
        self.h = self.v = 0.0
        self.cd = None

    def step(self, a_meas, h_meas, truth, brake):
        self.h, self.v = truth["h"], truth["v"]


KF_PORT_CP_SD = 0.03
CTRL_DELAY_S = 0.0
KF_NAIVE = True
PAD_TEMP_SD_K = 3.0  # how well the flight computer knows the air temperature
USE_PAD_TEMP = True


def load_estimator(name, dt, dT_est=0.0):
    if name == "truth":
        return TruthEstimator(dt)
    import sys
    sys.path.insert(0, str(Path(__file__).resolve().parent / "kf"))
    import kf_estimator
    if name == "kf4":
        sens = dict(SENSORS, port_cp_sd=KF_PORT_CP_SD, dT_K=dT_est, coast_min_time_s=COAST_GUARD_S)
        cls = kf_estimator.FlightKF4
    else:
        sens = dict(SENSORS, port_cp_sd=0.0 if KF_NAIVE else KF_PORT_CP_SD, dT_K=dT_est, coast_min_time_s=COAST_GUARD_S)
        cls = kf_estimator.FlightKF
    return cls(dt=dt, sensors=sens, cd0=CD_NOMINAL,
                                 dcd_brake=DCD_BRAKE, s_ref=S_REF_M2,
                                 mass_kg=W_BURNOUT_LB * 0.45359237)


def predict_apogee_ft(h_m, v_ms, cd, dT=0.0):
    m = W_BURNOUT_LB * 0.45359237
    vt2 = 2 * m * G0 / (isa_rho(h_m, dT) * S_REF_M2 * cd)
    return (h_m + vt2 / (2 * G0) * math.log(1 + max(v_ms, 0.0) ** 2 / vt2)) / FT


ALL_ERRORS = ("port", "scale", "quant", "lat", "temp", "att")


def run_case(factors, wind=None, est="truth", cd_mode="nominal", rng=None, log=False,
             ideal=False, errors=ALL_ERRORS):
    errors = set(() if ideal else errors)
    rng = rng or random.Random(0)
    sweep.build_case(factors, open_loop=False, wind=wind)
    f = jsbsim.FGFDMExec(str(sweep.JSBSIM_ROOT))
    f.set_debug_level(0)
    if not f.load_script(f"scripts/{sweep.CASE_NAME}.xml"):
        raise RuntimeError("SIL script failed to load")
    f["atmosphere/delta-T"] = 1.8 * factors.get("dT_K", 0.0)  # JSBSim's delta-T is in Rankine
    if not f.run_ic():
        raise RuntimeError("SIL initialisation failed")
    # SIL metrics and optional sensor traces are recorded below; the separate
    # 50 Hz FDM CSV is unused and would dominate repeated Monte Carlo I/O.
    f.disable_output()
    f["fcs/sil-mode"] = 1.0
    dt_fc = 1.0 / FC_RATE_HZ
    sensor_seed = rng.getrandbits(64)
    acc_rng = named_rng(sensor_seed, "accelerometer")
    baro_rng = named_rng(sensor_seed, "barometer")
    calibration_rng = named_rng(sensor_seed, "calibration")
    biases = [calibration_rng.gauss(0, SENSORS["acc_bias_sigma_g"])*G0 for _ in range(3)]
    if ideal:
        biases = [0.0]*3
    bias = biases[0]
    attitude = AttitudeSensor(f["attitude/phi-rad"], f["attitude/theta-rad"],
                              f["attitude/psi-rad"], dt_fc, sensor_seed, "att" in errors)
    port_cp = calibration_rng.gauss(REAL["port_cp_mean"], REAL["port_cp_sd"])
    acc_scale = calibration_rng.gauss(1.0, REAL["acc_scale_sd"])
    dT_true = factors.get("dT_K", 0.0)
    dT_est = dT_true + calibration_rng.gauss(0.0, PAD_TEMP_SD_K)
    if "port" not in errors:
        port_cp = 0.0
    if "scale" not in errors:
        acc_scale = 1.0
    if "temp" not in errors:
        dT_est = dT_true
    elif not USE_PAD_TEMP:
        dT_est = 0.0
    est_obj = load_estimator(est, dt_fc, dT_est)
    lsb = REAL["acc_lsb_g"] * G0
    lag = REAL["latency_cycles"] if "lat" in errors else 0
    p_pad = f["atmosphere/P-psf"] * PSF_TO_PA
    h_pad_baro = isa_alt_m(p_pad, dT_est)
    h_pad_agl = f["position/h-agl-ft"]*FT
    cd_true = CD_NOMINAL * factors.get("cd0", 1.0)
    steps = int(round(dt_fc / f.get_delta_t()))
    n, burnout, apo, rows, decel_count = 0, False, 0.0, [], 0
    h_err, v_err, cd_err, baro_err = [], [], [], []
    buf = []
    gated_h_err, gated_v_err, cd_low_err, att_err = [], [], [], []
    ascent_finished = False
    t_apogee = None
    while True:
        if not f.run():
            break
        n += 1
        if ascent_finished:
            f["fcs/sil-brake-cmd"] = 0.0
            continue
        h = f["position/h-agl-ft"] * FT
        if not math.isfinite(h):raise RuntimeError("Nonfinite plant height during SIL ascent")
        apo = max(apo, h)
        if n % steps:
            continue
        v = f["velocities/h-dot-fps"] * FT
        body_a = [f[f"accelerations/a-pilot-{axis}-ft_sec2"]*FT for axis in "xyz"]
        rates = [f[f"velocities/{axis}-rad_sec"] for axis in "pqr"]
        if not all(math.isfinite(x) for x in [h,v,*body_a,*rates]):raise RuntimeError("Nonfinite plant/sensor state during SIL ascent")
        attitude.step(rates)
        measured_a = []
        for a_body, axis_bias in zip(body_a, biases):
            white = acc_rng.gauss(0, SENSORS["acc_sigma_g"]*G0)
            raw = acc_scale*a_body + axis_bias + (0.0 if ideal else white)
            measured_a.append(round(raw/lsb)*lsb if "quant" in errors else raw)
        a_meas = attitude.specific_force_up(measured_a)
        q_truth = quaternion_from_euler(f["attitude/phi-rad"], f["attitude/theta-rad"], f["attitude/psi-rad"])
        qdot = abs(sum(a*b for a,b in zip(q_truth, attitude.q)))
        att_err.append(math.degrees(2*math.acos(min(1.0,qdot))))
        p_port = (f["atmosphere/P-psf"] - port_cp * f["aero/qbar-psf"]) * PSF_TO_PA
        h_meas = isa_alt_m(p_port, dT_est) - h_pad_baro + h_pad_agl + (0.0 if ideal else baro_rng.gauss(0, SENSORS["baro_sigma_m"]))
        baro_err.append(h_meas - h)
        buf.append((a_meas, h_meas))
        a_use, h_use = buf[-1 - lag] if len(buf) > lag else buf[0]
        brake = f["fcs/airbrake-pos-norm"]
        est_obj.step(a_use, h_use, dict(h=h, v=v), brake)
        if not all(math.isfinite(x) for x in [est_obj.h,est_obj.v,a_use,h_use]):raise RuntimeError("Nonfinite SIL measurement or estimated state")
        if hasattr(est_obj,'P'):
            import numpy as np
            if not np.isfinite(est_obj.x).all() or not np.isfinite(est_obj.P).all():raise RuntimeError("Nonfinite EKF state or covariance")
        h_err.append(est_obj.h - h)
        v_err.append(est_obj.v - v)
        t = f["simulation/sim-time-sec"]
        decel_count = decel_count + 1 if (t + 1e-10 >= COAST_GUARD_S and a_use < -0.2 * G0) else 0
        if not burnout and decel_count >= 3:
            burnout, t_bo = True, t
        cmd = 0.0
        if burnout and t - t_bo >= CTRL_DELAY_S and est_obj.v > 0 and est_obj.h / FT > GATE_FT:
            if cd_mode == "estimated" and est_obj.cd is not None:
                cd = min(CD_MAX, max(CD_MIN, est_obj.cd))
            elif cd_mode == "oracle":
                cd = cd_true
            else:
                cd = CD_NOMINAL
            pred = predict_apogee_ft(est_obj.h, est_obj.v, cd, dT_est)
            cmd = min(1.0, max(0.0, KP_PER_FT * (pred - TARGET_FT)))
            gated_h_err.append(est_obj.h - h)
            gated_v_err.append(est_obj.v - v)
            if est_obj.cd is not None:
                mach = f["velocities/mach"]
                machs = [0.0,0.5,0.8,0.95,1.05,1.3,2.0,3.0]
                cds = [0.45,0.45,0.48,0.72,0.78,0.70,0.58,0.48]
                import numpy as np
                body_cd = float(np.interp(mach,machs,cds))*factors.get("cd0",1.0)
                equivalent_cd = body_cd*W_BURNOUT_LB/f["inertia/weight-lbs"]
                cd_err.append(est_obj.cd-equivalent_cd)
                cd_low_err.append(est_obj.cd-cd_true)
        f["fcs/sil-brake-cmd"] = cmd
        if log:
            rows.append(dict(t=t, h=h, v=v, h_meas=h_meas, h_hat=est_obj.h, v_hat=est_obj.v,
                             cd_hat=est_obj.cd, brake=brake, cmd=cmd))
        if burnout and v < 0:
            ascent_finished = True
            t_apogee = t
    terminal = [f["position/h-agl-ft"], f["velocities/h-dot-fps"], f["simulation/sim-time-sec"]]
    if not all(math.isfinite(x) for x in terminal) or not ascent_finished or terminal[0] >= 1.1 or terminal[2] >= 599.9:
        raise RuntimeError("SIL flight did not reach apogee and complete landing")
    rms = lambda x: math.sqrt(sum(e * e for e in x) / len(x)) if x else float("nan")
    out = dict(apogee=apo / FT, h_rms_m=rms(h_err), v_rms_ms=rms(v_err), cd_rms=rms(cd_err),
               h_gate_rms_m=rms(gated_h_err), v_gate_rms_ms=rms(gated_v_err),
                h_max_abs_m=max(map(abs,h_err)), v_max_abs_ms=max(map(abs,v_err)),
                cd_low_mach_rms=rms(cd_low_err), attitude_rms_deg=rms(att_err),
                t_apogee=t_apogee, t_land=f["simulation/sim-time-sec"], ascent_sensor_samples=len(h_err),
                gate_sensor_samples=len(gated_h_err), cd_projected_updates=getattr(est_obj,"cd_projection_count",0),
                coast_guard_s=COAST_GUARD_S, controller_coast_detected_s=t_bo,
                filter_coast_detected_s=getattr(est_obj,"coast_detected_time_s",t_bo),
                baro_updates=getattr(est_obj,"measurement_updates",{}).get("baro",0),
                accel_updates=getattr(est_obj,"measurement_updates",{}).get("accel",0),
                baro_rejections=getattr(est_obj,"measurement_rejections",{}).get("baro",0),
                accel_rejections=getattr(est_obj,"measurement_rejections",{}).get("accel",0),
                baro_max_nis=getattr(est_obj,"max_nis",{}).get("baro",0.0),
                accel_max_nis=getattr(est_obj,"max_nis",{}).get("accel",0.0),
                cd_true=cd_true, dT_est=dT_est, port_cp=port_cp, acc_scale=acc_scale, acc_bias=bias / G0,
               baro_err_max_m=max(baro_err, key=abs) if baro_err else float("nan"))
    for key,value in out.items():
        undefined_cd=est_obj.cd is None and key in ['cd_rms','cd_low_mach_rms']
        undefined_gate=not gated_h_err and key in ['h_gate_rms_m','v_gate_rms_ms','cd_rms','cd_low_mach_rms']
        if not math.isfinite(value) and not (undefined_cd or undefined_gate):raise RuntimeError(f"Nonfinite applicable SIL output metric {key}")
    return (out, rows) if log else out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=["nominal", "mc"])
    ap.add_argument("n", nargs="?", type=int, default=500)
    ap.add_argument("--est", default="truth")
    ap.add_argument("--cd", default="nominal", dest="cd_mode")
    ap.add_argument("--root", default=".")
    ap.add_argument("--csv", default=None)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--target", type=float, default=None, help="target apogee, ft (default 10000)")
    ap.add_argument("--cd-sigma", type=float, default=None, dest="cd_sigma",
                    help="drag spread other than the assigned 0.20 scenario (as sweep.py)")
    ap.add_argument("--ideal", action="store_true",
                    help="ideal sensors: all measurement noise, biases and listed errors disabled")
    ap.add_argument("--errors", default=",".join(ALL_ERRORS),
                    help="comma list of v6 sensor errors to switch on (default all; '' = none)")
    ap.add_argument("--kf-q-aware", action="store_true", dest="kf_q",
                    help="the filter trusts the barometer less at high dynamic pressure")
    ap.add_argument("--no-pad-temp", action="store_true", dest="no_pad_temp",
                    help="the flight computer assumes ISA instead of the measured pad temperature")
    ap.add_argument("--pad-temp-sd", type=float, default=None, dest="pad_temp_sd",
                    help="sd of the flight computer's knowledge of the column temperature, K (default 3)")
    ap.add_argument("--fix-brake", action="store_true", dest="fix_brake",
                    help="set the true brake increment to nominal (isolates the brake's uncertainty)")
    ap.add_argument("--ctrl-delay", type=float, default=0.0, dest="ctrl_delay",
                    help="start control this many seconds after the end of boost is detected")
    a = ap.parse_args()
    global KF_NAIVE, USE_PAD_TEMP, CTRL_DELAY_S, PAD_TEMP_SD_K
    if a.pad_temp_sd is not None:
        PAD_TEMP_SD_K = a.pad_temp_sd
    KF_NAIVE = not a.kf_q
    USE_PAD_TEMP = not a.no_pad_temp
    CTRL_DELAY_S = a.ctrl_delay
    errs = tuple(e for e in a.errors.split(",") if e)
    assert all(e in ALL_ERRORS for e in errs), errs
    if a.target is not None:
        global TARGET_FT
        TARGET_FT = a.target
    for k, v in sweep._reroot(Path(a.root)).items():
        setattr(sweep, k, v)
    if a.mode == "nominal":
        out, rows = run_case({}, est=a.est, cd_mode=a.cd_mode, log=True, ideal=a.ideal, errors=errs)
        print(out)
        if a.csv:
            with open(a.csv, "w", newline="") as fh:
                w = csv.DictWriter(fh, fieldnames=list(rows[0]))
                w.writeheader(); w.writerows(rows)
        return
    rng = random.Random(a.seed)
    rows = []
    sig = dict(sweep.MC_SIGMA)
    if a.cd_sigma is not None:
        sig["cd0"] = a.cd_sigma
    for i in range(a.n):
        f, w = sweep.draw_inputs(rng, sig)
        if a.fix_brake:
            f["airbrake_dcd"] = 1.0
        noise = random.Random(a.seed * 1_000_003 + i)
        try:
            r = run_case(f, wind=w, est=a.est, cd_mode=a.cd_mode, rng=noise, ideal=a.ideal, errors=errs)
        except Exception as error:
            import json
            failure = dict(status="FAILED", requested=a.n, completed=len(rows), seed=a.seed,
                           sample=i, factors=f, wind=w, sensor_rng_seed=a.seed*1_000_003+i,
                           exception=repr(error))
            if a.csv:
                Path(a.csv+".failure.json").write_text(json.dumps(failure,indent=2),encoding="utf-8")
            raise RuntimeError(f"SIL Monte Carlo failed at sample {i}; no partial success CSV is written") from error
        rows.append(dict(sample=i, **{k: f[k] for k in sorted(f)}, wind_fps=w[0], wind_from=w[1], **r))
        if (i + 1) % 25 == 0:
            print(f"  {i + 1}/{a.n}", flush=True)
    if a.csv:
        with open(a.csv, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0]))
            w.writeheader(); w.writerows(rows)
        print("wrote", a.csv)
    ap_ = [r["apogee"] for r in rows]
    print(f"target {TARGET_FT:.0f}  mean {sum(ap_)/len(ap_):.0f}"
          f"  within100 {sum(abs(x-TARGET_FT)<=100 for x in ap_)/len(ap_):.3f}"
          f"  below {sum(x<TARGET_FT-100 for x in ap_)/len(ap_):.3f}")


if __name__ == "__main__":
    main()
