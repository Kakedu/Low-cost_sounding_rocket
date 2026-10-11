# Flight computer in the loop: JSBSim flies the rocket, this script reads
# simulated sensors every 20 ms, runs the Kalman filter + apogee predictor,
# and writes the brake command back.

import argparse
import csv
import math
import random
from pathlib import Path

import jsbsim

import sweep

FT = 0.3048
G0 = 9.80665
SENSORS = dict(baro_sigma_m=0.5, acc_sigma_g=0.04, acc_bias_sigma_g=0.02)
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


def isa_alt_m(p_pa, dT=0.0):
    return (288.15 + dT) / 0.0065 * (1.0 - (p_pa / 101325.0) ** (1.0 / 5.2559))


def isa_rho(h_m, dT=0.0):
    t0 = 288.15 + dT
    t = t0 - 0.0065 * h_m
    return 101325.0 * (t / t0) ** 5.2559 / (287.05 * t)


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
        sens = dict(SENSORS, port_cp_sd=KF_PORT_CP_SD, dT_K=dT_est)
        cls = kf_estimator.FlightKF4
    else:
        sens = dict(SENSORS, port_cp_sd=0.0 if KF_NAIVE else KF_PORT_CP_SD, dT_K=dT_est)
        cls = kf_estimator.FlightKF
    return cls(dt=dt, sensors=sens, cd0=CD_NOMINAL,
                                 dcd_brake=DCD_BRAKE, s_ref=S_REF_M2,
                                 mass_kg=W_BURNOUT_LB * 0.45359237)


def predict_apogee_ft(h_m, v_ms, cd, dT=0.0):
    m = W_BURNOUT_LB * 0.45359237
    vt2 = 2 * m * G0 / (isa_rho(h_m, dT) * S_REF_M2 * cd)
    return (h_m + vt2 / (2 * G0) * math.log(1 + max(v_ms, 0.0) ** 2 / vt2)) / FT


ALL_ERRORS = ("port", "scale", "quant", "lat", "temp")


def run_case(factors, wind=None, est="truth", cd_mode="nominal", rng=None, log=False,
             ideal=False, errors=ALL_ERRORS):
    errors = set(() if ideal else errors)
    rng = rng or random.Random(0)
    sweep.build_case(factors, open_loop=False, wind=wind)
    f = jsbsim.FGFDMExec(str(sweep.JSBSIM_ROOT))
    f.set_debug_level(0)
    f.load_script(f"scripts/{sweep.CASE_NAME}.xml")
    f["atmosphere/delta-T"] = 1.8 * factors.get("dT_K", 0.0)  # JSBSim's delta-T is in Rankine
    f.run_ic()
    f["fcs/sil-mode"] = 1.0
    dt_fc = 1.0 / FC_RATE_HZ
    bias = rng.gauss(0, SENSORS["acc_bias_sigma_g"]) * G0
    port_cp = rng.gauss(REAL["port_cp_mean"], REAL["port_cp_sd"])
    acc_scale = rng.gauss(1.0, REAL["acc_scale_sd"])
    dT_true = factors.get("dT_K", 0.0)
    dT_est = dT_true + rng.gauss(0.0, PAD_TEMP_SD_K)
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
    cd_true = CD_NOMINAL * factors.get("cd0", 1.0)
    steps = int(round(dt_fc / f.get_delta_t()))
    n, burnout, apo, rows, decel_count = 0, False, 0.0, [], 0
    h_err, v_err, cd_err, baro_err = [], [], [], []
    buf = []
    while True:
        if not f.run():
            break
        n += 1
        h = f["position/h-agl-ft"] * FT
        apo = max(apo, h)
        if n % steps:
            continue
        v = f["velocities/h-dot-fps"] * FT
        a_true = f["accelerations/a-pilot-x-ft_sec2"] * FT
        a_raw = acc_scale * a_true + bias + rng.gauss(0, SENSORS["acc_sigma_g"] * G0)
        a_meas = round(a_raw / lsb) * lsb if "quant" in errors else a_raw
        p_port = (f["atmosphere/P-psf"] - port_cp * f["aero/qbar-psf"]) * PSF_TO_PA
        h_meas = isa_alt_m(p_port, dT_est) - h_pad_baro + rng.gauss(0, SENSORS["baro_sigma_m"])
        baro_err.append(h_meas - h)
        buf.append((a_meas, h_meas))
        a_use, h_use = buf[-1 - lag] if len(buf) > lag else buf[0]
        brake = f["fcs/airbrake-pos-norm"]
        est_obj.step(a_use, h_use, dict(h=h, v=v), brake)
        t = f["simulation/sim-time-sec"]
        decel_count = decel_count + 1 if (t > 1.0 and a_use < -0.2 * G0) else 0
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
            h_err.append(est_obj.h - h)
            v_err.append(est_obj.v - v)
            if est_obj.cd is not None:
                cd_err.append(est_obj.cd - cd_true)
        f["fcs/sil-brake-cmd"] = cmd
        if log:
            rows.append(dict(t=t, h=h, v=v, h_meas=h_meas, h_hat=est_obj.h, v_hat=est_obj.v,
                             cd_hat=est_obj.cd, brake=brake, cmd=cmd))
        if burnout and v < 0:
            break
    rms = lambda x: math.sqrt(sum(e * e for e in x) / len(x)) if x else float("nan")
    out = dict(apogee=apo / FT, h_rms_m=rms(h_err), v_rms_ms=rms(v_err), cd_rms=rms(cd_err),
               cd_true=cd_true, dT_est=dT_est, port_cp=port_cp, acc_scale=acc_scale, acc_bias=bias / G0,
               baro_err_max_m=max(baro_err, key=abs) if baro_err else float("nan"))
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
                    help="drag spread other than the literature 0.20 (as sweep.py)")
    ap.add_argument("--ideal", action="store_true",
                    help="ideal sensors: no port error, scale error, quantisation or latency")
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
        r = run_case(f, wind=w, est=a.est, cd_mode=a.cd_mode, rng=noise, ideal=a.ideal, errors=errs)
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
