# Quick checks of the filter with ideal sensors.

import math
import random
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent))

import kf_estimator as K

G0 = 9.80665
SENS = dict(baro_sigma_m=0.5, acc_sigma_g=0.04, acc_bias_sigma_g=0.02)
S, M, DCD = 0.0845 * 0.09290304, 38.8 * 0.45359237, 2.0


def test_synthetic(cd_true=0.55, seed=3):
    rng = random.Random(seed)
    dt = 0.02
    kf = K.FlightKF(dt, SENS, 0.45, DCD, S, M)
    h, v = 580.0, 262.0
    # This synthetic test starts mid-flight with prescribed initial h/v. A
    # gate-enabled filter must not be tested as though that first observation
    # were a 580 m outlier from a launch-at-zero prior. Production SIL retains
    # its original launch prior and is checked separately below.
    kf.x[:2] = [h, v]
    kf.h, kf.v = h, v
    bias = rng.gauss(0, SENS["acc_bias_sigma_g"]) * G0
    t, err2s = 0.0, None
    while v > 0:
        a = -K.isa_rho(h) * v * v * S * cd_true / (2 * M)
        kf.step(a + bias + rng.gauss(0, SENS["acc_sigma_g"] * G0),
                h + rng.gauss(0, SENS["baro_sigma_m"]), None, 0.0)
        v += (a - G0) * dt
        h += v * dt
        t += dt
        if err2s is None and t >= 2.0 + 3 * dt:
            err2s = abs(kf.cd - cd_true) / cd_true
    ok = err2s is not None and err2s < 0.05
    print(f"[{'PASS' if ok else 'FAIL'}] synthetic coast: Cd error 2 s after burnout {100*err2s:.1f}% (limit 5%)")
    return ok


ROOT = Path(__import__("os").environ.get("JSBSIM_ROOT", HERE.parent))


def test_sil():
    import sil
    for k, v in sil.sweep._reroot(ROOT).items():
        setattr(sil.sweep, k, v)
    ref = sil.run_case({}, est="truth", cd_mode="nominal")
    out = sil.run_case({}, est="kf", cd_mode="estimated", rng=random.Random(5), ideal=True)
    ok2 = out["h_rms_m"] < 3 and out["v_rms_ms"] < 3 and abs(out["apogee"] - ref["apogee"]) < 60
    print(f"[{'PASS' if ok2 else 'FAIL'}] SIL nominal: h rms {out['h_rms_m']:.2f} m, v rms {out['v_rms_ms']:.2f} m/s, "
          f"apogee {out['apogee']:.0f} ft vs truth-estimator {ref['apogee']:.0f} ft")
    # An underpowered high-drag flight cannot be repaired by adding drag. Use
    # an overpowered low-drag case to test the estimator's control benefit.
    hi_nom = sil.run_case({"cd0": 0.8}, est="kf", cd_mode="nominal", rng=random.Random(6), ideal=True)
    hi_est = sil.run_case({"cd0": 0.8}, est="kf", cd_mode="estimated", rng=random.Random(6), ideal=True)
    ok3 = abs(hi_est["apogee"] - 10000) < abs(hi_nom["apogee"] - 10000)
    print(f"[{'PASS' if ok3 else 'FAIL'}] drag -20%: apogee error with estimated Cd {hi_est['apogee']-10000:+.0f} ft, "
          f"with nominal Cd {hi_nom['apogee']-10000:+.0f} ft")
    return ok2 and ok3


if __name__ == "__main__":
    a = test_synthetic()
    b = test_sil()
    print("ALL PASS" if a and b else "NOT YET")
    sys.exit(0 if a and b else 1)
