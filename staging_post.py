# Post-processing of the staging Monte Carlo: gyro drift on the tilt check
# and the separation-timer margin.

import argparse
import csv
import math
import random
import statistics as st
import json
from pathlib import Path

BIASES = (0.01, 0.05, 0.1, 0.3, 1.0)
SEP_DELAY = 0.2


def wilson(k, n, z=1.96):
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return max(0.0, c - h), min(1.0, c + h)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("csv")
    ap.add_argument("--tilt-limit", type=float, default=20.0, dest="lim")
    ap.add_argument("--draws", type=int, default=200)
    ap.add_argument("--seed", type=int, default=3)
    ap.add_argument("--json")
    a = ap.parse_args()
    rows = list(csv.DictReader(open(a.csv)))
    if not rows or a.draws < 1:
        raise ValueError("Nonempty staging data and positive draw count required")
    rng = random.Random(a.seed)
    gyro_records=[];timer_records=[]
    coasts = sorted({float(r["coast"]) for r in rows})

    print(f"GYRO DRIFT at the ignition decision (tilt limit {a.lim:.0f} deg, "
          f"{a.draws} error draws per flight)")
    print(f"  {'coast s':>7s} {'bias deg/s':>10s} {'true>lim':>9s} {'false lockout':>22s} {'false pass':>22s}"
          f" {'err sd at t_dec':>15s}")
    for c in coasts:
        sub = [r for r in rows if float(r["coast"]) == c]
        n_true_over = sum(float(r["dec_tilt"]) > a.lim for r in sub)
        for b in BIASES:
            fl = fp = n_ok = n_over = 0
            sds = []
            for r in sub:
                T = float(r["dec_tilt"])
                t_dec = float(r.get("decision_time_s") or (float(r["t_bo1"]) + SEP_DELAY + c))
                s = b * t_dec
                sds.append(s)
                for _ in range(a.draws):
                    ex, ey = rng.gauss(0, s), rng.gauss(0, s)
                    m = math.hypot(T + ex, ey)
                    if T <= a.lim:
                        n_ok += 1
                        fl += m > a.lim
                    else:
                        n_over += 1
                        fp += m <= a.lim
            gyro_records.append(dict(coast=c,bias_dps=b,true_over=n_true_over,false_lockout_draws=fl,safe_draws=n_ok,base_within_tilt_limit=len(sub)-n_true_over,false_pass_draws=fp,unsafe_draws=n_over,false_lockout_pct=100*fl/n_ok if n_ok else None,false_pass_pct=100*fp/n_over if n_over else None,mean_error_sd_deg=st.mean(sds),base_flights=len(sub)))
            flp = f"{100*fl/n_ok:6.2f}%" if n_ok else "     -"
            fpp = f"{100*fp/n_over:6.2f}%" if n_over else "     -"
            print(f"  {c:7.0f} {b:10.2f} {n_true_over:9d} {flp:>22s} {fpp:>22s} {st.mean(sds):12.2f} deg")

    print("\nSEPARATION TIMER vs BURNOUT DETECTION (booster burnout time from the Monte Carlo)")
    by_sample={}
    for row in rows: by_sample.setdefault(row['sample'],set()).add(float(row['t_bo1']))
    for sample in by_sample:
        sample_coasts=[float(row['coast']) for row in rows if row['sample']==sample]
        if len(sample_coasts)!=len(coasts) or set(sample_coasts)!=set(coasts):
            raise ValueError('Each sample must occur exactly once at each paired coast')
    if any(len(values)!=1 for values in by_sample.values()):
        raise ValueError('Booster burnout times disagree between paired coasts for one sample')
    tb=sorted((sample,next(iter(values))) for sample,values in by_sample.items())
    t = sorted(x[1] for x in tb)
    n = len(t)
    print(f"  burnout time: n {n}, mean {st.mean(t):.3f} s, sd {st.stdev(t):.3f} s, "
          f"min {t[0]:.3f}, max {t[-1]:.3f}")
    nominal = st.median(t)
    for m in (0.0, 0.1, 0.2, 0.3, 0.5):
        T = nominal + m
        k = sum(x > T for x in t)
        lo, hi = wilson(k, n)
        timer_records.append(dict(margin_s=m,timer_s=T,early_count=k,flights=n,pct=100*k/n,wilson_low_pct=100*lo,wilson_high_pct=100*hi))
        print(f"  timer {T:.3f} s (median + {m:.1f} s): separation before burnout in {k}/{n} flights "
              f"= {100*k/n:.1f}% [{100*lo:.1f}-{100*hi:.1f}]")

    if a.json:
        Path(a.json).write_text(json.dumps(dict(gyro=gyro_records,timers=timer_records,burnout=dict(n=n,median_s=nominal,mean_s=st.mean(t),sd_s=st.stdev(t),min_s=min(t),max_s=max(t)),conditional_denominator=f"Draws conditional on truth-within or truth-over the tilt limit; underlying independent flights are n={n}, not gyro draws",compatibility_fields="safe_draws and unsafe_draws mean truth-within and truth-over the numerical tilt limit only; they are not physical safety classifications"),indent=2))


if __name__ == "__main__":
    main()
