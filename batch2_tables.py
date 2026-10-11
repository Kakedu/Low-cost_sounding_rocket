# Builds the Monte Carlo / SIL summary tables (hit rates, error budget, McNemar tests).

import argparse
import csv
import math
import os
import random
import statistics as st

TOLS = (50.0, 75.0, 100.0)
TOL = 100.0
TARGET = 10000.0
D = "."


def load(name):
    p = os.path.join(D, name)
    if not os.path.exists(p):
        return None
    return {int(r["sample"]): r for r in csv.DictReader(open(p))}


def wilson(k, n, z=1.96):
    if n == 0:
        return (float("nan"),) * 2
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return c - h, c + h


METRIC = "geo"


def ap(r):
    a = float(r["apogee"])
    if METRIC == "baro":
        a *= 288.15 / (288.15 + float(r.get("dT_K", 0.0) or 0.0))
    return a


def stats(rows, target=TARGET, tol=TOL):
    a = [ap(r) for r in rows.values()]
    n = len(a)
    k_in = sum(abs(x - target) <= tol for x in a)
    k_lo = sum(x < target - tol for x in a)
    k_hi = sum(x > target + tol for x in a)
    return dict(n=n, mean=st.mean(a), sd=st.stdev(a), p_in=k_in / n, ci_in=wilson(k_in, n),
                p_lo=k_lo / n, ci_lo=wilson(k_lo, n), p_hi=k_hi / n, ci_hi=wilson(k_hi, n))


def fmt(s):
    return (f"n {s['n']:4d}  mean {s['mean']:7.0f}  sd {s['sd']:5.0f}  "
            f"hit {100*s['p_in']:5.1f}% [{100*s['ci_in'][0]:.1f}-{100*s['ci_in'][1]:.1f}]  "
            f"low {100*s['p_lo']:5.1f}% [{100*s['ci_lo'][0]:.1f}-{100*s['ci_lo'][1]:.1f}]  "
            f"high {100*s['p_hi']:5.1f}% [{100*s['ci_hi'][0]:.1f}-{100*s['ci_hi'][1]:.1f}]")


def common(*sets):
    keys = set(sets[0])
    for s in sets[1:]:
        keys &= set(s)
    return sorted(keys)


def mcnemar(a, b, target=TARGET, tol=TOL):
    ks = common(a, b)
    ha = [abs(ap(a[k]) - target) <= tol for k in ks]
    hb = [abs(ap(b[k]) - target) <= tol for k in ks]
    b10 = sum(x and not y for x, y in zip(ha, hb))
    b01 = sum(y and not x for x, y in zip(ha, hb))
    n = b10 + b01
    if n == 0:
        return b10, b01, 1.0
    k = min(b10, b01)
    p = sum(math.comb(n, i) for i in range(k + 1)) / 2 ** n
    return b10, b01, min(1.0, 2 * p)


def paired_abs_diff(a, b, target=TARGET, n_boot=2000, seed=7):
    ks = common(a, b)
    d = [abs(ap(a[k]) - target) - abs(ap(b[k]) - target) for k in ks]
    rng = random.Random(seed)
    bs = sorted(st.median(rng.choice(d) for _ in d) for _ in range(n_boot))
    return st.median(d), bs[int(0.025 * n_boot)], bs[int(0.975 * n_boot) - 1]


def budget(o, cases, tol=TOL):
    print(f"\nERROR BUDGET (misses of +/-{tol:.0f} ft; energy short = open loop below "
          f"{TARGET - tol:.0f} ft)")
    for name, rows in cases.items():
        ks = common(o, rows)
        short = {k: ap(o[k]) < TARGET - tol for k in ks}
        miss = {k: abs(ap(rows[k]) - TARGET) > tol for k in ks}
        e = sum(miss[k] and short[k] for k in ks)
        lo = sum(ap(rows[k]) < TARGET - tol and not short[k] for k in ks)
        hi = sum(ap(rows[k]) > TARGET + tol for k in ks)
        reach = [k for k in ks if not short[k]]
        hit_r = sum(not miss[k] for k in reach)
        print(f"  {name:34s} n {len(ks)}  misses {sum(miss.values()):3d} = energy short {e:3d}"
              f" + other low {lo:3d} + high {hi:3d}   hit when reachable "
              f"{100 * hit_r / len(reach):5.1f}% [{100*wilson(hit_r, len(reach))[0]:.1f}-"
              f"{100*wilson(hit_r, len(reach))[1]:.1f}]")


def src(rows, inputs):
    import numpy as np
    ks = sorted(rows)
    X = np.array([[float(rows[k][i]) for i in inputs] for k in ks])
    y = np.array([ap(rows[k]) for k in ks])
    keep = X.std(axis=0) > 0
    X, names = X[:, keep], [n for n, k in zip(inputs, keep) if k]
    Z = (X - X.mean(0)) / X.std(0)
    yz = (y - y.mean()) / y.std()
    A = np.column_stack([Z, np.ones(len(yz))])
    b, *_ = np.linalg.lstsq(A, yz, rcond=None)
    r2 = 1 - ((yz - A @ b) ** 2).sum() / (yz ** 2).sum()
    return sorted(zip(names, b[:-1]), key=lambda z: -abs(z[1])), r2


def main():
    global D
    ap_ = argparse.ArgumentParser()
    ap_.add_argument("--dir", default=".")
    ap_.add_argument("--open-nominal", type=float, default=10554.5, dest="open_nominal",
                     help="nominal open-loop apogee, ft (sweep.py nominal --open)")
    a = ap_.parse_args()
    D = a.dir
    o = load("mc_open_v6.csv")
    if o is None:
        raise SystemExit("mc_open_v6.csv not found in " + D)
    cases = {}
    for tag, f in (("XML, true state, nominal Cd", "mc_closed_v6.csv"),
                   ("XML, true state, oracle Cd", "mc_oracle_v6.csv"),
                   ("SIL KF3 ideal sensors, est Cd", "mc_sil_kf_est_ideal_v6.csv"),
                   ("SIL KF3, est Cd", "mc_sil_kf_est_v6.csv"),
                   ("SIL KF3, nominal Cd", "mc_sil_kf_nom_v6.csv"),
                   ("SIL KF4 ideal sensors, est Cd", "mc_sil_kf4_est_ideal_v6.csv"),
                   ("SIL KF4, est Cd", "mc_sil_kf4_est_v6.csv"),
                   ("SIL KF4, nominal Cd", "mc_sil_kf4_nom_v6.csv"),
                   ("SIL KF4, oracle Cd", "mc_sil_kf4_orc_v6.csv"),
                   ("SIL KF4, est Cd, temp known", "mc_sil_kf4_est_notemp_v6.csv"),
                   ("SIL KF4, est Cd, FC assumes ISA", "mc_sil_kf4_est_isa_v6.csv")):
        r = load(f)
        if r:
            cases[tag] = r

    global METRIC
    for METRIC in ("geo", "baro"):
        print(f"\n################ METRIC: {'geometric apogee' if METRIC == 'geo' else 'barometric altimeter reading (ISA)'}")
        report(o, cases, a)
    METRIC = "geo"
    rest(o, cases, a)


def report(o, cases, a):
    print("APOGEE STATISTICS, target 10,000 ft, +/-100 ft")
    print(f"  {'open loop':34s} {fmt(stats(o))}")
    for name, rows in cases.items():
        print(f"  {name:34s} {fmt(stats(rows))}")
        r0 = next(iter(rows.values()))
        if "h_rms_m" in r0:
            g = lambda k: st.mean(float(r[k]) for r in rows.values() if r[k] not in ("", "nan"))
            print(f"  {'':34s} mean rms errors: h {g('h_rms_m'):.2f} m, v {g('v_rms_ms'):.2f} m/s, "
                  f"Cd {g('cd_rms'):.3f}")

    print("\nTOLERANCE (P(hit) at +/-50 / 75 / 100 ft)")
    for name, rows in [("open loop", o)] + list(cases.items()):
        print(f"  {name:34s} " + "   ".join(
            f"{100*stats(rows, tol=t)['p_in']:5.1f}%" for t in TOLS))

    for t in TOLS:
        budget(o, cases, t)

    print("\nPAIRED COMPARISONS (McNemar exact, +/-100 ft; median paired difference of"
          " |error|, ft, negative = first is better)")
    pairs = [("SIL KF4, est Cd", "SIL KF3, est Cd"),
             ("SIL KF4, est Cd", "SIL KF4, nominal Cd"),
             ("SIL KF4, oracle Cd", "SIL KF4, est Cd"),
             ("SIL KF4, est Cd, temp known", "SIL KF4, est Cd"),
             ("SIL KF4, est Cd, FC assumes ISA", "SIL KF4, est Cd"),
             ("SIL KF4 ideal sensors, est Cd", "SIL KF4, est Cd, temp known"),
             ("SIL KF3 ideal sensors, est Cd", "XML, true state, oracle Cd"),
             ("SIL KF4 ideal sensors, est Cd", "SIL KF3 ideal sensors, est Cd"),
             ("XML, true state, oracle Cd", "XML, true state, nominal Cd")]
    for p1, p2 in pairs:
        if p1 in cases and p2 in cases:
            b10, b01, p = mcnemar(cases[p1], cases[p2])
            md, lo, hi = paired_abs_diff(cases[p1], cases[p2])
            print(f"  {p1} vs {p2}\n      hit only in first {b10}, only in second {b01}, p = {p:.3g};"
                  f"  median d|err| {md:+.0f} ft [{lo:+.0f}, {hi:+.0f}]")

    print("\nSENSOR ERRORS ONE AT A TIME (SIL 3-state KF, estimated Cd; flights 0-299)")
    ideal = cases.get("SIL KF3 ideal sensors, est Cd")
    for e in ("port", "scale", "quant", "lat", "temp"):
        r = load(f"mc_sil_kf_est_only_{e}_v6.csv")
        if r:
            line = f"  only {e:6s} {fmt(stats(r))}"
            if ideal:
                b10, b01, p = mcnemar(r, ideal)
                line += f"  vs ideal: hit only here {b10}, only ideal {b01}, p {p:.3g}"
            print(line)
    for lab, f in (("all, q-aware KF", "mc_sil_kf_est_qaware_v6.csv"),):
        r = load(f)
        if r:
            print(f"  {lab:11s} {fmt(stats(r))}")


def rest(o, cases, a):
    print("\nMARGIN CURVE (margin = nominal open-loop apogee - target)")
    for lab, pat, base in (("XML, nominal Cd", "mc_closed_t{t}_v6.csv", "mc_closed_v6.csv"),
                           ("SIL KF3 ideal sensors, est Cd", "mc_sil_kf_est_ideal_t{t}_v6.csv",
                            "mc_sil_kf_est_ideal_v6.csv"),
                           ("SIL KF3, est Cd", "mc_sil_kf_est_t{t}_v6.csv", "mc_sil_kf_est_v6.csv"),
                           ("SIL KF4, est Cd", "mc_sil_kf4_est_t{t}_v6.csv", "mc_sil_kf4_est_v6.csv"),
                           ("SIL KF4, est Cd, temp known", "mc_sil_kf4_est_notemp_t{t}_v6.csv",
                            "mc_sil_kf4_est_notemp_v6.csv")):
        print(f"  {lab}")
        for t in (9600, 9800, 10000, 10200, 10400):
            rows = load(base if t == 10000 else pat.format(t=t))
            if rows:
                m = a.open_nominal - t
                print(f"  {t:7d} {m:7.0f} ft {100*m/t:5.1f}%  {fmt(stats(rows, t))}")

    print("\nINPUT SPREAD (drag CV and impulse sd)")
    for lab, fo, fc, fk in (("Cd CV 0.10", "mc_open_s10_v6.csv", "mc_closed_s10_v6.csv", "mc_sil_kf4_est_s10_v6.csv"),
                            ("Cd CV 0.20", "mc_open_v6.csv", "mc_closed_v6.csv", "mc_sil_kf4_est_v6.csv"),
                            ("Cd CV 0.30", "mc_open_s30_v6.csv", "mc_closed_s30_v6.csv", "mc_sil_kf4_est_s30_v6.csv"),
                            ("impulse sd 0.01", "mc_open_i01_v6.csv", "mc_closed_i01_v6.csv", None),
                            ("impulse sd 0.05", "mc_open_i05_v6.csv", "mc_closed_i05_v6.csv", None)):
        for k, f in (("open", fo), ("XML closed", fc), ("SIL KF4 est", fk)):
            r = load(f) if f else None
            if r:
                print(f"  {lab:16s} {k:10s} {fmt(stats(r))}")

    for lab in ("open", "closed"):
        p = os.path.join(D, f"morris_{lab}_v6_summary.csv")
        if not os.path.exists(p):
            continue
        print(f"\nMORRIS, {lab} loop (ft per unit of the 2.5-97.5 % quantile range)")
        for r in sorted(csv.DictReader(open(p)), key=lambda z: -float(z["mu_star"])):
            mus = float(r["mu_star"])
            print(f"  {r['input']:12s} mu* {mus:7.0f} [{float(r['ci_lo']):6.0f}, {float(r['ci_hi']):6.0f}]"
                  f"  mu {float(r['mu']):7.0f}  sigma {float(r['sigma']):7.0f}")

    inputs = ["prop_mass", "cd0", "cn_alpha", "cp_x", "dry_mass", "launch_elev", "wind_fps",
              "airbrake_dcd", "dT_K"]
    for name, rows in (("open loop", o), ("XML closed", cases.get("XML, true state, nominal Cd")),
                       ("SIL KF4 est", cases.get("SIL KF4, est Cd"))):
        if rows is None:
            continue
        r0 = next(iter(rows.values()))
        if "wind_fps" not in r0:
            continue
        res, r2 = src(rows, [i for i in inputs if i in r0])
        print(f"\nSRC, {name} (R^2 of the linear fit {r2:.3f})")
        for n, b in res:
            print(f"  {n:12s} {b:+.3f}")


if __name__ == "__main__":
    main()
