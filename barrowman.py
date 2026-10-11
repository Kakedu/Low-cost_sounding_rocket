# Normal-force slope, CP and pitch damping from the fin/nose geometry (Barrowman),
# with a boundary-layer factor on the fins and the body crossflow coefficient.
# Writes barrowman_v6.json.

import json
import math

IN = 0.0254
D = 3.94 * IN
R = D / 2
L1 = 98.4 * IN
LN = 4.0 * D
FIN = dict(N=4, Cr=0.16, Ct=0.08, S=0.0426, Xt=0.08)
BOOST_FIN = dict(N=4, Cr=0.20, Ct=0.10, S=0.0, Xt=0.10)
STACK_CG_IN = 74.04
CLUSTER_FIN = dict(N=4, Cr=0.16, Ct=0.08, S=0.0, Xt=0.08)
L_BOOST = (196.8 - 98.4) * IN
MACHS = [0.0, 0.3, 0.5, 0.6, 0.8, 0.9, 1.0, 1.2]


NU_1KM, V_REP = 1.58e-5, 0.5 * 336.4
ETA_CF, CDC = 0.7, 1.2


def bl_thickness(x_m):
    re = V_REP * x_m / NU_1KM
    return 0.37 * x_m / re ** 0.2


def bl_factor(span, x_m):
    d = bl_thickness(x_m)
    if span <= d:
        return 7.0 / 9.0 * (span / d) ** (2.0 / 7.0)
    return (d * 7.0 / 9.0 + (span - d)) / span


def ogive_r(x, Ln, Rb):
    if x >= Ln:
        return Rb
    rho = (Rb * Rb + Ln * Ln) / (2 * Rb)
    return math.sqrt(rho ** 2 - (Ln - x) ** 2) + Rb - rho


def crossflow(length_m):
    n, A, M = 20000, 0.0, 0.0
    for i in range(n):
        x = (i + 0.5) * length_m / n
        w = 2 * ogive_r(x, LN, R) * length_m / n
        A += w
        M += w * x
    return dict(coef=ETA_CF * CDC * A / (math.pi * R * R), x_plan_in=M / A / IN,
                aplan_over_aref=A / (math.pi * R * R))


def fin_set(f, x_root_le):
    N, Cr, Ct, S, Xt = f["N"], f["Cr"], f["Ct"], f["S"], f["Xt"]
    Lf = math.sqrt(S ** 2 + (Xt + Ct / 2 - Cr / 2) ** 2)
    K = 1 + R / (S + R)
    cna0 = K * 4 * N * (S / D) ** 2 / (1 + math.sqrt(1 + (2 * Lf / (Cr + Ct)) ** 2))
    cna0 *= bl_factor(S, x_root_le)
    x = x_root_le + Xt / 3 * (Cr + 2 * Ct) / (Cr + Ct) + (Cr + Ct - Cr * Ct / (Cr + Ct)) / 6
    Af = 0.5 * (Cr + Ct) * S
    AR = 2 * S ** 2 / Af
    gam = math.atan((Xt + Ct / 2 - Cr / 2) / S)

    def dieder(M):
        b = math.sqrt(max(1e-6, abs(1 - M * M)))
        return AR / (2 + math.sqrt(4 + (b * AR / math.cos(gam)) ** 2))

    return cna0, x, (lambda M: dieder(min(M, 0.8)) / dieder(0.0))


def config(name, fin_sets, cg_launch_in, cg_burnout_in):
    cna_n, x_n = 2.0, 0.466 * LN
    fins = [fin_set(f, x) for f, x in fin_sets]
    rows = []
    for M in MACHS:
        parts = [(cna_n, x_n)] + [(c0 * fac(M), x) for c0, x, fac in fins]
        cna = sum(c for c, _ in parts)
        xcp = sum(c * x for c, x in parts) / cna
        rows.append(dict(mach=M, cna=cna, xcp_in=xcp / IN))
    out = dict(name=name, rows=rows)
    for tag, cg in (("launch", cg_launch_in), ("burnout", cg_burnout_in)):
        out[f"margin_{tag}_cal"] = (rows[0]["xcp_in"] - cg) / (D / IN)
    cg_mid = 0.5 * (cg_launch_in + cg_burnout_in) * IN
    parts0 = [(cna_n, x_n)] + [(c0, x) for c0, x, _ in fins]
    out["cmq"] = -2 * sum(c * ((x - cg_mid) / D) ** 2 for c, x in parts0)
    out["cg_launch_in"], out["cg_burnout_in"] = cg_launch_in, cg_burnout_in
    return out


if __name__ == "__main__":
    lo, hi = 0.01, 0.3
    for _ in range(50):
        FIN["S"] = (lo + hi) / 2
        c = config("single", [(FIN, L1 - FIN["Cr"])], 44.39, 38.47)
        lo, hi = (lo, FIN["S"]) if c["margin_launch_cal"] > 2.0 else (FIN["S"], hi)
    FIN["S"] = round(FIN["S"], 4)
    single = config("single", [(FIN, L1 - FIN["Cr"])], 44.39, 38.47)
    lo, hi = 0.01, 0.4
    for _ in range(50):
        BOOST_FIN["S"] = (lo + hi) / 2
        c = config("stack", [(FIN, L1 - FIN["Cr"]), (BOOST_FIN, L1 + L_BOOST - BOOST_FIN["Cr"])],
                   STACK_CG_IN, STACK_CG_IN)
        lo, hi = (lo, BOOST_FIN["S"]) if c["margin_launch_cal"] > 2.0 else (BOOST_FIN["S"], hi)
    BOOST_FIN["S"] = round(BOOST_FIN["S"], 4)
    stack = config("stack", [(FIN, L1 - FIN["Cr"]), (BOOST_FIN, L1 + L_BOOST - BOOST_FIN["Cr"])],
                   STACK_CG_IN, STACK_CG_IN)
    lo, hi = 0.01, 0.3
    for _ in range(50):
        CLUSTER_FIN["S"] = (lo + hi) / 2
        c = config("cluster", [(CLUSTER_FIN, L1 - CLUSTER_FIN["Cr"])], 52.16, 52.16)
        lo, hi = (lo, CLUSTER_FIN["S"]) if c["margin_launch_cal"] > 2.0 else (CLUSTER_FIN["S"], hi)
    CLUSTER_FIN["S"] = round(CLUSTER_FIN["S"], 4)
    cluster = config("cluster", [(CLUSTER_FIN, L1 - CLUSTER_FIN["Cr"])], 52.16, 45.07)
    res = dict(geometry=dict(D_m=D, L1_m=L1, nose="tangent ogive", Ln_m=LN,
                             fin=FIN, booster_fin=BOOST_FIN, cluster_fin=CLUSTER_FIN,
                             L_booster_m=L_BOOST),
               single=single, stack=stack, cluster=cluster,
               bl=dict(delta_single_mm=1e3 * bl_thickness(L1 - FIN["Cr"]),
                       factor_single=bl_factor(FIN["S"], L1 - FIN["Cr"]),
                       factor_booster=bl_factor(BOOST_FIN["S"], L1 + L_BOOST - BOOST_FIN["Cr"]),
                       factor_cluster=bl_factor(CLUSTER_FIN["S"], L1 - CLUSTER_FIN["Cr"])),
               crossflow=dict(single=crossflow(L1), stack=crossflow(L1 + L_BOOST)))
    for c in (single, stack, cluster):
        print(f"== {c['name']}: margin launch {c['margin_launch_cal']:.2f} cal, "
              f"burnout {c['margin_burnout_cal']:.2f} cal, Cmq {c['cmq']:.0f}")
        for r in c["rows"]:
            print(f"   M {r['mach']:.1f}  CNa {r['cna']:.2f}  Xcp {r['xcp_in']:.2f} in")
    print("fin spans (m): single", FIN["S"], "booster", BOOST_FIN["S"], "cluster", CLUSTER_FIN["S"])
    print("boundary layer", {k: round(v, 3) for k, v in res["bl"].items()})
    print("crossflow", res["crossflow"])
    json.dump(res, open("barrowman_v6.json", "w"), indent=1)
