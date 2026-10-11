# Static FE check of one brake flap (scikit-fem), plus its first modes.

import json
import math

import numpy as np
from scipy.sparse.linalg import eigsh
from skfem import (Basis, ElementHex2, ElementVector, FacetBasis, MeshHex, asm, condense, solve,
                   BilinearForm, LinearForm)
from skfem.helpers import dot, ddot, sym_grad, eye, trace
from skfem.models.elasticity import lame_parameters, linear_elasticity

B, L, T = 0.060, 0.055, 0.003
E, NU, RHO, SY = 68.9e9, 0.33, 2700.0, 276e6
F_FLAP = 186.0
import os as _os
if _os.path.exists("results_v6/loads_v6.json"):
    _c = json.load(open("results_v6/loads_v6.json"))["cases"]
    F_FLAP = max(c["brake_per_flap_N"] for c in _c if c["case"].startswith("impulse"))
P = F_FLAP / (B * L)
LAM, MU = lame_parameters(E, NU)


def run(nx, ny, nz, modes=False):
    m = MeshHex.init_tensor(np.linspace(0, B, nx + 1), np.linspace(0, L, ny + 1),
                            np.linspace(0, T, nz + 1)).with_defaults()
    e = ElementVector(ElementHex2())
    ib = Basis(m, e, intorder=4)
    K = asm(linear_elasticity(LAM, MU), ib)
    top = m.facets_satisfying(lambda x: np.isclose(x[2], T))
    fb = FacetBasis(m, e, facets=top, intorder=4)

    @LinearForm
    def pressure(v, w):
        return -P * v[2]

    f = asm(pressure, fb)
    clamped = ib.get_dofs(lambda x: np.isclose(x[1], 0.0)).all()
    u = solve(*condense(K, f, D=clamped))
    w_max = -u[ib.nodal_dofs[2]].min()
    du = ib.interpolate(u)
    eps = sym_grad(du)
    sig = 2 * MU * eps + LAM * eye(trace(eps), 3)
    s = sig
    vm = np.sqrt(0.5 * ((s[0, 0] - s[1, 1]) ** 2 + (s[1, 1] - s[2, 2]) ** 2 + (s[2, 2] - s[0, 0]) ** 2)
                 + 3 * (s[0, 1] ** 2 + s[1, 2] ** 2 + s[0, 2] ** 2))
    y_q = ib.global_coordinates().value[1]
    away = y_q > 2 * T
    out = dict(mesh=f"{nx}x{ny}x{nz}", dofs=int(K.shape[0]), w_max_mm=float(1e3 * w_max),
               vm_max_MPa=float(vm.max()) / 1e6, vm_max_away_MPa=float(vm[away].max()) / 1e6)
    if modes:
        @BilinearForm
        def mass(u_, v_, w):
            return RHO * dot(u_, v_)
        M = asm(mass, ib)
        I = ib.complement_dofs(clamped)
        vals, _ = eigsh(K[I][:, I], k=3, M=M[I][:, I], sigma=0.0)
        out["f_Hz"] = [float(math.sqrt(v) / (2 * math.pi)) for v in sorted(vals)]
    return out


ST_PLATE = 0.146
FACTORS = dict(dcd_2sd=1.0 + 2 * 0.20, ultimate=1.5)


def dynamics(f1_plate):
    I_h = RHO * B * T * L ** 3 / 3
    hinge = {f"k={k:g} N.m/rad": math.sqrt(k / I_h) / (2 * math.pi) for k in (1, 10, 100, 1000)}
    k_req = I_h * (2 * math.pi * 2 * f1_plate) ** 2
    band = {}
    for v in (100.0, 200.0, 270.0):
        band[f"V={v:.0f} m/s"] = {f"h={h * 1e3:.0f} mm": ST_PLATE * v / h for h in (0.010, 0.025, 0.055)}
    return dict(I_hinge_kgm2=I_h, hinge_mode_Hz=hinge, k_for_2x_f1_Nm_per_rad=k_req, shedding_Hz=band)


def beam_theory():
    I = B * T ** 3 / 12
    M_root = F_FLAP * L / 2
    return dict(sigma_root_MPa=M_root * (T / 2) / I / 1e6,
                w_tip_mm=1e3 * F_FLAP * L ** 3 / (8 * E * I),
                f1_Hz=(1.8751 ** 2 / (2 * math.pi)) * math.sqrt(E * I / (RHO * B * T * L ** 4)))


if __name__ == "__main__":
    print(f"pressure {P / 1e3:.1f} kPa on {B * 1e3:.0f} x {L * 1e3:.0f} x {T * 1e3:.0f} mm plate")
    res = [run(6, 6, 1), run(12, 11, 2), run(24, 22, 3), run(30, 28, 4, modes=True)]
    for r in res:
        print(r)
    bt = beam_theory()
    print("beam theory", {k: round(v, 2) for k, v in bt.items()})
    fine = res[-1]
    summary = dict(geometry_mm=[B * 1e3, L * 1e3, T * 1e3], material="6061-T6", load_N=F_FLAP,
                   pressure_kPa=P / 1e3, meshes=res, beam_theory=bt,
                   sf_yield_beam=SY / 1e6 / bt["sigma_root_MPa"],
                   sf_yield_fe_away=SY / 1e6 / fine["vm_max_away_MPa"])
    print(f"safety factor on yield: beam-theory root stress {summary['sf_yield_beam']:.2f}, "
          f"FE away from the clamp corner {summary['sf_yield_fe_away']:.2f}")
    dyn = dynamics(fine["f_Hz"][0])
    lf = FACTORS["dcd_2sd"] * FACTORS["ultimate"]
    summary.update(dynamics=dyn, load_factor=lf, factors=FACTORS,
                   sf_ultimate_beam=summary["sf_yield_beam"] / lf,
                   margin_of_safety_beam=summary["sf_yield_beam"] / lf - 1)
    print("hinge (rigid-body) mode, Hz:", {k: round(v) for k, v in dyn["hinge_mode_Hz"].items()},
          f" stiffness for a hinge mode above 2 x f1: {dyn['k_for_2x_f1_Nm_per_rad']:.0f} N.m/rad")
    print("shedding frequency St V / h, Hz:", {k: {kk: round(vv) for kk, vv in v.items()}
                                               for k, v in dyn["shedding_Hz"].items()})
    print(f"load factor {lf:.2f} (dCd 2 sd x ultimate 1.5): yield SF at factored load "
          f"{summary['sf_ultimate_beam']:.2f}, margin of safety {summary['margin_of_safety_beam']:+.2f}")
    json.dump(summary, open("results_v6/flap_fea_v6.json", "w"), indent=1)
