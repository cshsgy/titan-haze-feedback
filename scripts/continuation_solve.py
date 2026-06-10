#!/usr/bin/env python3
"""Damped branch-tracking solve of the coupled fixed point T = F(haze(T)).

The deferred Step-3 "continuation" solve.  Naive coupled iteration overshoots
(the composite map is steeply decreasing through the transition), so it hops
between attractors instead of converging.  Here we track ONE branch:

  T_{k+1} = (1-alpha) T_k + alpha * RT_eq( micro(T_k), init=T_k )

i.e. the haze is recomputed from the current iterate every pass (full coupling,
unlike the fixed-haze tests), the RT relaxation is *initialized from the current
iterate* (so when the RT is multivalued we stay on the branch continuously
connected to the start), and the update is under-relaxed by alpha.  A converged
chain satisfies T = F(haze(T)) to within `TOL_FP` -- a genuine coupled fixed
point.  Run from a warm and a cool start: two distinct converged fixed points =
the coupled system is bistable; the same fixed point = monostable.

    .rtenv/bin/python scripts/continuation_solve.py <mono|bi> <sigma_f> <warm|cool>

Writes writing/figs/continuation_<tag>.npz (T/P/z profile + stratopause history).
"""

import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import numpy as np

from microphysics import Atmosphere, DEFAULT, solve_bvp_profile
from microphysics.scaling_law_bimodal import solve_bimodal_kzero, to_micro
from rt.column import Column
from rt.optics import OpticsParams
from rt.energy_balance import radiative_equilibrium

NLYR = 100
N_ITER = 2500          # inner RT relaxation per outer pass
ALPHA = float(os.environ.get("CONT_ALPHA", "0.3"))   # outer under-relaxation
                       # (the khare-LW map was steep, |slope|>>1 -> damping; the
                       # obs-LW map is nearly flat, slope ~0.1 -> CONT_ALPHA=1)
TOL_FP = 1.0           # [K] fixed-point residual max|F(T)-T| to declare converged
P_FP_MIN = 3.0         # [Pa] exclude the top-two-layer 1 Pa zig-zag (known
                       # boundary artifact, +-7 K/day) from the fp residual
MAX_OUTER = 25
SIGMA_S = 1.5


def scaled_start(base, scale):
    z = base.z
    T = base.temperature(z)
    Tc = float(T.min())
    above = z > z[int(np.argmin(T))]
    Tn = T.copy()
    Tn[above] = Tc + scale * (T[above] - Tc)
    return Tn


def main():
    closure, sf, start = sys.argv[1], float(sys.argv[2]), sys.argv[3]
    tag = f"{closure}{sf:g}_{start}" if closure == "bi" else f"{closure}_{start}"

    base = Atmosphere.titan_reference()
    z = base.z
    op = OpticsParams(haze_mode="rdg")
    T = scaled_start(base, {"warm": 1.35, "cool": 0.60}[start])

    def micro(atm):
        if closure == "mono":
            return solve_bvp_profile(atm, DEFAULT, n_nodes=200)
        return to_micro(solve_bimodal_kzero(atm, DEFAULT, SIGMA_S, sf,
                                            n_out=200, n_quad=14))

    hist = []
    T_star, fp_resid = T, np.inf
    print(f"[{tag}] damped branch-tracking, alpha={ALPHA}, nlyr={NLYR}, "
          f"n_iter={N_ITER}", flush=True)
    for k in range(MAX_OUTER):
        atm_k = Atmosphere.from_profile(z, T, P_surf=base.P_surf)
        haze = micro(atm_k)
        col = Column.from_atmosphere(atm_k, nlyr=NLYR, p_top=1.0)   # RT init = T_k
        eq, h, _ = radiative_equilibrium(col, haze, op=op, n_iter=N_ITER)
        T_star = np.interp(z, eq.z, eq.T)                            # F(T_k) on atm grid
        mask = base.pressure(z) >= P_FP_MIN
        fp_resid = float(np.max(np.abs(T_star - T)[mask]))
        strat_in, strat_out = float(T.max()), float(T_star.max())
        hist.append((k, strat_in, strat_out, fp_resid, float(h[-1])))
        print(f"[{tag}] outer {k:2d}: strat {strat_in:6.1f} -> F(T) {strat_out:6.1f} K"
              f"   max|F(T)-T| {fp_resid:6.2f} K   rt_resid {h[-1]:.2f}", flush=True)
        if fp_resid < TOL_FP:
            print(f"[{tag}] CONVERGED coupled fixed point: stratopause "
                  f"{strat_out:.1f} K (fp resid {fp_resid:.2f} K)", flush=True)
            break
        T = (1.0 - ALPHA) * T + ALPHA * T_star
    else:
        print(f"[{tag}] NOT converged in {MAX_OUTER} outer iterations "
              f"(last fp resid {fp_resid:.2f} K)", flush=True)

    out = ROOT / "writing" / "figs" / f"continuation_{tag}.npz"
    np.savez(out, z=z, T=T, T_star=T_star, P=base.pressure(z),
             hist=np.array(hist), converged=float(fp_resid < TOL_FP))
    print(f"[{tag}] wrote {out.relative_to(ROOT)}", flush=True)


if __name__ == "__main__":
    main()
