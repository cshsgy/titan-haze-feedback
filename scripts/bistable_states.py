#!/usr/bin/env python3
"""The converged equilibria of the coupled system: bistable (monodisperse) vs
monostable (polydisperse).

Replaces the iteration-history view.  Instead of plotting the non-converging
coupled loop, we show the *stable states each model settles to*: hold the haze
fixed (the transition-state haze of each model) and relax DISORT to radiative-
convective equilibrium from a WARM and a COOL initial profile.

  - monodisperse haze: the two starts land on two distinct equilibria (bistable);
  - bimodal (polydisperse) haze: both starts converge to one equilibrium
    (monostable), even at the thick, nearly-observed sigma_F=1.2 limit.

    .rtenv/bin/python scripts/bistable_states.py
Writes writing/figs/bistable_states.png
"""

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
N_ITER = 4000
WARM, COOL = 1.35, 0.60


def scaled_atm(base, scale):
    z = base.z
    T = base.temperature(z)
    Tcold = float(T.min())
    above = z > z[int(np.argmin(T))]
    Tnew = T.copy()
    Tnew[above] = Tcold + scale * (T[above] - Tcold)
    return Atmosphere.from_profile(z, Tnew, P_surf=base.P_surf)


def equilibrate(base, micro, scale):
    atm = scaled_atm(base, scale)
    col = Column.from_atmosphere(atm, nlyr=NLYR, p_top=1.0)
    eq, hist, _ = radiative_equilibrium(col, micro, op=OpticsParams(haze_mode="rdg"),
                                        n_iter=N_ITER)
    return eq, float(hist[-1])


def max_profile_split(eqW, eqC):
    """Largest |warm-eq - cool-eq| over the shared pressure grid."""
    oW, oC = np.argsort(eqW.P), np.argsort(eqC.P)
    dT = np.abs(np.interp(eqW.P[oW], eqC.P[oC], eqC.T[oC]) - eqW.T[oW])
    return float(dT.max())


def main():
    base = Atmosphere.titan_reference()
    # transition-state haze for each model (input T scaled to the ~137 K state),
    # the same fixed haze used in the monodisperse RT-multiplicity diagnosis.
    atm_t = scaled_atm(base, 0.66)
    hazes = {
        "monodisperse": solve_bvp_profile(atm_t, DEFAULT, n_nodes=200),
        "bimodal_sf1.2": to_micro(solve_bimodal_kzero(atm_t, DEFAULT, 1.5, 1.2,
                                                      n_out=200, n_quad=14)),
        "bimodal_sf2.0": to_micro(solve_bimodal_kzero(atm_t, DEFAULT, 1.5, 2.0,
                                                      n_out=200, n_quad=14)),
    }

    cache = {}
    for name, haze in hazes.items():
        eqW, rW = equilibrate(base, haze, WARM)
        eqC, rC = equilibrate(base, haze, COOL)
        dT = abs(eqW.T.max() - eqC.T.max())
        dTmax = max_profile_split(eqW, eqC)
        cache[f"{name}_W_T"] = eqW.T; cache[f"{name}_W_P"] = eqW.P
        cache[f"{name}_C_T"] = eqC.T; cache[f"{name}_C_P"] = eqC.P
        print(f"{name:14s}: warm {eqW.T.max():.1f} K (resid {rW:.2f})  "
              f"cool {eqC.T.max():.1f} K (resid {rC:.2f})  -> stratopause split "
              f"{dT:.0f} K, max-profile split {dTmax:.0f} K  "
              f"{'BISTABLE' if dT > 10 else 'monostable'}")
    cache["base_T"] = base.temperature(base.z); cache["base_P"] = base.pressure(base.z)
    p = ROOT / "writing" / "figs" / "bistable_states.npz"
    np.savez(p, **cache)
    print(f"\nwrote {p.relative_to(ROOT)}  (plot with scripts/plot_bistable.py)")


if __name__ == "__main__":
    main()
