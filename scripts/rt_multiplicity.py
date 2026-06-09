#!/usr/bin/env python3
"""Confirm the ~137 K transition is on the RT side (Step 3 diagnosis).

scripts/diagnose_transition.py showed the microphysics map T->haze is smooth
across 137 K, so the coupled loop's bistability must come from the radiative
transfer.  Here we test that directly and engine-independently: hold the haze
FIXED (the transition-state microphysics haze) and relax DISORT to radiative-
convective equilibrium from a WARM and a COOL initial temperature.  If they land
on different equilibria, the RT itself is bistable for that haze.

    .rtenv/bin/python scripts/rt_multiplicity.py
"""

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import numpy as np

from microphysics import Atmosphere, DEFAULT, solve_bvp_profile
from rt.column import Column
from rt.optics import OpticsParams
from rt.energy_balance import radiative_equilibrium


def scaled_atm(base, scale):
    z = base.z
    T = base.temperature(z)
    Tcold = float(T.min())
    above = z > z[int(np.argmin(T))]
    Tnew = T.copy()
    Tnew[above] = Tcold + scale * (T[above] - Tcold)
    return Atmosphere.from_profile(z, Tnew, P_surf=base.P_surf)


def main():
    base = Atmosphere.titan_reference()
    # FIXED haze: the microphysics haze at the transition (~137 K) state
    micro_fixed = solve_bvp_profile(scaled_atm(base, 0.66), DEFAULT, n_nodes=200)
    op = OpticsParams(haze_mode="rdg")

    def equilibrate(scale, label):
        atm = scaled_atm(base, scale)
        col = Column.from_atmosphere(atm, nlyr=100, p_top=1.0)
        eq, hist, _ = radiative_equilibrium(col, micro_fixed, op=op, n_iter=4000)
        print(f"  {label}: initial stratopause {col.T.max():.0f} K "
              f"-> equilibrium {eq.T.max():.1f} K @ {eq.P[np.argmax(eq.T)]:.1f} Pa "
              f"(resid {hist[-1]:.2f})")
        return eq

    print("Same fixed haze, two initial states (DISORT radiative-convective eq):")
    eqW = equilibrate(1.35, "warm start")
    eqC = equilibrate(0.60, "cool start")

    o = np.argsort(eqW.P)
    dT = np.abs(np.interp(eqW.P[o], eqC.P[np.argsort(eqC.P)],
                          eqC.T[np.argsort(eqC.P)]) - eqW.T[o])
    print(f"\n  |warm-eq - cool-eq| : max {dT.max():.1f} K, "
          f"stratopause {eqW.T.max():.0f} vs {eqC.T.max():.0f} K")
    if dT.max() > 10:
        print("  => RT IS BISTABLE for this haze: the two starts settle on "
              "different equilibria.\n     The ~137 K transition is a genuine "
              "radiative feedback (haze-altitude gain>1), not a microphysics or\n"
              "     mapping artifact.")
    else:
        print("  => DISORT is MONOSTABLE for this haze (both starts converge "
              "together).\n     The Fortran loop's jump is then path/cold-start "
              "dependence in that engine,\n     not intrinsic RT multiplicity -- "
              "still RT-side, not microphysics.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
