#!/usr/bin/env python3
"""Phase 1 measurement: stratopause effect of the obs-calibrated haze LW opacity.

DISORT radiative-convective equilibrium with the SAME mono BVP haze, switching
only the LW haze model: 'khare' (RDG with Khare lab tholin k, the old behaviour,
20-90x the observational LW absorptivity) vs 'obs' (chi(wn,P) from the
observational preschaze tables).  Target: the prescribed-haze stratopause ~183 K.

    .rtenv/bin/python scripts/lw_fix_effect.py
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

atm = Atmosphere.titan_reference()
mono = solve_bvp_profile(atm, DEFAULT, n_nodes=200)
out = {}
for lw in ("khare", "obs"):
    col = Column.from_atmosphere(atm, nlyr=100, p_top=1.0)
    eq, h, _ = radiative_equilibrium(col, mono, op=OpticsParams(haze_mode="rdg",
                                                                lw_haze=lw),
                                     n_iter=4000)
    out[lw] = eq
    print(f"lw_haze={lw:6s}: stratopause {eq.T.max():.1f} K "
          f"@ {eq.P[np.argmax(eq.T)]:.1f} Pa (resid {h[-1]:.2f})", flush=True)
print(f"\nprescribed-haze reference: ~183-200 K; "
      f"LW fix moves the stratopause {out['obs'].T.max()-out['khare'].T.max():+.1f} K")
