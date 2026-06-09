#!/usr/bin/env python3
"""Convergence check: is the bimodal transition-state warm/cool split real or
under-iterated?  Re-run the sigma_F=2.0 bimodal equilibria at increasing n_iter."""

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import numpy as np
from microphysics import Atmosphere, DEFAULT
from microphysics.scaling_law_bimodal import solve_bimodal_kzero, to_micro
from rt.column import Column
from rt.optics import OpticsParams
from rt.energy_balance import radiative_equilibrium


def scaled_atm(base, scale):
    z = base.z; T = base.temperature(z); Tc = float(T.min())
    above = z > z[int(np.argmin(T))]
    Tn = T.copy(); Tn[above] = Tc + scale * (T[above] - Tc)
    return Atmosphere.from_profile(z, Tn, P_surf=base.P_surf)


base = Atmosphere.titan_reference()
haze = to_micro(solve_bimodal_kzero(scaled_atm(base, 0.66), DEFAULT, 1.5, 2.0,
                                    n_out=200, n_quad=14))
op = OpticsParams(haze_mode="rdg")

for ni in (4000, 8000, 16000):
    res = {}
    for label, scale in (("warm", 1.35), ("cool", 0.60)):
        col = Column.from_atmosphere(scaled_atm(base, scale), nlyr=100, p_top=1.0)
        eq, h, _ = radiative_equilibrium(col, haze, op=op, n_iter=ni)
        res[label] = (eq.T.max(), float(h[-1]))
    dT = abs(res["warm"][0] - res["cool"][0])
    print(f"n_iter={ni:6d}: warm {res['warm'][0]:.1f} K (resid {res['warm'][1]:.3f})  "
          f"cool {res['cool'][0]:.1f} K (resid {res['cool'][1]:.3f})  split {dT:.1f} K")
