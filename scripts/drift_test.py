#!/usr/bin/env python3
"""Stiffness check for the continuation solve: relax the WARM-tracked coupled
fixed point with a long, tight-tolerance RT run (tol 0.3 -> 0.05 K/day) and see
whether the radiatively-sluggish lower stratosphere (50-200 km), which holds
10-18 K of initial-condition memory at the standard tolerance, drifts toward the
cool-tracked fixed point (artifact) or stays put (real multiplicity).

    .rtenv/bin/python scripts/drift_test.py <mono|bi2>
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

pre = sys.argv[1]
F = ROOT / "writing" / "figs"
w = np.load(F / f"continuation_{pre}_warm.npz")
c = np.load(F / f"continuation_{pre}_cool.npz")
base = Atmosphere.titan_reference()
atm = Atmosphere.from_profile(base.z, w["T"], P_surf=base.P_surf)
haze = (solve_bvp_profile(atm, DEFAULT, n_nodes=200) if pre == "mono"
        else to_micro(solve_bimodal_kzero(atm, DEFAULT, 1.5, 2.0,
                                          n_out=200, n_quad=14)))
col = Column.from_atmosphere(atm, nlyr=100, p_top=1.0)
eq, h, _ = radiative_equilibrium(col, haze, op=OpticsParams(haze_mode="rdg"),
                                 n_iter=20000, tol=0.05)
Tn = np.interp(base.z, eq.z, eq.T)
print(f"[{pre}] long tight relax from WARM fp (resid {h[-1]:.3f}, {len(h)} iters):",
      flush=True)
for zkm in (76, 90, 120):
    i = int(np.argmin(np.abs(base.z / 1e3 - zkm)))
    print(f"  z={zkm:3d} km: warm fp {w['T'][i]:6.1f} -> {Tn[i]:6.1f}  "
          f"(cool fp {c['T'][i]:6.1f})", flush=True)
band = (base.z / 1e3 > 50) & (base.z / 1e3 < 200)
d0 = np.abs(w["T"] - c["T"])[band].max()
d1 = np.abs(Tn - c["T"])[band].max()
print(f"  50-200 km band: max|warm-cool| {d0:.1f} -> {d1:.1f} K after tight relax",
      flush=True)
print(f"  stratopause: {w['T'].max():.1f} -> {Tn.max():.1f} "
      f"(cool fp {c['T'].max():.1f})", flush=True)
