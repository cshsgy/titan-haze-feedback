#!/usr/bin/env python3
"""Effect of Coulomb coagulation inhibition (n_e=15 e-/um, use_charge=True) on
the haze profile, its tau(P) match to the observational haze, and the DISORT
equilibrium (obs-calibrated LW).

    .rtenv/bin/python scripts/charge_effect.py
"""

import sys
from dataclasses import replace
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import numpy as np

from microphysics import Atmosphere, DEFAULT, solve_bvp_profile
from microphysics.scaling_law_bimodal import solve_bimodal_kzero, to_micro
from coupling.presc_haze import observational_haze
from coupling.presc_haze import microphysics_haze
from rt.column import Column
from rt.optics import OpticsParams
from rt.energy_balance import radiative_equilibrium

atm = Atmosphere.titan_reference()
obs_sw, _ = observational_haze()
jv = int(np.argmin(np.abs(obs_sw.wn - 20000.0)))
PC = replace(DEFAULT, use_charge=True)

cases = {
    "mono":         solve_bvp_profile(atm, DEFAULT, n_nodes=200),
    "mono+q":       solve_bvp_profile(atm, PC, n_nodes=200),
    "bi1.2":        to_micro(solve_bimodal_kzero(atm, DEFAULT, 1.5, 1.2, n_out=200, n_quad=14)),
    "bi1.2+q":      to_micro(solve_bimodal_kzero(atm, PC, 1.5, 1.2, n_out=200, n_quad=14)),
}

print("=== cumulative visible tau from TOA (obs in parens) ===")
taus = {k: microphysics_haze(m, atm)[0] for k, m in cases.items()}
hdr = "P [Pa]    obs    " + "  ".join(f"{k:>9s}" for k in cases)
print(hdr)
for P0 in (1.0, 10.0, 100.0, 1000.0, 1e4, 1.4e5):
    i = int(np.argmin(np.abs(obs_sw.pl - P0)))
    row = "  ".join(f"{taus[k].tau[i, jv]:9.3f}" for k in cases)
    print(f"{obs_sw.pl[i]:9.1f} {obs_sw.tau[i, jv]:6.3f}  {row}")

print("\n=== Nbar at selected P (charge keeps particles smaller) ===", flush=True)
for k in ("mono", "mono+q"):
    m = cases[k]
    P = atm.pressure(m.z)
    vals = [f"{m.Nbar[int(np.argmin(np.abs(P-P0)))]:9.0f}" for P0 in (10., 100., 1000.)]
    print(f"  {k:8s} Nbar @ 10/100/1000 Pa: " + " ".join(vals))

print("\n=== DISORT equilibrium (obs LW), charge off vs on ===", flush=True)
for k in ("mono", "mono+q", "bi1.2", "bi1.2+q"):
    col = Column.from_atmosphere(atm, nlyr=100, p_top=1.0)
    eq, h, _ = radiative_equilibrium(col, cases[k], op=OpticsParams(haze_mode="rdg"),
                                     n_iter=6000)
    print(f"  {k:8s}: stratopause {eq.T.max():.1f} K @ {eq.P[np.argmax(eq.T)]:.1f} Pa "
          f"(resid {h[-1]:.2f})", flush=True)
print("\nreference: prescribed-haze ~183 K; uncharged mono one-shot 194.2 K")
