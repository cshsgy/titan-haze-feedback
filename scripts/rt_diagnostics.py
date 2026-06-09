#!/usr/bin/env python3
"""Find discrepancies between our DISORT RT and the reference Fortran RT.

Compares the two engines at the SAME state: takes the Fortran's converged
pressure/temperature profile and (prescribed) haze, computes our shortwave and
longwave heating + fluxes on it, and compares to the Fortran's own diagnostics
at that snapshot.  Plus targeted checks for the methodological differences the
code analysis surfaced (Rayleigh scattering; CIA exp-sum fit vs HITRAN).

    .rtenv/bin/python scripts/rt_diagnostics.py
"""

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import numpy as np

from microphysics import Atmosphere, DEFAULT, solve_bvp_profile, constants as C
from microphysics.atmosphere import gravity
from rt.column import Column, TITAN_DAY
from rt.optics import OpticsParams
from rt.correlated_k import CorrelatedKSW, CorrelatedKLW
from rt.cia import CIABands, Composition
from rt import energy_balance as eb

FORT = ROOT / "src" / "example_bowen_fort"
DATA = FORT / "INPUT" / "DATA"


def fortran_state():
    """Last valid (non-NaN) T snapshot + its sw/lw heating [K/s] per layer."""
    T = np.loadtxt(FORT / "temperatures.txt")
    P = T[0]
    sw = np.loadtxt(FORT / "sw.txt")
    lw = np.loadtxt(FORT / "lw.txt")
    valid = [i for i in range(1, T.shape[0]) if np.all(np.isfinite(T[i]))]
    k = valid[-1]
    return P, T[k], sw[k - 1], lw[k - 1]


def build_column(P, T):
    """Build a Column from Fortran (P[Pa], T[K]) layer values (ascending altitude)."""
    o = np.argsort(P)[::-1]              # high P (surface) -> low P (top)
    P, T = P[o], T[o]
    # hydrostatic z from the surface up: dz = R T / (M g) dlnP
    z = np.zeros_like(P)
    for i in range(1, P.size):
        Tm = 0.5 * (T[i] + T[i - 1])
        g = gravity(z[i - 1])
        z[i] = z[i - 1] + (C.R_GAS * Tm / (0.028 * g)) * np.log(P[i - 1] / P[i])
    g = gravity(z)
    return Column(z, T, P, g)            # ascending altitude


def main():
    P, Tf, sw_f, lw_f = fortran_state()
    print(f"Using Fortran snapshot (P {P.min():.2f}-{P.max():.1e} Pa, "
          f"T {Tf.min():.0f}-{Tf.max():.0f} K)\n")

    col = build_column(P, Tf)
    micro = solve_bvp_profile(Atmosphere.titan_reference(), DEFAULT, n_nodes=200)
    op = OpticsParams(prescribed_haze=True)          # same haze the Fortran uses
    ck = CorrelatedKSW(sma_au=9.58)
    cklw = CorrelatedKLW()
    fx = eb.compute_fluxes(col, micro, op=op, ck=ck, ck_lw=cklw)

    # our heating in K/s on the same layers (Fortran heating is K/s)
    sw_ours = fx.sw_heating / TITAN_DAY
    lw_ours = fx.lw_heating / TITAN_DAY

    # --- 1. column energetics ---
    print("=== TOA / column energetics (same T, prescribed haze) ===")
    print(f"  ours : SW absorbed {fx.sw_net[-1]-fx.sw_net[0]:6.3f}  OLR {-fx.lw_net[-1]:6.3f} W/m^2")
    # Fortran integrated heating -> flux divergence (rough): sum(rho cp H dz)
    print(f"  Fortran sw/lw heating (K/s) ranges: SW {sw_f.max():.2e}  LW {lw_f.min():.2e}")

    # --- 2. heating-rate profile comparison (interp to common pressure) ---
    Pc = np.logspace(0, np.log10(1.3e5), 12)
    o = np.argsort(col.P_mid)
    sw_o = np.interp(Pc, col.P_mid[o], sw_ours[o])
    lw_o = np.interp(Pc, col.P_mid[o], lw_ours[o])
    of = np.argsort(P)
    sw_ff = np.interp(Pc, P[of], sw_f[of])
    lw_ff = np.interp(Pc, P[of], lw_f[of])
    print("\n=== heating rate [K/Titan-day] vs pressure (same T) ===")
    print(f"  {'P[Pa]':>9} {'SW ours':>9} {'SW fort':>9} {'LW ours':>9} {'LW fort':>9}")
    for i in range(len(Pc)):
        print(f"  {Pc[i]:9.1f} {sw_o[i]*TITAN_DAY:9.2f} {sw_ff[i]*TITAN_DAY:9.2f}"
              f" {lw_o[i]*TITAN_DAY:9.2f} {lw_ff[i]*TITAN_DAY:9.2f}")

    # --- 3. Rayleigh scattering (Fortran has it, we don't) ---
    ray = np.loadtxt(DATA / "Rayleigh.txt", skiprows=2)
    bc = ck.bands
    tauray_band = np.interp(bc, ray[:, 0], ray[:, 1])          # per-band coeff (1/Pa-ish)
    print("\n=== Rayleigh scattering (present in Fortran, ABSENT in ours) ===")
    print("  Rayleigh coeff rises steeply to short wavelength (~nu^4):")
    print(f"    band {bc[5]:.0f}: {tauray_band[5]:.2e}   band {bc[25]:.0f}: {tauray_band[25]:.2e}"
          f"   band {bc[-3]:.0f}: {tauray_band[-3]:.2e}")
    print("  -> a pure-scattering source in the blue/UV that our SW omits")

    # --- 4. CIA: our HITRAN band-avg vs Fortran exp-sum transmission fit ---
    cia = CIABands()
    tau_cia_ours = cia.optical_depth(col, Composition()).sum(axis=1)   # per band
    bc_lw = 0.5 * (cia.band_lo + cia.band_hi)
    print("\n=== CIA column optical depth: ours (HITRAN) vs band ~100/300 cm^-1 ===")
    for nu in (100, 300, 500):
        j = int(np.argmin(np.abs(bc_lw - nu)))
        print(f"  nu~{nu:3d} cm-1: ours column tau = {tau_cia_ours[j]:.2f}")
    print("  (Fortran uses exp-sum transmission fits trans_*.txt -- different data;")
    print("   compare magnitudes/shape to assess the LW continuum discrepancy)")


if __name__ == "__main__":
    main()
