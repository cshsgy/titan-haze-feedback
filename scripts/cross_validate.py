#!/usr/bin/env python3
"""Cross-validate the microphysics scaling law / BVP against published Titan
haze constraints.

Checks:
  1. Haze extinction scale height above 80 km vs Tomasko et al. (2008): H ~ 65 km.
  2. Characteristic particle radius in the main haze vs de Trenquelleon et al.
     (2025): r ~ 0.5 um (and Tomasko monomer/aggregate sizes).
  3. Global mass balance: surface deposition flux == column production.
  4. Consistency of the eddy-diffusion BVP with the K->0 master ODE.

Run from the repo root:
    python3 scripts/cross_validate.py
Writes writing/figs/bvp_validation.png
"""

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from microphysics import Atmosphere, DEFAULT, solve_scaling_law, solve_bvp_profile


def extinction(res):
    """Geometric extinction coefficient ~ n * projected area ~ n * r_a^2.

    Wavelength dependence and the absolute Q_ext cancel in a scale-height fit,
    so this suffices for comparing the *shape* of the opacity profile.
    """
    return res.n * res.r_a**2


def scale_height(z, k, z_lo, z_hi):
    """Exponential scale height H [m] of k(z) fit over [z_lo, z_hi]."""
    m = (z >= z_lo) & (z <= z_hi) & (k > 0)
    slope = np.polyfit(z[m], np.log(k[m]), 1)[0]
    return -1.0 / slope


def main():
    atm = Atmosphere.titan_reference()
    master = solve_scaling_law(atm, DEFAULT)
    bvp = solve_bvp_profile(atm, DEFAULT, n_nodes=300)

    print(f"BVP converged: {bvp.success}  ({bvp.z.size} nodes)\n")

    # --- 1. extinction scale height (Tomasko 2008: ~65 km above 80 km) ---
    print("1. Extinction scale height above 80 km  (Tomasko 2008: ~65 km)")
    for name, res in (("master (K->0)", master), ("BVP (eddy)", bvp)):
        H = scale_height(res.z, extinction(res), 80e3, 300e3) / 1e3
        print(f"   {name:14s}: H = {H:5.1f} km")

    # --- 2. characteristic radius in the main haze (dT25: ~0.5 um) ---
    print("\n2. Mass radius r(z)  (dT25 main-haze aerosols ~0.5 um; r_c=0.46 um)")
    print(f"   {'z[km]':>6}{'master r[um]':>14}{'BVP r[um]':>12}")
    for zk in (300, 200, 150, 100, 40):
        rm = np.interp(zk * 1e3, master.z[::-1], master.r[::-1]) * 1e6
        rb = np.interp(zk * 1e3, bvp.z, bvp.r) * 1e6
        print(f"   {zk:6d}{rm:14.3f}{rb:12.3f}")

    # --- 3. global mass balance ---
    P = DEFAULT.P_flux
    dep = bvp.omega[np.argmin(bvp.z)] * bvp.M[np.argmin(bvp.z)]  # surface settling flux
    print("\n3. Global mass balance (surface deposition / column production)")
    print(f"   surface settling flux / P = {dep / P:.4f}  (should be ~1)")
    print(f"   monomer flux Phi_M / P in "
          f"[{bvp.flux_M.min()/P:.4f}, {bvp.flux_M.max()/P:.4f}]")

    # --- 4. BVP vs master agreement in the (settling-dominated) lower haze ---
    zc = np.linspace(0, 250e3, 50)
    rb = np.interp(zc, bvp.z, bvp.r)
    rm = np.interp(zc, master.z[::-1], master.r[::-1])
    rel = np.abs(rb - rm) / rm
    print("\n4. BVP vs master ODE, mass radius, 0-250 km")
    print(f"   mean |rel. diff| = {rel.mean()*100:.1f}%  (max {rel.max()*100:.1f}%)")

    # --- figure ---
    fig, ax = plt.subplots(1, 3, figsize=(13, 5), sharey=True)
    zkm_m, zkm_b = master.z / 1e3, bvp.z / 1e3

    ax[0].plot(master.r * 1e6, zkm_m, label="master (K->0)")
    ax[0].plot(bvp.r * 1e6, zkm_b, "--", label="BVP (eddy diff.)")
    ax[0].axvspan(0.4, 0.5, color="green", alpha=0.15, label="dT25 ~0.5 um")
    ax[0].set_xscale("log"); ax[0].set_xlabel(r"mass radius [$\mu$m]")
    ax[0].set_ylabel("altitude [km]"); ax[0].set_title("particle size"); ax[0].legend()

    ax[1].plot(master.n * 1e-6, zkm_m, label="master")
    ax[1].plot(bvp.n * 1e-6, zkm_b, "--", label="BVP")
    ax[1].set_xscale("log"); ax[1].set_xlabel(r"number density [cm$^{-3}$]")
    ax[1].set_title("number density"); ax[1].legend()

    # extinction with the Tomasko 65 km reference slope anchored at 80 km
    kb = extinction(bvp)
    ax[2].plot(kb / kb.max(), zkm_b, "--", color="C1", label="BVP n·$r_a^2$")
    km = extinction(master)
    ax[2].plot(km / km.max(), zkm_m, color="C0", label="master n·$r_a^2$")
    z_ref = np.linspace(80e3, 300e3, 20)
    k80 = np.interp(80e3, bvp.z, kb)
    ax[2].plot((k80 * np.exp(-(z_ref - 80e3) / 65e3)) / kb.max(), z_ref / 1e3,
               ":", color="k", label="Tomasko H=65 km")
    ax[2].set_xscale("log"); ax[2].set_xlabel("normalized extinction")
    ax[2].set_title("opacity shape"); ax[2].legend()

    for a in ax:
        a.grid(alpha=0.3); a.set_ylim(0, 420)
    fig.suptitle("Cross-validation of the Titan haze microphysics against published constraints")
    fig.tight_layout()
    out = ROOT / "writing" / "figs" / "bvp_validation.png"
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=130)
    print(f"\nwrote {out.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
