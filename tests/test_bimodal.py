#!/usr/bin/env python3
"""Tests for the K->0 bimodal master ODE (polydisperse piece 2).

    python3 tests/test_bimodal.py
"""

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import numpy as np

from microphysics import Atmosphere, DEFAULT
from microphysics import moments as mm
from microphysics.scaling_law_bimodal import solve_bimodal_kzero

_P = _F = 0


def chk(name, cond, info=""):
    global _P, _F
    if cond:
        _P += 1; print(f"PASS  {name}")
    else:
        _F += 1; print(f"FAIL  {name}  {info}")


def test_runs_finite():
    r = solve_bimodal_kzero(n_out=100, n_quad=14)
    chk("integration reaches surface", r.z.min() < 1e3)
    chk("all profiles finite", np.all(np.isfinite(r.rho_h))
        and np.all(np.isfinite(r.M0F)) and np.all(np.isfinite(r.M3F)))
    chk("positive moments", np.all(r.M0F >= 0) and np.all(r.M3F >= 0))
    return r


def test_volume_flux_conserved(r):
    """K->0: coagulation conserves volume, so the total downward volume flux
    Phi3 = <w>3^S M3^S + <w>3^F M3^F must equal the production flux P3 at all z."""
    p = DEFAULT
    P3 = (3.0 / (4.0 * np.pi)) * p.Q_mass / p.rho_p
    flux = np.empty_like(r.z)
    for i, z in enumerate(r.z):
        _, w3s = mm.settling_velocities(r.M0S[i], r.M3S[i], r.sigma_s, 3.0, z, r.atm, p)
        _, w3f = mm.settling_velocities(r.M0F[i], r.M3F[i], r.sigma_f, p.D_f, z, r.atm, p)
        flux[i] = w3s * r.M3S[i] + w3f * r.M3F[i]
    # below the production region (exclude the top few km where it's injected)
    deep = r.z < 0.95 * p.z0
    rel = np.abs(flux[deep] - P3) / P3
    chk("total volume flux == production (mass conservation)",
        rel.max() < 0.05, f"max rel dev {rel.max():.3f}")


def test_bimodal_handoff(r):
    """Production seeds the spherical mode at the top; the fractal mode takes
    over below as aggregates form (the S->F transfer)."""
    zk = r.z / 1e3
    top = int(np.argmax(zk))                       # ~z0
    fS = r.M3S / (r.M3S + r.M3F + 1e-300)           # S mass fraction
    chk("spherical mode dominates at production", fS[top] > 0.9, f"fS={fS[top]:.2f}")
    surf = int(np.argmin(zk))
    chk("fractal mode dominates at the surface", fS[surf] < 0.05, f"fS={fS[surf]:.3f}")
    chk("surface aggregate radius ~ micron", 0.3e-6 < r.r0F[surf] < 10e-6,
        f"r0F={r.r0F[surf]*1e6:.2f}um")


def test_mass_bottom_heavy(r):
    """Haze mass density increases downward (rains down), as for the monodisperse
    model and Titan observations."""
    zk = r.z / 1e3
    rho_top = r.rho_h[zk > 350].mean()
    rho_low = r.rho_h[(zk > 80) & (zk < 150)].mean()
    chk("rho_h bottom-heavy", rho_low > 3 * rho_top, f"low/top={rho_low/rho_top:.1f}")


if __name__ == "__main__":
    r = test_runs_finite()
    test_volume_flux_conserved(r)
    test_bimodal_handoff(r)
    test_mass_bottom_heavy(r)
    print(f"\n{_P}/{_P + _F} passed")
    sys.exit(1 if _F else 0)
