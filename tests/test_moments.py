#!/usr/bin/env python3
"""Unit tests for the log-normal moment machinery (polydisperse Step 2, piece 1a).

    python3 tests/test_moments.py
"""

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import numpy as np

from microphysics import Atmosphere, DEFAULT
from microphysics import moments as mm

_P = _F = 0


def chk(name, cond, info=""):
    global _P, _F
    if cond:
        _P += 1; print(f"PASS  {name}")
    else:
        _F += 1; print(f"FAIL  {name}  {info}")


def test_alpha():
    chk("alpha(0)=1", np.isclose(mm.alpha(0.0, 2.0), 1.0))
    chk("alpha(3)=1", np.isclose(mm.alpha(3.0, 2.0), 1.0))
    chk("alpha(k,sigma=1)=1", np.allclose([mm.alpha(k, 1.0) for k in (-1, 1, 2, 5)], 1.0))
    # broader sigma -> larger high-order moments (alpha grows for k>3 and k<0)
    chk("alpha(6) grows with sigma", mm.alpha(6.0, 2.0) > mm.alpha(6.0, 1.3) > 1.0)


def test_geometry():
    M0, M3 = 1e6, 1e6 * (0.3e-6) ** 3            # r0 = 0.3 um
    chk("mean_radius", np.isclose(mm.mean_radius(M0, M3), 0.3e-6))
    chk("E=1 for spheres", np.isclose(mm.E_factor(3.0, DEFAULT.d_mono), 1.0))
    chk("N_mean", np.isclose(mm.N_mean(M0, M3, DEFAULT.d_mono),
                             (0.3e-6 / DEFAULT.d_mono) ** 3))


def test_settling_monodisperse_limit():
    """sigma -> 1 must reproduce the single-particle settling at the mean size."""
    atm = Atmosphere.titan_reference()
    p = DEFAULT
    for Df, Nbar, z in [(2.0, 5e3, 250e3), (2.0, 1e2, 150e3), (3.0, 50.0, 300e3)]:
        r0 = p.d_mono * Nbar ** (1.0 / 3.0)
        M0 = 1e7
        M3 = M0 * r0 ** 3
        w0, w3 = mm.settling_velocities(M0, M3, 1.0, Df, z, atm, p)
        # reference: the existing N-based law with the SAME fixed Df
        pref = 2.0 * p.rho_p * atm.gravity(z) * p.d_mono ** 2 / (9.0 * atm.viscosity(z))
        lam = atm.mfp(z)
        w_ref = pref * (Nbar ** ((Df - 1) / Df)
                        + (p.A_slip * lam / p.d_mono) * Nbar ** ((Df - 2) / Df))
        chk(f"sigma->1 <w>_0 == single-particle (Df={Df},N={Nbar:.0e})",
            np.isclose(w0, w_ref, rtol=1e-9), f"{w0:.3e} vs {w_ref:.3e}")
        chk(f"sigma->1 <w>_0 == <w>_3 (Df={Df},N={Nbar:.0e})", np.isclose(w0, w3, rtol=1e-9))


def test_gravitational_sorting():
    """For sigma>1 the volume settles faster than the number: <w>_3 > <w>_0."""
    atm = Atmosphere.titan_reference()
    p = DEFAULT
    r0 = p.d_mono * (5e3) ** (1.0 / 3.0)
    M0 = 1e7; M3 = M0 * r0 ** 3
    for sigma in (1.5, 2.0, 2.5):
        w0, w3 = mm.settling_velocities(M0, M3, sigma, 2.0, 250e3, atm, p)
        chk(f"<w>_3 > <w>_0 at sigma={sigma}", w3 > w0 > 0, f"w0={w0:.3e} w3={w3:.3e}")
    # the ratio grows with sigma (more sorting)
    r = []
    for sigma in (1.3, 2.0, 2.6):
        w0, w3 = mm.settling_velocities(M0, M3, sigma, 2.0, 250e3, atm, p)
        r.append(w3 / w0)
    chk("sorting ratio increases with sigma", r[0] < r[1] < r[2], str([round(x, 2) for x in r]))


if __name__ == "__main__":
    test_alpha()
    test_geometry()
    test_settling_monodisperse_limit()
    test_gravitational_sorting()
    print(f"\n{_P}/{_P + _F} passed")
    sys.exit(1 if _F else 0)
