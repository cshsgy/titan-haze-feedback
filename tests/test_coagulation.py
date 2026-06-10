#!/usr/bin/env python3
"""Unit tests for the bimodal coagulation moment tendencies (piece 1b).

    python3 tests/test_coagulation.py
"""

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import numpy as np

from microphysics import Atmosphere, DEFAULT
from microphysics import coagulation as cg
from microphysics import transport as tr

_P = _F = 0


def chk(name, cond, info=""):
    global _P, _F
    if cond:
        _P += 1; print(f"PASS  {name}")
    else:
        _F += 1; print(f"FAIL  {name}  {info}")


def _state(N0S=0.3, M0S=1e8, N0F=5e3, M0F=1e6, d=None):
    d = d or DEFAULT.d_mono
    M3S = M0S * (d * N0S ** (1 / 3)) ** 3
    M3F = M0F * (d * N0F ** (1 / 3)) ** 3
    return M0S, M3S, M0F, M3F


def test_volume_conservation():
    """Coagulation must conserve total monomer volume: dM3S + dM3F = 0."""
    atm = Atmosphere.titan_reference()
    for z in (300e3, 200e3, 100e3, 30e3):
        s = _state()
        d = cg.coag_tendencies(*s, z, atm, DEFAULT, 1.5, 2.0, n=24)
        dM3_tot = d[1] + d[3]
        scale = abs(d[1]) + abs(d[3]) + 1e-300
        chk(f"volume conserved at z={z/1e3:.0f}km",
            abs(dM3_tot) / scale < 1e-10, f"rel {abs(dM3_tot)/scale:.1e}")


def test_number_decreases():
    """Coagulation reduces total particle number: dM0S + dM0F <= 0."""
    atm = Atmosphere.titan_reference()
    s = _state()
    d = cg.coag_tendencies(*s, 200e3, atm, DEFAULT, 1.5, 2.0, n=24)
    chk("total number decreases", d[0] + d[2] <= 0, f"dN={d[0]+d[2]:.2e}")


def test_SS_transfer_feeds_F():
    """With only the spherical mode, S+S coagulation should LOSE S number and
    (since aggregated products exceed a monomer) FEED the fractal mode."""
    atm = Atmosphere.titan_reference()
    # spherical mode of fairly large spheres so products exceed the monomer size
    M0S, M3S, _, _ = _state(N0S=3.0, M0S=1e8)
    d = cg.coag_tendencies(M0S, M3S, 0.0, 0.0, 200e3, atm, DEFAULT, 1.5, 2.0, n=24)
    chk("S number decreases", d[0] < 0, f"dM0S={d[0]:.2e}")
    chk("F number fed by S+S", d[2] > 0, f"dM0F={d[2]:.2e}")
    chk("volume moved S->F conserved", abs(d[1] + d[3]) / (abs(d[1]) + 1e-300) < 1e-10)


def test_monodisperse_limit():
    """A single narrow (sigma->1) fractal mode must recover the monodisperse
    self-coagulation number loss -0.5 beta(N0) M0^2 (the existing coag_kernel)."""
    atm = Atmosphere.titan_reference()
    p = DEFAULT
    z = 150e3
    N0, M0 = 5e3, 1e6
    M3 = M0 * (p.d_mono * N0 ** (1 / 3)) ** 3
    d = cg.coag_tendencies(0.0, 0.0, M0, M3, z, atm, p, 1.5, 1.0 + 1e-6, n=40)
    beta = tr.coag_kernel(N0, z, atm, p)
    ref = -0.5 * beta * M0 ** 2
    chk("monodisperse FF number loss == -0.5 beta n^2",
        np.isclose(d[2], ref, rtol=2e-3), f"{d[2]:.4e} vs {ref:.4e}")


def test_charge_factor():
    """Coulomb inhibition: W=1 at n_e=0, decreases with size, ~0.3 at 0.1 um."""
    W0 = tr.charge_factor(1e-7, 1e-7, 160.0, 0.0)
    W1 = tr.charge_factor(1e-7, 1e-7, 160.0, 15.0)
    W5 = tr.charge_factor(5e-7, 5e-7, 160.0, 15.0)
    chk("charge factor unity at n_e=0", np.isclose(W0, 1.0))
    chk("charge factor ~0.3 at r_a=0.1 um", 0.2 < W1 < 0.45, f"{W1:.3f}")
    chk("charge factor decreases with size", W5 < 0.01 * W1, f"{W5:.2e} vs {W1:.3f}")


def test_charge_inhibits_both_kernels():
    """use_charge slows mono kernel and bimodal tendencies; volume conserved."""
    from dataclasses import replace
    atm = Atmosphere.titan_reference()
    pc = replace(DEFAULT, use_charge=True)
    z = 150e3
    b0 = tr.coag_kernel(5e3, z, atm, DEFAULT)
    b1 = tr.coag_kernel(5e3, z, atm, pc)
    chk("mono kernel reduced by charge", b1 < 0.1 * b0, f"{b1:.3e} vs {b0:.3e}")
    s = _state()
    d0 = cg.coag_tendencies(*s, z, atm, DEFAULT, 1.5, 2.0)
    d1 = cg.coag_tendencies(*s, z, atm, pc, 1.5, 2.0)
    chk("bimodal number loss reduced by charge", abs(d1[2]) < abs(d0[2]),
        f"{d1[2]:.3e} vs {d0[2]:.3e}")
    chk("bimodal volume conserved with charge",
        abs(d1[1] + d1[3]) <= 1e-10 * max(abs(d1[1]), abs(d1[3]), 1e-300),
        f"{d1[1]:.3e}+{d1[3]:.3e}")


if __name__ == "__main__":
    test_volume_conservation()
    test_number_decreases()
    test_SS_transfer_feeds_F()
    test_monodisperse_limit()
    test_charge_factor()
    test_charge_inhibits_both_kernels()
    print(f"\n{_P}/{_P + _F} passed")
    sys.exit(1 if _F else 0)
