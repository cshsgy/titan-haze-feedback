"""Structural checks for the Step 1 DISORT energy balance.

Requires the pydisort runtime; run with the RT venv:
    .rtenv/bin/python tests/test_rt.py
"""

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import numpy as np

from microphysics import Atmosphere, DEFAULT, solve_bvp_profile
from rt.column import Column
from rt.energy_balance import compute_fluxes, radiative_equilibrium


def _setup(nlyr=30):
    atm = Atmosphere.titan_reference()
    micro = solve_bvp_profile(atm, DEFAULT, n_nodes=150)
    col = Column.from_atmosphere(atm, nlyr=nlyr, z_top=450e3)
    return atm, micro, col


def test_disort_imports():
    import pydisort  # noqa: F401
    import torch     # noqa: F401


def test_fluxes_finite():
    _, micro, col = _setup()
    fx = compute_fluxes(col, micro)
    assert np.all(np.isfinite(fx.sw_net)) and np.all(np.isfinite(fx.lw_net))
    assert np.all(np.isfinite(fx.net_heating))


def test_shortwave_absorbed_fraction_sane():
    """SW net down is positive, decreases downward, and 0 < absorbed <= TOA."""
    _, micro, col = _setup()
    fx = compute_fluxes(col, micro)
    toa, surf = fx.sw_net[-1], fx.sw_net[0]
    assert toa > 0
    assert -1e-6 <= surf <= toa + 1e-6
    # monotone non-increasing downward (allow tiny numerical wiggle)
    assert np.all(np.diff(fx.sw_net) >= -1e-6)


def test_outgoing_longwave_positive():
    """OLR (net upward LW at TOA) must be positive for a warm atmosphere."""
    _, micro, col = _setup()
    fx = compute_fluxes(col, micro)
    olr = -fx.lw_net[-1]
    assert olr > 0, olr


def test_haze_raises_shortwave_albedo():
    """More haze => more shortwave scattered back => smaller TOA net down."""
    _, micro, col = _setup()
    fx_full = compute_fluxes(col, micro)
    # thin the haze 10x by scaling number density
    import copy
    thin = copy.copy(micro)
    thin.n = micro.n * 0.1
    fx_thin = compute_fluxes(col, thin)
    assert fx_full.sw_net[-1] < fx_thin.sw_net[-1], \
        (fx_full.sw_net[-1], fx_thin.sw_net[-1])


def test_equilibrium_reduces_residual():
    """Relaxation should cut the bulk net-heating residual substantially."""
    _, micro, col = _setup()
    eq, hist = radiative_equilibrium(col, micro, n_iter=200)
    assert hist[-50:].mean() < 0.25 * hist[0], (hist[0], hist[-50:].mean())
    # equilibrium profile stays physical
    assert np.all((eq.T > 50) & (eq.T < 300))


if __name__ == "__main__":
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    failed = 0
    for t in tests:
        try:
            t()
            print(f"PASS  {t.__name__}")
        except AssertionError as e:
            failed += 1
            print(f"FAIL  {t.__name__}: {e}")
    print(f"\n{len(tests) - failed}/{len(tests)} passed")
    sys.exit(1 if failed else 0)
