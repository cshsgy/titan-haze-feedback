"""Sanity + cross-validation checks for the eddy-diffusion BVP.

Run from the repo root:
    python3 tests/test_bvp.py
"""

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import numpy as np

from microphysics import Atmosphere, DEFAULT, solve_scaling_law, solve_bvp_profile


def _bvp():
    return solve_bvp_profile(Atmosphere.titan_reference(), DEFAULT, n_nodes=300)


def test_bvp_converges():
    assert _bvp().success


def test_monomer_flux_constant_equals_P():
    """No interior source => downward monomer flux is constant and equals P."""
    res = _bvp()
    P = DEFAULT.P_flux
    assert np.allclose(res.flux_M / P, 1.0, rtol=1e-3), \
        (res.flux_M.min() / P, res.flux_M.max() / P)


def test_surface_deposition_balances_production():
    """Global mass balance: surface settling flux == column production."""
    res = _bvp()
    i = int(np.argmin(res.z))
    dep = res.omega[i] * res.M[i]
    assert np.isclose(dep / DEFAULT.P_flux, 1.0, rtol=1e-3)


def test_bvp_agrees_with_master_in_lower_haze():
    """Eddy diffusion is a modest correction below 250 km (settling-dominated)."""
    atm = Atmosphere.titan_reference()
    bvp = solve_bvp_profile(atm, DEFAULT, n_nodes=300)
    master = solve_scaling_law(atm, DEFAULT)
    z = np.linspace(0, 250e3, 40)
    rb = np.interp(z, bvp.z, bvp.r)
    rm = np.interp(z, master.z[::-1], master.r[::-1])
    rel = np.abs(rb - rm) / rm
    assert rel.mean() < 0.25, f"mean rel diff {rel.mean():.2f} too large"


def test_extinction_scale_height_matches_tomasko():
    """Extinction scale height above 80 km should be near Tomasko (2008) 65 km."""
    res = _bvp()
    k = res.n * res.r_a**2
    m = (res.z >= 80e3) & (res.z <= 300e3) & (k > 0)
    H = -1.0 / np.polyfit(res.z[m], np.log(k[m]), 1)[0] / 1e3
    assert 45.0 < H < 90.0, f"extinction scale height {H:.1f} km off Tomasko ~65 km"


def test_characteristic_radius_submicron():
    """Main-haze characteristic radius near dT25's ~0.5 um."""
    res = _bvp()
    r300 = np.interp(300e3, res.z, res.r) * 1e6
    assert 0.1 < r300 < 1.5, f"r(300 km) = {r300:.2f} um out of range"


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
