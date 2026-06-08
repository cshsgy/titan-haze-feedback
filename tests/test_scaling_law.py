"""Sanity checks for the microphysics scaling law.

Run from the repo root:
    python3 -m pytest tests/        (if pytest available)
    python3 tests/test_scaling_law.py   (standalone)
"""

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import numpy as np

from microphysics import Atmosphere, AerosolParams, DEFAULT, solve_scaling_law
from microphysics import transport as tr


def test_geometry_switch_continuous_at_N1():
    """r and r_a coincide at N=1 (both equal the monomer radius d)."""
    p = DEFAULT
    assert np.isclose(tr.mass_radius(1.0, p.d_mono), p.d_mono)
    assert np.isclose(tr.mobility_radius(1.0, p.d_mono, p.D_f), p.d_mono)
    # spheres (N<1) have r_a == r; fractals (N>1) are more open (r_a > r for D<3)
    assert np.isclose(tr.mobility_radius(0.1, p.d_mono, p.D_f),
                      tr.mass_radius(0.1, p.d_mono))
    assert tr.mobility_radius(1e3, p.d_mono, p.D_f) > tr.mass_radius(1e3, p.d_mono)


def test_fractal_dimension_two_phase():
    assert tr.fractal_dimension(0.5, 2.0) == 3.0
    assert tr.fractal_dimension(5.0, 2.0) == 2.0


def test_mass_flux_conserved():
    """rho_h * omega must equal the imposed mass production flux everywhere."""
    res = solve_scaling_law(Atmosphere.titan_reference(), DEFAULT)
    ratio = res.rho_h * res.omega / DEFAULT.Q_mass
    assert np.allclose(ratio, 1.0, rtol=1e-6), (ratio.min(), ratio.max())


def test_size_grows_downward():
    """Aggregates grow as they fall: Nbar increases with decreasing altitude."""
    res = solve_scaling_law(Atmosphere.titan_reference(), DEFAULT)
    order = np.argsort(res.z)        # ascending altitude
    Nbar_asc = res.Nbar[order]
    assert np.all(np.diff(Nbar_asc) <= 1e-9), "Nbar should be monotone in altitude"
    assert res.Nbar[np.argmin(res.z)] > res.Nbar[np.argmax(res.z)]


def test_reaches_micron_sizes():
    """Surface aggregates should reach ~0.1-10 um (cf. dT25: aerosols ~0.5 um)."""
    res = solve_scaling_law(Atmosphere.titan_reference(), DEFAULT)
    r_surf = res.r[np.argmin(res.z)]
    assert 1e-8 < r_surf < 1e-5, f"surface mass radius {r_surf:.2e} m out of range"


def test_free_molecular_exponent():
    """Near the source (free-molecular, D=2) the master ODE slope d ln Nbar/dz
    should match the analytic exponent: dNbar/dz ~ -Nbar^{1/2}."""
    from microphysics.scaling_law import _dlnN_dz
    atm = Atmosphere.titan_reference()
    p = DEFAULT
    z = 380e3  # just below the source, free-molecular regime
    # local slope of dNbar/dz vs Nbar in log-log should be ~ +1/2 (s=1/2)
    Ns = np.array([2.0, 4.0, 8.0])
    dN = []
    for N in Ns:
        omega = tr.settling_velocity(N, z, atm, p)
        beta = tr.coag_kernel(N, z, atm, p)
        dN.append(beta * p.P_flux / (2.0 * omega**2))  # |dNbar/dz|
    slope = np.polyfit(np.log(Ns), np.log(dN), 1)[0]
    assert abs(slope - 0.5) < 0.15, f"FM exponent {slope:.3f} != ~0.5"


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
