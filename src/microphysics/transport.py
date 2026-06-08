"""Size-dependent transport coefficients: fractal geometry, settling velocity
omega(N, z), and the Brownian coagulation kernel beta(N, N; z).

All expressions follow Burgalat & Rannou (2017) [BR17] and the derivation in
``docs/scaling_law.md``.  N is a *continuous monomer-volume count* (N < 1 is a
sub-monomer sphere), and the fractal dimension is two-phase: D = 3 for N < 1
(compact spheres), D = D_f for N >= 1 (fractal aggregates).
"""

from __future__ import annotations

import numpy as np

from . import constants as C


def fractal_dimension(N, D_f: float):
    """Two-phase fractal dimension D(N): 3 for N < 1, D_f for N >= 1."""
    return np.where(np.asarray(N) < 1.0, 3.0, D_f)


def mass_radius(N, d: float):
    """Volume-equivalent (mass) radius r = d N^{1/3} [m]."""
    return d * np.asarray(N, dtype=float) ** (1.0 / 3.0)


def mobility_radius(N, d: float, D_f: float):
    """Mobility (drag) radius r_a = d N^{1/D(N)} [m], two-phase D."""
    D = fractal_dimension(N, D_f)
    return d * np.asarray(N, dtype=float) ** (1.0 / D)


def settling_velocity(N, z, atm, p):
    """Stokes settling velocity with first-order Cunningham slip, fractal D.

    omega = (2 rho_p g d^2 / 9 eta) [ N^{(D-1)/D} + (A lambda / d) N^{(D-2)/D} ]
    (BR17 Eq. 29).  Returns the downward speed [m/s] (>0).
    """
    N = np.asarray(N, dtype=float)
    D = fractal_dimension(N, p.D_f)
    g = atm.gravity(z)
    eta = atm.viscosity(z)
    lam = atm.mfp(z)

    pref = 2.0 * p.rho_p * g * p.d_mono**2 / (9.0 * eta)
    stokes = N ** ((D - 1.0) / D)
    slip = (p.A_slip * lam / p.d_mono) * N ** ((D - 2.0) / D)
    return pref * (stokes + slip)


def coag_kernel(N, z, atm, p):
    """Equal-size Brownian coagulation kernel beta(N, N; z) [m^3/s].

    Fuchs harmonic bridge of the continuum (BR17 Eq. 10) and free-molecular
    (BR17 Eq. 13) forms, evaluated for two identical aggregates of size N:

      beta_CO = (8 k_B T / 3 eta) (1 + A Kn),     Kn = lambda / r_a
      beta_FM = 4 sqrt(2) (6 k_B T / rho_p)^{1/2} r_a^2 r^{-3/2}
      beta    = beta_CO beta_FM / (beta_CO + beta_FM)

    The optional charge-inhibition factor (BR17 Eq. 27) is deferred.
    """
    N = np.asarray(N, dtype=float)
    T = atm.temperature(z)
    eta = atm.viscosity(z)
    lam = atm.mfp(z)

    r = mass_radius(N, p.d_mono)
    r_a = mobility_radius(N, p.d_mono, p.D_f)
    Kn = lam / r_a

    beta_co = (8.0 * C.K_B * T / (3.0 * eta)) * (1.0 + p.A_slip * Kn)
    beta_fm = 4.0 * np.sqrt(2.0) * np.sqrt(6.0 * C.K_B * T / p.rho_p) * r_a**2 * r ** (-1.5)
    return beta_co * beta_fm / (beta_co + beta_fm)
