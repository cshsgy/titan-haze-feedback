"""Log-normal moment machinery for the bimodal 2-moment microphysics (Step 2+).

Polydisperse upgrade of the monodisperse closure: at each altitude every mode is
a log-normal in bulk radius ``r`` carried by two moments -- number ``M0`` and
volume ``M3`` (``M_n = integral r^n n(r) dr``).  See ``docs/polydisperse_scheme.md``
for the full Burgalat & Rannou (2017) scheme.  This module provides the pieces
that follow in closed form: the inter-moment function ``alpha(k)``, the mode
geometry, and the moment-averaged settling velocities (the size-dependent
gravitational sorting the monodisperse model could not represent).  The
coagulation closures (gamma coefficients, b_k bridging, inter-mode transfer) come
in a companion module.

Conventions (B&R 2017): bulk radius ``r = d_mono N^(1/3)``; apparent/mobility
radius ``r_a = E r^(3/Df)`` with ``E = d_mono^((Df-3)/Df)`` (so spheres Df=3 give
E=1, r_a=r); mean radius ``r0 = (M3/M0)^(1/3)`` and ``M_k = M0 r0^k alpha(k)``.
"""

from __future__ import annotations

import numpy as np


def alpha(k, sigma):
    """Log-normal inter-moment factor ``M_k / M0 = r0^k alpha(k)``.

    For a log-normal of geometric standard deviation ``sigma`` about the volume
    mean radius ``r0 = (M3/M0)^(1/3)``::

        alpha(k) = exp[ (k^2 - 3k) ln(sigma)^2 / 2 ],

    so ``alpha(0) = alpha(3) = 1`` and ``sigma -> 1`` (monodisperse) gives
    ``alpha = 1`` for all k.
    """
    s2 = np.log(np.asarray(sigma, float)) ** 2
    k = np.asarray(k, float)
    return np.exp(0.5 * (k * k - 3.0 * k) * s2)


def E_factor(Df, d_mono):
    """Apparent-to-bulk conversion factor ``E = d_mono^((Df-3)/Df)`` (B&R Eq. 8)."""
    return d_mono ** ((Df - 3.0) / Df)


def mean_radius(M0, M3):
    """Volume-mean bulk radius ``r0 = (M3 / M0)^(1/3)`` [m]."""
    return np.cbrt(np.asarray(M3, float) / np.maximum(np.asarray(M0, float), 1e-300))


def N_mean(M0, M3, d_mono):
    """Mean monomer count of the mode, ``N0 = (r0 / d_mono)^3``."""
    return (mean_radius(M0, M3) / d_mono) ** 3


def settling_velocities(M0, M3, sigma, Df, z, atm, p):
    """Moment-averaged settling speeds ``(<w>_0, <w>_3)`` [m/s, downward >0].

    From the single-particle Stokes+slip law (B&R Eq. 29) integrated over the
    log-normal (Eq. 30) and reduced via ``alpha``::

        <w>_k = (2 rho_p g)/(9 eta E)
                [ r0^(3-3/Df) alpha(k+3-3/Df-3+... ) + (A lam/E) r0^(3-6/Df) alpha(...) ]

    with the moment orders ``3-3/Df`` (Stokes) and ``3-6/Df`` (slip) in r0, and
    the alpha weights ``(3,6)-3/Df`` for k=(0,3) Stokes and ``(3,6)-6/Df`` for the
    slip term.  At ``sigma -> 1`` both collapse to the single-particle speed at
    r0, matching ``transport.settling_velocity``.
    """
    Df = float(Df)
    E = E_factor(Df, p.d_mono)
    g = atm.gravity(z)
    eta = atm.viscosity(z)
    lam = atm.mfp(z)
    r0 = mean_radius(M0, M3)

    pref = 2.0 * p.rho_p * g / (9.0 * eta * E)
    rs = r0 ** (3.0 - 3.0 / Df)            # Stokes radius power
    rp = r0 ** (3.0 - 6.0 / Df)            # slip radius power
    slip = (p.A_slip * lam / E)

    def wk(k):
        # alpha weights: Stokes term order (k+3)-3/Df reduced by /M_k -> 3+k-3/Df ... see doc
        a_stokes = alpha(k + 3.0 - 3.0 / Df, sigma) / alpha(k, sigma)
        a_slip = alpha(k + 3.0 - 6.0 / Df, sigma) / alpha(k, sigma)
        return pref * (rs * a_stokes + slip * rp * a_slip)

    return wk(0.0), wk(3.0)
