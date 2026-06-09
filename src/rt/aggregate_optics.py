"""Mean-field (Rayleigh-Debye-Gans) optics for fractal haze aggregates.

The gray ``Q_ext * pi * r_a^2`` cross-section (mobility radius) overestimates the
extinction of large fractal aggregates by ~8x.  RDG is the standard treatment for
Titan haze: each monomer absorbs independently, so the aggregate absorption is

    C_abs(N, lambda) = N * C_abs,monomer(lambda),
    C_abs,monomer    = (8 pi^2 d^3 / lambda) * Im{(m^2-1)/(m^2+2)},

i.e. it scales with the number of monomers (mass), NOT with r_a^2.  The monomer
is in the Rayleigh regime; ``m(lambda)`` is the tholin refractive index (Khare et
al. 1984).  We obtain the extinction from the absorption and an observational
single-scattering albedo, ``C_ext = C_abs / (1 - omega0)`` -- avoiding the RDG-FA
scattering structure factor while fixing the magnitude and wavelength dependence.
"""

from __future__ import annotations

import numpy as np

# Tholin complex refractive index m = n + i k (Khare et al. 1984, representative).
# wavelength [um], n, k.
_LAM = np.array([0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.70,
                 0.80, 0.90, 1.00, 1.25, 1.50, 2.00, 2.50, 3.00, 3.50, 4.00,
                 5.00, 6.00, 8.00, 10.0, 15.0, 20.0, 30.0, 50.0])
_N = np.array([1.65, 1.78, 1.69, 1.66, 1.65, 1.64, 1.63, 1.62, 1.61, 1.60,
               1.59, 1.58, 1.58, 1.57, 1.57, 1.57, 1.56, 1.56, 1.62, 1.61,
               1.62, 1.64, 1.68, 1.72, 1.78, 1.82, 1.88, 1.95])
_K = np.array([3.85e-1, 2.05e-1, 9.28e-2, 6.23e-2, 3.86e-2, 2.34e-2, 1.36e-2,
               7.38e-3, 3.87e-3, 1.48e-3, 1.02e-3, 8.00e-4, 6.92e-4, 9.0e-4,
               1.48e-3, 5.18e-3, 1.66e-2, 6.82e-2, 4.66e-2, 2.07e-2, 4.31e-2,
               6.48e-2, 8.5e-2, 1.10e-1, 1.6e-1, 2.0e-1, 2.6e-1, 3.4e-1])


def tholin_index(lam_um):
    """Tholin m(lambda) = n + i k, interpolated (log-wavelength)."""
    lam = np.clip(np.asarray(lam_um, float), _LAM[0], _LAM[-1])
    ll = np.log(lam)
    n = np.interp(ll, np.log(_LAM), _N)
    k = np.interp(ll, np.log(_LAM), _K)
    return n + 1j * k


def monomer_cabs(lam_um, d):
    """Rayleigh monomer absorption cross-section [m^2] (d = monomer radius [m])."""
    lam = np.asarray(lam_um, float) * 1e-6                 # m
    m = tholin_index(lam_um)
    F = (m**2 - 1.0) / (m**2 + 2.0)
    return (8.0 * np.pi**2 * d**3 / lam) * np.abs(F.imag)


def band_lambda_um(band_centers_cm):
    """Wavelength [um] of band centres given in (true) cm^-1."""
    return 1.0e4 / np.asarray(band_centers_cm, float)


def aggregate_haze_layer_tau(column, micro, band_centers_cm, omega0_band, d_mono,
                             pure_absorber=False):
    """RDG haze optical depth tau[nband, nlyr] from the microphysics profile.

    tau(band, z) = n(z) * Nbar(z) * C_abs,monomer(lambda_band) / (1 - omega0) * dz,
    i.e. absorption is additive over monomers (~mass) and extinction follows from
    the (observational) single-scattering albedo.  Set ``pure_absorber`` for the
    thermal IR (omega0 ~ 0), where extinction == absorption.
    """
    lam = band_lambda_um(band_centers_cm)                 # (nband,)
    cabs = monomer_cabs(lam, d_mono)                      # (nband,) per monomer
    o = np.argsort(micro.z)
    zc = column.z_mid
    n = np.interp(zc, micro.z[o], micro.n[o])
    Nbar = np.interp(zc, micro.z[o], micro.Nbar[o])
    monomer_col = (n * Nbar)[None, :]                     # (1, nlyr)
    if pure_absorber:
        denom = 1.0
    else:
        denom = np.clip(1.0 - np.asarray(omega0_band, float), 1e-3, 1.0)[:, None]
    ext = monomer_col * cabs[:, None] / denom             # (nband, nlyr) [1/m]
    return np.maximum(ext * column.dz[None, :], 0.0)
