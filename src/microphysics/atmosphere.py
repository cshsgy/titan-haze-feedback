"""Background atmosphere: T(z), P(z), g(z), and the gas properties the
microphysics needs -- dynamic viscosity eta(T), mean free path lambda(T,P),
and the eddy diffusivity K(z).

For standalone microphysics runs a crude HASI-like Titan T(z) is provided and
P(z) is obtained by hydrostatic integration.  In the coupled model (Step 3) the
temperature profile is replaced by the DISORT-derived T(z): construct an
``Atmosphere`` from arbitrary (z, T) arrays via :meth:`from_profile`.
"""

from __future__ import annotations

import numpy as np

from . import constants as C

# Crude Titan reference temperature profile (HASI/Yelle-like anchors).
# altitude [m]  -> temperature [K].  Surface 94 K, tropopause ~71 K near 44 km,
# stratosphere warming to ~170 K, upper-atmosphere ~170 K plateau.
_Z_REF = np.array([0.0, 44e3, 100e3, 150e3, 200e3, 300e3, 415e3, 540e3])
_T_REF = np.array([94.0, 71.0, 110.0, 150.0, 170.0, 175.0, 170.0, 160.0])


def gravity(z: np.ndarray | float) -> np.ndarray | float:
    """Altitude-dependent gravity g(z) [m/s^2]."""
    return C.G_SURF * (C.R_TITAN / (C.R_TITAN + np.asarray(z))) ** 2


def viscosity_N2(T: np.ndarray | float) -> np.ndarray | float:
    """Dynamic viscosity of N2 [Pa s] via Sutherland's law.

    eta = eta0 (T/T0)^{3/2} (T0 + S)/(T + S), with N2 constants.  Swappable for a
    kinetic-theory expression if higher accuracy at <100 K is needed.
    """
    eta0, T0, S = 1.781e-5, 300.55, 111.0
    T = np.asarray(T, dtype=float)
    return eta0 * (T / T0) ** 1.5 * (T0 + S) / (T + S)


def mean_free_path(T: np.ndarray | float, P: np.ndarray | float) -> np.ndarray | float:
    """Gas mean free path lambda = k_B T / (sqrt(2) pi d_g^2 P) [m]."""
    T = np.asarray(T, dtype=float)
    P = np.asarray(P, dtype=float)
    return C.K_B * T / (np.sqrt(2.0) * np.pi * C.D_GAS**2 * P)


def eddy_diffusivity(z: np.ndarray | float) -> np.ndarray | float:
    """Eddy diffusion coefficient K(z) [m^2/s] (placeholder Titan profile).

    Small in the lower stratosphere, rising steeply toward the homopause.  Order
    of magnitude only; replace with a Vuitton (2019)-type profile when needed.
    """
    z = np.asarray(z, dtype=float)
    # ~4e2 m^2/s low down, ramping to ~1e5 m^2/s above ~500 km
    return 4.0e2 * np.exp(z / 90e3)


class Atmosphere:
    """Container for the background profiles on a fixed altitude grid.

    Holds T(z) and P(z) (P from hydrostatic integration of T) and exposes
    interpolated T, P, eta, lambda, g, K at arbitrary altitude.
    """

    def __init__(self, z: np.ndarray, T: np.ndarray, P_surf: float = 1.47e5,
                 molar_mass: float = C.M_N2):
        order = np.argsort(z)
        self.z = np.asarray(z, dtype=float)[order]
        self.T = np.asarray(T, dtype=float)[order]
        self.P_surf = float(P_surf)
        self.molar_mass = float(molar_mass)
        self.P = self._hydrostatic_pressure()

    @classmethod
    def titan_reference(cls, n: int = 200, z_top: float = 540e3) -> "Atmosphere":
        """Default crude Titan column on a uniform grid 0..z_top."""
        z = np.linspace(0.0, z_top, n)
        T = np.interp(z, _Z_REF, _T_REF)
        return cls(z, T)

    @classmethod
    def from_profile(cls, z: np.ndarray, T: np.ndarray, P_surf: float = 1.47e5) -> "Atmosphere":
        """Build from an externally supplied (z, T) profile, e.g. DISORT output."""
        return cls(z, T, P_surf=P_surf)

    def _hydrostatic_pressure(self) -> np.ndarray:
        """Integrate dP/dz = -P M g(z) / (R T) upward from the surface."""
        z, T = self.z, self.T
        g = np.asarray(gravity(z), dtype=float)
        # d(lnP)/dz = -M g / (R T); trapezoidal cumulative integral.
        integrand = -self.molar_mass * g / (C.R_GAS * T)
        dz = np.diff(z)
        avg = 0.5 * (integrand[1:] + integrand[:-1])
        ln_ratio = np.concatenate([[0.0], np.cumsum(avg * dz)])
        return self.P_surf * np.exp(ln_ratio)

    # --- interpolated accessors ---
    def temperature(self, z):
        return np.interp(z, self.z, self.T)

    def pressure(self, z):
        # interpolate in log-pressure for smoothness
        return np.exp(np.interp(z, self.z, np.log(self.P)))

    def viscosity(self, z):
        return viscosity_N2(self.temperature(z))

    def mfp(self, z):
        return mean_free_path(self.temperature(z), self.pressure(z))

    def gravity(self, z):
        return gravity(z)

    def K(self, z):
        return eddy_diffusivity(z)
