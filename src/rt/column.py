"""Layered 1-D column for radiative transfer.

Builds DISORT-ready level/layer structure from a microphysics
:class:`~microphysics.atmosphere.Atmosphere`.  DISORT orders levels top-down
(index 0 = top of atmosphere, index nlyr = surface); this class stores levels
bottom-up (ascending altitude, matching the microphysics convention) and exposes
helpers to hand top-down arrays to the driver.
"""

from __future__ import annotations

import numpy as np

# N2-dominated atmosphere thermodynamics
CP = 1040.0          # specific heat at constant pressure [J/kg/K] (N2)
M_AIR = 0.0280       # mean molar mass [kg/mol] (N2 + ~few % CH4)
R_GAS = 8.314462618
TITAN_DAY = 15.945 * 86400.0   # seconds in a Titan solar day


class Column:
    """Plane-parallel layered column (levels ascending in altitude)."""

    def __init__(self, z_lev, T_lev, P_lev, g_lev):
        self.z = np.asarray(z_lev, float)       # level altitudes [m], ascending
        self.T = np.asarray(T_lev, float)       # level temperatures [K]
        self.P = np.asarray(P_lev, float)       # level pressures [Pa]
        self.g = np.asarray(g_lev, float)       # level gravity [m/s^2]
        self.nlvl = self.z.size
        self.nlyr = self.nlvl - 1

    @classmethod
    def from_atmosphere(cls, atm, nlyr: int = 50, z_top: float = 500e3):
        """Build a column on a pressure-spaced grid from an Atmosphere."""
        # levels equally spaced in log-pressure from surface to z_top
        P_surf = float(atm.pressure(0.0))
        P_top = float(atm.pressure(z_top))
        P_lev = np.logspace(np.log10(P_surf), np.log10(P_top), nlyr + 1)
        # invert P(z) by interpolation on the atmosphere grid
        z_grid = atm.z
        P_grid = atm.pressure(z_grid)
        # P decreases with z: interpolate z as a function of log P (ascending logP)
        order = np.argsort(np.log(P_grid))
        z_lev = np.interp(np.log(P_lev), np.log(P_grid)[order], z_grid[order])
        T_lev = atm.temperature(z_lev)
        g_lev = atm.gravity(z_lev)
        return cls(z_lev, T_lev, P_lev, g_lev)

    # --- layer-centred quantities (size nlyr) ---
    @property
    def dz(self):
        """Layer geometric thickness [m] (positive)."""
        return np.abs(np.diff(self.z))

    @property
    def dP(self):
        """Layer pressure thickness [Pa] (positive, = P_lower - P_upper)."""
        return np.abs(np.diff(self.P))

    @property
    def z_mid(self):
        return 0.5 * (self.z[1:] + self.z[:-1])

    @property
    def T_mid(self):
        return 0.5 * (self.T[1:] + self.T[:-1])

    @property
    def P_mid(self):
        return np.sqrt(self.P[1:] * self.P[:-1])

    @property
    def g_mid(self):
        return 0.5 * (self.g[1:] + self.g[:-1])

    def number_column(self):
        """Gas number column per layer [molecules/m^2] from hydrostatic dP/g."""
        # mass column = dP/g; number column = mass_col * N_A / M_air
        from microphysics.constants import N_A
        return (self.dP / self.g_mid) * N_A / M_AIR

    # --- DISORT ordering helpers (top-down) ---
    def T_levels_topdown(self):
        return self.T[::-1].copy()

    def heating_rate(self, dFnet_down):
        """Convert a per-layer net-downward-flux convergence to a heating rate.

        ``dFnet_down`` is F_net_down(top of layer) - F_net_down(bottom), the
        energy deposited per unit area in the layer [W/m^2].  Returns K per
        Titan day, layer-centred (ascending order).
        """
        dFnet_down = np.asarray(dFnet_down, float)
        # mass per unit area in layer = dP/g
        mass_col = self.dP / self.g_mid           # [kg/m^2]
        rate_per_s = dFnet_down / (mass_col * CP)  # [K/s]
        return rate_per_s * TITAN_DAY
