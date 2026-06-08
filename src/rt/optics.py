"""Per-layer optical properties for the two-band (shortwave / longwave) RT.

Two contributions are summed per layer:

* **Haze** -- from the microphysics profile (number density and mobility radius).
  Extinction uses a geometric cross-section ``Q_ext * pi * r_a^2`` (Q_ext ~ 2);
  single-scattering albedo and asymmetry are taken from the Huygens/DISR-based
  literature (Doose 2016; Tomasko 2008): scattering + absorbing in the
  shortwave, near-pure absorber in the thermal IR.

* **Gas** -- a deliberately simple GRAY placeholder standing in for the
  correlated-k CH4 + collision-induced-absorption opacity of the real model
  (Lombardo 2023).  Shortwave gas opacity is distributed by mass column (CH4-like
  absorption); longwave gas opacity scales as number-density squared (CIA-like).
  Column totals are tunable and currently set to give Tomasko-2008-like absorbed
  fractions.  THIS IS THE PRINCIPAL APPROXIMATION in the present RT and the first
  thing to replace with real k-tables.

Layer optical properties are returned as (tau, ssa, g) arrays, ascending in
altitude; the driver reverses them to DISORT top-down order.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass
class OpticsParams:
    # haze single-scattering properties
    Qext: float = 2.0            # geometric extinction efficiency
    ssa_haze_sw: float = 0.92    # shortwave single-scattering albedo (Doose 2016)
    g_haze_sw: float = 0.60      # shortwave asymmetry parameter
    ssa_haze_lw: float = 0.0     # IR: pure absorber (Tomasko 2008)
    g_haze_lw: float = 0.0
    # gray-gas placeholders (column totals)
    tau_gas_sw: float = 2.0      # total shortwave gas (CH4-like) optical depth
    ssa_gas_sw: float = 0.0      # gas absorbs in SW (Rayleigh scattering ignored)
    tau_gas_lw_gray: float = 0.0  # optional extra gray LW continuum (lines); CIA is explicit


def haze_extinction_per_length(n, r_a, p: OpticsParams):
    """Haze extinction coefficient [1/m] = n * Q_ext * pi * r_a^2."""
    return n * p.Qext * np.pi * r_a**2


def layer_haze_tau(column, micro, p: OpticsParams):
    """Haze optical depth per layer, interpolated onto the column layers."""
    # interpolate microphysics (number density, mobility radius) to layer mids
    zc = column.z_mid
    # micro arrays may be descending (master) or ascending (bvp); sort
    order = np.argsort(micro.z)
    zz = micro.z[order]
    n = np.interp(zc, zz, micro.n[order])
    r_a = np.interp(zc, zz, micro.r_a[order])
    # no haze above the microphysics domain (avoid flat-extrapolating the top)
    n = np.where(zc > zz[-1], 0.0, n)
    kext = haze_extinction_per_length(n, r_a, p)
    return kext * column.dz


def layer_gas_tau_sw(column, p: OpticsParams):
    """Shortwave gray gas optical depth per layer, distributed by mass column."""
    mass_col = column.dP / column.g_mid          # [kg/m^2] per layer
    frac = mass_col / mass_col.sum()
    return p.tau_gas_sw * frac




def combine(tau_a, ssa_a, g_a, tau_b, ssa_b, g_b):
    """Mix two species' (tau, ssa, g) into layer-effective values."""
    tau = tau_a + tau_b
    sca = ssa_a * tau_a + ssa_b * tau_b           # scattering optical depth
    ssa = np.where(tau > 0, sca / np.maximum(tau, 1e-300), 0.0)
    gnum = g_a * ssa_a * tau_a + g_b * ssa_b * tau_b
    g = np.where(sca > 0, gnum / np.maximum(sca, 1e-300), 0.0)
    return tau, np.clip(ssa, 0.0, 0.999999), np.clip(g, 0.0, 0.95)


def shortwave_optics(column, micro, p: OpticsParams):
    """Return (tau, ssa, g) per layer for the shortwave band (ascending)."""
    th = layer_haze_tau(column, micro, p)
    tg = layer_gas_tau_sw(column, p)
    return combine(th, p.ssa_haze_sw, p.g_haze_sw, tg, p.ssa_gas_sw, 0.0)


def longwave_optics_spectral(column, micro, cia, p: OpticsParams, comp=None):
    """Per-band longwave (tau, ssa, g) arrays, shape (nband, nlyr), ascending.

    Opacity = collision-induced absorption (explicit, per band, from
    :class:`rt.cia.CIABands`) + haze (gray, near-pure absorber) + an optional
    extra gray continuum.  Both gas CIA and the IR haze are pure absorbers, so
    the single-scattering albedo and asymmetry are zero everywhere.
    """
    tau_cia = cia.optical_depth(column, comp)            # (nband, nlyr)
    tau_haze = layer_haze_tau(column, micro, p)          # (nlyr,)
    # optional extra gray LW continuum (e.g. unmodeled rotational lines)
    if p.tau_gas_lw_gray > 0:
        mass = column.dP / column.g_mid
        tau_gray = p.tau_gas_lw_gray * mass / mass.sum()
    else:
        tau_gray = 0.0
    tau = tau_cia + (tau_haze + tau_gray)[None, :]
    ssa = np.zeros_like(tau)
    g = np.zeros_like(tau)
    return tau, ssa, g
