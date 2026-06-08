"""Shortwave heating, longwave cooling, and radiative-equilibrium temperature.

Couples the layered :class:`~rt.column.Column`, the haze microphysics profile,
and the two-band DISORT driver into heating/cooling rates and a relaxation
solver for the radiative-equilibrium temperature profile.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .column import Column
from .optics import OpticsParams, shortwave_optics, longwave_optics_spectral
from .cia import CIABands, Composition
from . import disort_driver as dd

# Solar constant at Titan (~9.58 AU): 1361 / 9.58^2 W/m^2
S0_TITAN = 1361.0 / 9.58**2          # ~14.8 W/m^2


@dataclass
class SolarForcing:
    fbeam: float = S0_TITAN          # beam flux on a surface normal to the beam
    umu0: float = 0.35               # cos(solar zenith); ~ Huygens-site average
    albedo: float = 0.1              # surface albedo


@dataclass
class Fluxes:
    z_lev: np.ndarray
    sw_net: np.ndarray               # shortwave net downward flux at levels [W/m^2]
    lw_net: np.ndarray               # longwave net downward flux at levels [W/m^2]
    sw_heating: np.ndarray           # per-layer [K/Titan-day]
    lw_heating: np.ndarray           # per-layer (negative = cooling)
    net_heating: np.ndarray          # per-layer
    z_mid: np.ndarray


def compute_fluxes(column: Column, micro, op: OpticsParams | None = None,
                   solar: SolarForcing | None = None, nstr: int = 8,
                   cia: CIABands | None = None, comp: Composition | None = None) -> Fluxes:
    """Run both bands once and return fluxes + heating rates.

    The longwave is spectral, with collision-induced absorption (CIA) as the
    explicit gas opacity (N2-N2, N2-CH4, CH4-CH4, N2-H2).
    """
    op = op or OpticsParams()
    solar = solar or SolarForcing()
    cia = cia or CIABands()

    tau_sw, ssa_sw, g_sw = shortwave_optics(column, micro, op)
    tau_lw, ssa_lw, g_lw = longwave_optics_spectral(column, micro, cia, op, comp)

    sw_net = dd.solve_shortwave(tau_sw, ssa_sw, g_sw,
                                fbeam=solar.fbeam, umu0=solar.umu0,
                                albedo=solar.albedo, nstr=nstr)
    lw_net = dd.solve_longwave_spectral(tau_lw, ssa_lw, g_lw, column.T,
                                        cia.band_lo, cia.band_hi, albedo=0.0,
                                        nstr=nstr)

    # per-layer deposited energy = net-downward-flux convergence (ascending)
    sw_dF = np.diff(sw_net)
    lw_dF = np.diff(lw_net)
    sw_h = column.heating_rate(sw_dF)
    lw_h = column.heating_rate(lw_dF)
    return Fluxes(column.z, sw_net, lw_net, sw_h, lw_h, sw_h + lw_h, column.z_mid)


def _levels_from_layers(column: Column, T_lyr, T_surf):
    """Reconstruct level temperatures from layer-mean temperatures.

    Monotonic interpolation of the layer-centred temperatures (at z_mid) onto the
    levels (at z), with the surface level pinned to ``T_surf``.  Because this map
    is interpolation (smoothing), not layer<->level aliasing, it does not seed the
    checkerboard mode that a layer->level->layer round trip does.
    """
    T_lev = np.interp(column.z, column.z_mid, T_lyr)
    T_lev[0] = T_surf
    return T_lev


def radiative_equilibrium(column: Column, micro, op: OpticsParams | None = None,
                          solar: SolarForcing | None = None, nstr: int = 8,
                          n_iter: int = 500, dT_max: float = 1.0,
                          dt_days: float = 0.01, relax: float = 0.15,
                          fix_surface: bool = True, tol: float = 0.2,
                          verbose: bool = False):
    """Relax toward radiative equilibrium in LAYER-mean temperature.

    The state is the layer-mean temperature ``T_lyr`` (one value per layer), which
    is in 1:1 correspondence with the per-layer net heating rate -- so the update
    has no layer<->level aliasing and does not zig-zag.  Each iteration:
    reconstruct level temperatures (for the Planck source), run both bands, and
    nudge ``T_lyr`` along the net heating rate (under-relaxed, per-layer clipped at
    ``dT_max``).  The surface temperature is held fixed.

    Returns (column_eq, history) where history is max |net heating| in the
    resolved bulk [K/Titan-day] per iteration.
    """
    op = op or OpticsParams()
    solar = solar or SolarForcing()
    cia = CIABands()                       # build the CIA table once
    T_surf = column.T[0]
    T_lyr = column.T_mid.copy()            # layer-mean temperature is the state
    history = []
    col = column

    for it in range(n_iter):
        T_lev = _levels_from_layers(column, T_lyr, T_surf)
        col = Column(column.z, T_lev, column.P, column.g)
        fx = compute_fluxes(col, micro, op, solar, nstr=nstr, cia=cia)
        # residual measured in the resolved bulk (exclude the near-massless top)
        bulk = col.P_mid > 1.0         # Pa
        resid = np.max(np.abs(fx.net_heating[bulk])) if bulk.any() else \
            np.max(np.abs(fx.net_heating))
        history.append(resid)
        if verbose and (it % 25 == 0 or it == n_iter - 1):
            print(f"  iter {it:4d}  max|net heating|(bulk) = {resid:8.3f} K/Titan-day")
        if resid < tol:
            break
        # per-layer clipped, under-relaxed step (no level<->layer aliasing)
        dT = relax * np.clip(fx.net_heating * dt_days, -dT_max, dT_max)
        T_lyr = np.clip(T_lyr + dT, 40.0, 400.0)

    T_lev = _levels_from_layers(column, T_lyr, T_surf) if fix_surface else col.T
    return Column(column.z, T_lev, column.P, column.g), np.array(history)
