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
from .cia import CIABands, CIAExpSum, Composition
from . import disort_driver as dd

# Solar constant at Titan (~9.5 AU, reference-model value): 1361 / 9.5^2 W/m^2
S0_TITAN = 1361.0 / 9.5**2           # ~15.1 W/m^2


@dataclass
class SolarForcing:
    fbeam: float = S0_TITAN          # beam flux on a surface normal to the beam
    umu0: float = 1.0 / np.pi        # diurnally-averaged cos(zenith) at lat=0, Ls=0
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
                   cia: CIABands | None = None, comp: Composition | None = None,
                   ck=None, ck_lw=None) -> Fluxes:
    """Run both bands once and return fluxes + heating rates.

    Longwave: collision-induced absorption + (if ``ck_lw`` given) correlated-k
    gas lines (CH4/C2H2/C2H6/C2H4/HCN) + IR haze.  Shortwave: correlated-k CH4
    (``ck``) + spectral haze if ``ck`` given, else the gray placeholder.

    ``ck`` and ``ck_lw`` must be passed together (both or neither): the
    correlated-k shortwave and the gas-line longwave form one energy balance, and
    ``radiative_equilibrium`` always uses both.  Mixing them (e.g. ck shortwave
    with CIA-only longwave) silently produces a profile whose net heating is
    enormous where the model actually converged to ~0 -- the kind of footgun that
    makes a *converged* run look wildly out of balance when re-plotted.  Refuse it.
    """
    if (ck is None) != (ck_lw is None):
        raise ValueError(
            "compute_fluxes: pass `ck` and `ck_lw` together (both or neither). "
            "Correlated-k shortwave with CIA-only longwave (or vice versa) is an "
            "inconsistent energy balance; radiative_equilibrium uses both.")
    op = op or OpticsParams()
    solar = solar or SolarForcing()
    # match-Fortran longwave continuum: the reference's exponential-sum
    # transmission CIA (trans_*.txt) rather than our band-averaged HITRAN table,
    # which overestimates the far-IR continuum by ~1.75x
    if cia is None:
        cia = CIAExpSum() if op.match_fortran else CIABands()

    # match-Fortran longwave boundaries: warm downwelling top + reflecting
    # surface (eps = 0.95); otherwise cold-space top + black surface
    lw_alb = 0.05 if op.match_fortran else 0.0
    Plev = column.P if op.match_fortran else None
    if ck_lw is not None:
        from .optics import longwave_optics_ck
        tau_lw, ssa_lw, blo, bhi, wlw = longwave_optics_ck(column, micro, cia,
                                                           ck_lw, op, comp)
        lw_net = dd.solve_longwave_spectral(tau_lw, ssa_lw, None, column.T,
                                            blo, bhi, albedo=lw_alb, nstr=nstr,
                                            weights=wlw, P_levels_ascending=Plev)
    else:
        tau_lw, ssa_lw, g_lw = longwave_optics_spectral(column, micro, cia, op, comp)
        lw_net = dd.solve_longwave_spectral(tau_lw, ssa_lw, g_lw, column.T,
                                            cia.band_lo, cia.band_hi, albedo=0.0,
                                            nstr=nstr)

    if ck is not None:
        from .optics import shortwave_optics_ck
        # match-Fortran surface: visible albedo 0.15 (albv); else the default
        sw_alb = 0.15 if op.match_fortran else solar.albedo
        tau_sw, ssa_sw, g_haze, fbeam, weights = shortwave_optics_ck(
            column, micro, ck, op, comp)
        sw_net = dd.solve_shortwave_spectral(tau_sw, ssa_sw, g_haze, fbeam,
                                             weights, umu0=solar.umu0,
                                             albedo=sw_alb, nstr=nstr)
    else:
        tau_sw, ssa_sw, g_sw = shortwave_optics(column, micro, op)
        sw_net = dd.solve_shortwave(tau_sw, ssa_sw, g_sw,
                                    fbeam=solar.fbeam, umu0=solar.umu0,
                                    albedo=solar.albedo, nstr=nstr)

    # per-layer deposited energy = net-downward-flux convergence (ascending)
    sw_dF = np.diff(sw_net)
    lw_dF = np.diff(lw_net)
    sw_h = column.heating_rate(sw_dF)
    lw_h = column.heating_rate(lw_dF)
    return Fluxes(column.z, sw_net, lw_net, sw_h, lw_h, sw_h + lw_h, column.z_mid)


def convective_adjust(T_lyr, z_mid, dP, gamma_crit=1.0, T_surf=None, z_surf=0.0,
                      n_sweep=60):
    """Dry/critical-lapse convective adjustment of layer-mean temperatures.

    Where the radiative lapse rate exceeds ``gamma_crit`` [K/km], relax that
    region to a constant critical lapse rate, conserving column enthalpy
    (mass-weighted, weight = dP).  ``gamma_crit`` represents the convective limit;
    for Titan ~1 K/km (sub-dry-adiabatic, methane-moderated).  If ``T_surf`` is
    given, the lowest layer is also checked against the surface.
    """
    T = np.asarray(T_lyr, float).copy()
    z = np.asarray(z_mid, float) / 1e3        # km, ascending
    w = np.asarray(dP, float)                  # mass weight per layer
    n = T.size
    for _ in range(n_sweep):
        # surface-driven instability: warm the lowest layer toward the surface adiabat
        if T_surf is not None:
            t_adia = T_surf - gamma_crit * (z[0] - z_surf / 1e3)
            if T[0] < t_adia - 1e-6:
                T[0] = t_adia
        # adjacent-pair instability: lapse between i (low) and i+1 (high) too steep
        lapse = (T[:-1] - T[1:]) / (z[1:] - z[:-1])   # K/km, >0 = T drops w/ height
        bad = np.where(lapse > gamma_crit + 1e-6)[0]
        if bad.size == 0:
            break
        for i in bad:
            # mix pair (i, i+1) to the critical lapse, conserving enthalpy
            dz = z[i + 1] - z[i]
            wi, wj = w[i], w[i + 1]
            # T_i = a, T_{i+1} = a - gamma*dz ; conserve wi*T_i + wj*T_{i+1}
            a = (wi * T[i] + wj * T[i + 1] + wj * gamma_crit * dz) / (wi + wj)
            T[i] = a
            T[i + 1] = a - gamma_crit * dz
    return T


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
                          n_iter: int = 2000, dT_max: float = 2.0,
                          dt_days: float = 0.02, relax: float = 0.2,
                          gamma_crit: float = 1.0, convective: bool = True,
                          fix_surface: bool = True, tol: float = 0.3,
                          use_ck: bool = True, verbose: bool = False):
    """Relax toward radiative-CONVECTIVE equilibrium in layer-mean temperature.

    The state is the layer-mean temperature ``T_lyr`` (1:1 with the per-layer net
    heating rate, so no layer<->level aliasing / zig-zag).  Each iteration:
    reconstruct level temperatures, run both bands, nudge ``T_lyr`` along the net
    heating rate (under-relaxed, per-layer clipped at ``dT_max``), then apply
    convective adjustment to a critical lapse rate ``gamma_crit`` [K/km].  Without
    the convective step the pure-radiative steady state of this gray-shortwave
    model is nearly isothermal; convection sets the troposphere.  The surface
    temperature is held fixed.

    ``tol`` is on the FULL net-heating residual (incl. the radiative top), so the
    reported convergence reflects top-budget closure, not just the bulk.
    """
    op = op or OpticsParams()
    solar = solar or SolarForcing()
    cia = CIAExpSum() if op.match_fortran else CIABands()   # build the CIA table once
    ck = ck_lw = None
    if use_ck:
        from .correlated_k import CorrelatedKSW, CorrelatedKLW
        ck = CorrelatedKSW()                   # Titan orbit (default sma=9.5 AU)
        ck_lw = CorrelatedKLW()                # IR gas lines (C2H2/HCN/C2H6/...)
    T_surf = column.T[0]
    T_lyr = column.T_mid.copy()            # layer-mean temperature is the state
    z_mid, dP = column.z_mid, column.dP
    history = []
    col = column

    for it in range(n_iter):
        T_lev = _levels_from_layers(column, T_lyr, T_surf)
        col = Column(column.z, T_lev, column.P, column.g)
        fx = compute_fluxes(col, micro, op, solar, nstr=nstr, cia=cia, ck=ck, ck_lw=ck_lw)
        # residual on convectively-stable (radiative) layers -- convective layers
        # are set by the adjustment, not by local radiative balance
        if convective:
            lapse = np.concatenate([[0.0],
                     (T_lyr[:-1] - T_lyr[1:]) / (np.diff(z_mid) / 1e3)])
            radiative = lapse < gamma_crit - 1e-3
        else:
            radiative = np.ones_like(T_lyr, dtype=bool)
        resid = np.max(np.abs(fx.net_heating[radiative])) if radiative.any() \
            else np.max(np.abs(fx.net_heating))
        history.append(resid)
        if verbose and (it % 100 == 0 or it == n_iter - 1):
            print(f"  iter {it:4d}  max|net heating|(radiative) = {resid:8.3f} K/Titan-day")
        if resid < tol:
            break
        # per-layer clipped, under-relaxed radiative step
        dT = relax * np.clip(fx.net_heating * dt_days, -dT_max, dT_max)
        T_lyr = np.clip(T_lyr + dT, 40.0, 400.0)
        # convective adjustment to the critical lapse rate
        if convective:
            T_lyr = convective_adjust(T_lyr, z_mid, dP, gamma_crit, T_surf=T_surf)

    T_lev = _levels_from_layers(column, T_lyr, T_surf)
    eq = Column(column.z, T_lev, column.P, column.g)
    # Return the fluxes the solver ACTUALLY converged with (same ck/ck_lw/op/
    # solar).  Callers must plot/diagnose with these rather than recomputing --
    # recomputing with a different optics combination (e.g. forgetting ck_lw)
    # silently yields a wildly different, un-converged-looking net heating.
    fx = compute_fluxes(eq, micro, op, solar, nstr=nstr, cia=cia, ck=ck, ck_lw=ck_lw)
    return eq, np.array(history), fx
