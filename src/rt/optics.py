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
from functools import lru_cache
from pathlib import Path

import numpy as np

_DATA = (Path(__file__).resolve().parents[2]
         / "src" / "example_bowen_fort" / "INPUT" / "DATA")


def _parse_preschaze(field: str):
    """Read preschazev_<field>.txt -> (wn[nw], pl[npl], data[npl, nw])."""
    L = (_DATA / f"preschazev_{field}.txt").read_text().split("\n")
    nw, npl = int(L[1]), int(L[2])
    wn = np.array(L[3].split(), float)
    pl = np.array(L[4].split(), float)
    data = np.array([L[5 + i].split() for i in range(npl)], float)
    return wn, pl, data


def prescribed_haze_layer_tau(column, band_centers, kind="v"):
    """Per-layer haze optical depth from the observational prescribed haze.

    The preschaze<kind>_tau.txt files give the *cumulative* haze optical depth
    (top-down) on a (wavenumber, pressure) grid; we interpolate to the model band
    centres and level pressures, then difference to per-layer tau[nband, nlyr].
    ``kind`` is 'v' (visible/SW) or 'i' (IR/LW).
    """
    wn, pl, data = _parse_preschaze("tau") if kind == "v" else _parse_preschaze_i("tau")
    bc = np.asarray(band_centers, float)
    Plev = column.P
    o = np.argsort(pl)
    # cumulative tau(band, pl): interpolate the (npl, nw) table to band centres
    tau_bp = np.vstack([np.interp(bc, wn, data[ip]) for ip in range(pl.size)]).T  # (nband, npl)
    # interpolate to the model level pressures
    cum = np.vstack([np.interp(Plev, pl[o], tau_bp[b][o]) for b in range(bc.size)])  # (nband, nlvl)
    # per-layer = difference of cumulative tau between bounding levels
    return np.maximum(np.abs(np.diff(cum, axis=1)), 0.0)                            # (nband, nlyr)


def _parse_preschaze_i(field: str):
    L = (_DATA / f"preschazei_{field}.txt").read_text().split("\n")
    nw, npl = int(L[1]), int(L[2])
    wn = np.array(L[3].split(), float)
    pl = np.array(L[4].split(), float)
    data = np.array([L[5 + i].split() for i in range(npl)], float)
    return wn, pl, data


@lru_cache(maxsize=4)
def _rayleigh_table():
    d = np.loadtxt(_DATA / "Rayleigh.txt", skiprows=2)
    return d[:, 0], d[:, 1]


def rayleigh_band_tau(column, band_centers):
    """Rayleigh optical depth tau[nband, nlyr] (pure scattering).

    Fortran form tau_ray(k) = dp * tauray(band), with tauray ~ nu^4; using the
    layer pressure thickness in mbar gives a physical Titan column tau (~0.1 in
    the visible, rising steeply into the UV).
    """
    wn, coef = _rayleigh_table()
    tr = np.interp(np.asarray(band_centers, float), wn, coef)   # per band
    dp_mbar = column.dP / 100.0                                 # Pa -> mbar
    return np.maximum(tr[:, None] * dp_mbar[None, :], 0.0)      # (nband, nlyr)


@lru_cache(maxsize=8)
def spectral_haze_sw(band_centers_key):
    """Observational haze single-scattering albedo and asymmetry per SW band.

    From the Huygens/DISR-derived prescribed haze (Doose 2016 / Tomasko 2008):
    omega0(lambda) and g(lambda), interpolated to the model band centres.  These
    properties are nearly pressure-independent, so a representative stratospheric
    level is used.  omega0 falls to ~0 (pure absorber) in the UV -- the dominant
    stratospheric shortwave heating that a gray albedo misses.
    """
    band_centers = np.array(band_centers_key, float)
    wn_s, pl, ssa = _parse_preschaze("ssa")
    wn_g, _, gg = _parse_preschaze("g")
    ip = int(np.argmin(np.abs(pl - 200.0)))        # ~200 Pa (stratosphere)
    omega0 = np.interp(band_centers, wn_s, ssa[ip])
    g = np.interp(band_centers, wn_g, gg[ip])
    return np.clip(omega0, 0.0, 0.999999), np.clip(g, 0.0, 0.95)


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
    prescribed_haze: bool = False  # back-compat: True => haze_mode='prescribed'
    # haze cross-section model: 'rdg' (mean-field aggregate optics, physical),
    # 'gray' (Q_ext pi r_a^2), or 'prescribed' (observational optical depths)
    haze_mode: str = "rdg"
    # scale on the RDG monomer absorption, anchoring the (uncertain) tholin
    # absorption to the observed haze opacity (Titan tholins absorb more than
    # the Khare lab values); 3.0 -> visible column tau ~ Doose's ~8
    haze_abs_scale: float = 3.0
    # match the reference Fortran RT: warm LW top boundary, LW surface reflection
    # (eps=0.95), SW surface albedo 0.15, Rayleigh scattering in the shortwave
    match_fortran: bool = True


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


def haze_band_tau(column, micro, p: OpticsParams, band_centers, kind,
                  omega0_band=None):
    """Per-band haze optical depth tau[nband, nlyr] for the selected haze model.

    ``kind`` is 'v' (shortwave) or 'i' (longwave).  Modes: 'rdg' (mean-field
    aggregate optics from the microphysics), 'gray' (Q_ext pi r_a^2), 'prescribed'
    (observational optical depths).
    """
    nb = len(band_centers)
    mode = "prescribed" if p.prescribed_haze else p.haze_mode
    if mode == "prescribed":
        return prescribed_haze_layer_tau(column, band_centers, kind=kind)
    if mode == "rdg":
        from .aggregate_optics import aggregate_haze_layer_tau
        d = micro.params.d_mono
        return p.haze_abs_scale * aggregate_haze_layer_tau(
            column, micro, band_centers, omega0_band, d, pure_absorber=(kind == "i"))
    return np.broadcast_to(layer_haze_tau(column, micro, p), (nb, column.nlyr))


def longwave_optics_ck(column, micro, cia, ck_lw, p: OpticsParams, comp=None):
    """Correlated-k longwave optics: per (band, Gauss-point) wave.

    Sums the gas correlated-k optical depth (CH4/C2H2/C2H6/C2H4/HCN rotational
    lines, the stratospheric coolants), the collision-induced absorption, and the
    (gray, pure-absorber) IR haze.  Everything is absorption, so ssa=0.  Returns
    (tau, ssa, band_lo, band_hi, weights) with tau/ssa (nwave, nlyr) ascending,
    band edges repeated per Gauss-point, weights the Gauss weights.  CIA and the
    gas k-table share the same band grid.
    """
    tau_gas = ck_lw.layer_tau(column)                 # (nband, ngauss, nlyr)
    tau_cia = cia.optical_depth(column, comp)         # (nband, nlyr)
    nb, ng, nlyr = tau_gas.shape
    haze_band = haze_band_tau(column, micro, p, ck_lw.bands, "i")   # (nb, nlyr)
    nwave = nb * ng
    tau = np.empty((nwave, nlyr))
    for b in range(nb):
        base = tau_cia[b] + haze_band[b]              # (nlyr,) absorbers shared over g
        for gi in range(ng):
            tau[b * ng + gi] = tau_gas[b, gi] + base
    ssa = np.zeros((nwave, nlyr))
    band_lo = np.repeat(ck_lw.band_lo, ng)
    band_hi = np.repeat(ck_lw.band_hi, ng)
    weights = np.tile(ck_lw.gw, nb)
    return tau, ssa, band_lo, band_hi, weights


def shortwave_optics_ck(column, micro, ck, p: OpticsParams, comp=None):
    """Correlated-k shortwave optics: per (band, Gauss-point) wave.

    Combines spectral haze single-scattering properties (omega0(lambda),
    g(lambda) from the observational prescribed haze) with spectral CH4
    absorption from ``rt.correlated_k.CorrelatedKSW``.  Returns
    (tau, ssa, g_wave, fbeam, weights): tau/ssa are (nwave, nlyr) ascending,
    g_wave[nwave] the per-wave asymmetry, fbeam[nwave] the per-band TOA solar
    flux, weights[nwave] the Gauss weights.  The gas is a pure absorber, so it
    lowers ``ssa`` without changing the (haze) phase function.
    """
    tau_gas = ck.layer_tau(column, comp)              # (nband, ngauss, nlyr)
    nb, ng, nlyr = tau_gas.shape
    omega0_b, g_b = spectral_haze_sw(tuple(np.round(ck.bands, 3)))   # per band
    haze_band = haze_band_tau(column, micro, p, ck.bands, "v", omega0_b)  # (nb, nlyr)
    # Rayleigh: pure-scattering molecular optical depth (Fortran has it, ~nu^4).
    tau_ray = (rayleigh_band_tau(column, ck.bands) if p.match_fortran
               else np.zeros((nb, nlyr)))              # (nb, nlyr)
    nwave = nb * ng
    tau = np.empty((nwave, nlyr))
    ssa = np.empty((nwave, nlyr))
    g_wave = np.empty(nwave)
    for b in range(nb):
        tau_haze = haze_band[b]                        # (nlyr,)
        # scattering = haze (g=g_b) + Rayleigh (g=0); keep the per-wave g as the
        # scattering-weighted asymmetry (Rayleigh dilutes the haze forward peak)
        sca_haze = omega0_b[b] * tau_haze             # haze scattering opt depth (nlyr)
        sca = sca_haze + tau_ray[b]                    # total scattering (nlyr)
        sca_tot = float(sca.sum())
        g_wave_b = (g_b[b] * float(sca_haze.sum()) / sca_tot) if sca_tot > 0 else g_b[b]
        for gi in range(ng):
            w = b * ng + gi
            tt = tau_haze + tau_ray[b] + tau_gas[b, gi]
            tau[w] = tt
            ssa[w] = np.clip(np.where(tt > 0, sca / np.maximum(tt, 1e-300), 0.0),
                             0.0, 0.999999)
            g_wave[w] = g_wave_b
    fbeam = np.repeat(ck.solar, ng)                   # per-band solar, repeated over g
    weights = np.tile(ck.gw, nb)                      # Gauss weight per wave
    return tau, ssa, g_wave, fbeam, weights


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
