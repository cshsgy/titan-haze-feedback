"""Thin wrapper around :mod:`pydisort` for single-band flux calculations.

Handles the bookkeeping: float64 tensors, top-down layer ordering, Henyey-
Greenstein phase-function moments, and extraction of the net downward flux at
levels.  Two entry points: :func:`solve_shortwave` (collimated solar beam) and
:func:`solve_longwave` (thermal emission via the Planck source).
"""

from __future__ import annotations

import numpy as np
import torch

import pydisort
from pydisort import DisortOptions, Disort, scattering_moments

torch.set_default_dtype(torch.float64)

# damping on the above-model-top optical depth for the warm LW top boundary,
# tuned so the top-layer cooling matches the reference model
LW_TOP_SCALE = 0.02


def _build(nlyr, nstr, flags, wave_lower=None, wave_upper=None):
    op = DisortOptions().flags(flags).nwave(1).ncol(1)
    op.ds().nlyr = nlyr
    op.ds().nstr = nstr
    op.ds().nmom = nstr
    op.ds().nphase = nstr
    if wave_lower is not None:
        op.wave_lower([float(wave_lower)])
        op.wave_upper([float(wave_upper)])
    return Disort(op)


def _prop_tensor(tau, ssa, g, nstr):
    """Assemble the (nwave=1, ncol=1, nlyr, 2+nmom) optical-property tensor.

    Layers are passed ascending (surface..top) and reversed here to DISORT
    top-down order.
    """
    tau = np.asarray(tau, float)[::-1]
    ssa = np.asarray(ssa, float)[::-1]
    g = np.asarray(g, float)[::-1]
    nlyr = tau.size
    nprop = 2 + nstr
    prop = torch.zeros((1, 1, nlyr, nprop))
    prop[0, 0, :, 0] = torch.from_numpy(tau.copy())
    prop[0, 0, :, 1] = torch.from_numpy(ssa.copy())
    for i in range(nlyr):
        prop[0, 0, i, 2:] = scattering_moments(nstr, "henyey-greenstein", float(g[i]))
    return prop


def _net_down(ds):
    """Net downward flux at each level (top-down), [W/m^2]."""
    f = ds.gather_flx()[0, 0].numpy()
    return f[:, pydisort.kIRFLDIR] + f[:, pydisort.kIFLDN] - f[:, pydisort.kIFLUP]


def solve_shortwave(tau, ssa, g, fbeam, umu0, albedo=0.1, nstr=8):
    """Solar-beam DISORT solve. Returns net downward flux at levels (ascending)."""
    ds = _build(tau.size, nstr, "onlyfl,lamber")
    prop = _prop_tensor(tau, ssa, g, nstr)
    ds.forward(prop,
               umu0=torch.tensor([float(umu0)]),
               phi0=torch.tensor([0.0]),
               fbeam=torch.tensor([[float(fbeam)]]),
               albedo=torch.tensor([[float(albedo)]]),
               fisot=torch.tensor([[0.0]]))
    return _net_down(ds)[::-1].copy()       # back to ascending


def solve_shortwave_spectral(tau, ssa, g_wave, fbeam, weights, umu0,
                             albedo=0.1, nstr=8):
    """Multi-wave (correlated-k) solar-beam DISORT solve.

    ``tau, ssa`` are (nwave, nlyr) ascending-in-altitude, one wave per
    (band, Gauss-point); ``g_wave[nwave]`` is the per-wave (haze) asymmetry,
    ``fbeam[nwave]`` the per-band TOA solar flux, ``weights[nwave]`` the Gauss
    weights.  Scattering is haze-only (Henyey-Greenstein, asymmetry varying by
    band); the gas is a pure absorber, lowering only ``ssa``.  Returns the
    weight-summed net downward flux at levels (ascending).
    """
    tau = np.asarray(tau, float)
    nwave, nlyr = tau.shape
    op = DisortOptions().flags("onlyfl,lamber").nwave(nwave).ncol(1)
    op.ds().nlyr = nlyr
    op.ds().nstr = nstr
    op.ds().nmom = nstr
    op.ds().nphase = nstr
    ds = Disort(op)

    nprop = 2 + nstr
    prop = torch.zeros((nwave, 1, nlyr, nprop))
    prop[:, 0, :, 0] = torch.from_numpy(tau[:, ::-1].copy())
    prop[:, 0, :, 1] = torch.from_numpy(np.asarray(ssa, float)[:, ::-1].copy())
    # Henyey-Greenstein moments per wave: g_w^l for l = 1..nstr
    g_wave = np.asarray(g_wave, float)
    moments = g_wave[:, None] ** np.arange(1, nstr + 1)[None, :]   # (nwave, nstr)
    prop[:, 0, :, 2:] = torch.from_numpy(moments.copy())[:, None, :]

    fbeam = np.asarray(fbeam, float).reshape(nwave, 1)
    ds.forward(prop,
               umu0=torch.tensor([float(umu0)]),
               phi0=torch.tensor([0.0]),
               fbeam=torch.from_numpy(fbeam.copy()),
               albedo=torch.full((nwave, 1), float(albedo)),
               fisot=torch.zeros((nwave, 1)))
    f = ds.gather_flx()[:, 0].numpy()                 # (nwave, nlvl, 8), top-down
    net_td = (f[:, :, pydisort.kIRFLDIR] + f[:, :, pydisort.kIFLDN]
              - f[:, :, pydisort.kIFLUP])
    w = np.asarray(weights, float)[:, None]
    return (net_td * w).sum(axis=0)[::-1].copy()      # weight-summed, ascending


def solve_longwave_spectral(tau, ssa, g, T_levels_ascending, band_lo, band_hi,
                            albedo=0.0, nstr=8, weights=None, P_levels_ascending=None):
    """Multi-wave thermal DISORT solve with the Planck source.

    ``tau, ssa, g`` are (nwave, nlyr) ascending-in-altitude; each wave has its
    own ``[band_lo, band_hi]`` Planck limits.  For correlated-k a wave is a
    (band, Gauss-point) pair: pass ``band_lo/band_hi`` repeated per Gauss-point
    and ``weights`` the Gauss weights (default: all 1, i.e. a pure band sum, as
    for CIA-only).  Gas + IR haze are pure absorbers (ssa=0, isotropic).  Returns
    the weight-summed net downward flux at levels (ascending); thermal cooling is
    the divergence of the net upward flux.
    """
    tau = np.asarray(tau, float)
    nband, nlyr = tau.shape
    band_lo = np.asarray(band_lo, float)
    band_hi = np.asarray(band_hi, float)
    w = np.ones(nband) if weights is None else np.asarray(weights, float)

    op = DisortOptions().flags("onlyfl,lamber,planck").nwave(nband).ncol(1)
    op.ds().nlyr = nlyr
    op.ds().nstr = nstr
    op.ds().nmom = nstr
    op.ds().nphase = nstr
    op.wave_lower(band_lo)
    op.wave_upper(band_hi)
    ds = Disort(op)

    nprop = 2 + nstr
    prop = torch.zeros((nband, 1, nlyr, nprop))
    prop[:, 0, :, 0] = torch.from_numpy(tau[:, ::-1].copy())          # top-down
    prop[:, 0, :, 1] = torch.from_numpy(np.asarray(ssa, float)[:, ::-1].copy())

    T_td = np.asarray(T_levels_ascending, float)[::-1]
    temf = torch.from_numpy(T_td.copy()).reshape(1, -1)
    # Top boundary.  Cold space (ttemp=2.7, temis=0) unless P_levels given, in
    # which case mimic the reference model's warm downwelling btop: the
    # atmosphere above the model top emits at the top-level temperature with an
    # effective emissivity 1-exp(-tau_top*ratio/mu_bar), ratio = P_top/(P_2 below
    # the top - P_top), mu_bar=0.5 (the IR diffusivity).
    if P_levels_ascending is None:
        ttemp, temis = 2.7, torch.zeros((nband, 1))
    else:
        Pa = np.asarray(P_levels_ascending, float)
        # effective optical depth of the (truncated) atmosphere above the model
        # top ~ local extinction-per-pressure * pressure above; LW_TOP_SCALE damps
        # it to the reference model's btop magnitude
        dptop = max(Pa[-2] - Pa[-1], 1e-30)
        tau_top = tau[:, -1]                          # top layer optical depth (per wave)
        tau_above = LW_TOP_SCALE * tau_top * (Pa[-1] / dptop)
        eps = 1.0 - np.exp(-tau_above / 0.5)
        ttemp = float(T_levels_ascending[-1])
        temis = torch.from_numpy(np.clip(eps, 0.0, 1.0).reshape(nband, 1).copy())
    ds.forward(prop, temf=temf,
               btemp=torch.tensor([float(T_levels_ascending[0])]),
               ttemp=torch.tensor([float(ttemp)]),
               temis=temis,
               albedo=torch.full((nband, 1), float(albedo)))
    f = ds.gather_flx()[:, 0].numpy()                 # (nwave, nlvl, 8), top-down
    net_td = (f[:, :, pydisort.kIRFLDIR] + f[:, :, pydisort.kIFLDN]
              - f[:, :, pydisort.kIFLUP])
    return (net_td * w[:, None]).sum(axis=0)[::-1].copy()   # weight-summed, ascending
