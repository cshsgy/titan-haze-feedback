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


def solve_shortwave_spectral(tau, ssa, g_haze, fbeam, weights, umu0,
                             albedo=0.1, nstr=8):
    """Multi-wave (correlated-k) solar-beam DISORT solve.

    ``tau, ssa`` are (nwave, nlyr) ascending-in-altitude, one wave per
    (band, Gauss-point); ``fbeam[nwave]`` is the per-band TOA solar flux and
    ``weights[nwave]`` the Gauss weights.  Scattering is haze-only, so the phase
    function is a single Henyey-Greenstein with asymmetry ``g_haze`` for every
    wave (gas is a pure absorber, varying only ``ssa``).  Returns the
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
    prop[:, 0, :, 2:] = scattering_moments(nstr, "henyey-greenstein", float(g_haze))

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
                            albedo=0.0, nstr=8):
    """Multi-band thermal DISORT solve with the Planck source.

    ``tau, ssa, g`` are (nband, nlyr) ascending-in-altitude.  Each band uses its
    own ``[band_lo, band_hi]`` wavenumber limits for the band-integrated Planck
    source.  Gas CIA and IR haze are pure absorbers, so ``ssa = g = 0`` and the
    phase-function moments are left at zero (isotropic, irrelevant when ssa=0).
    Returns the band-summed net downward flux at levels (ascending); thermal
    cooling is the divergence of the net upward flux.
    """
    tau = np.asarray(tau, float)
    nband, nlyr = tau.shape
    band_lo = np.asarray(band_lo, float)
    band_hi = np.asarray(band_hi, float)

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
    # Top boundary = cold space: emit at ~0 K (no spurious downwelling thermal).
    # btemp = surface temperature; the atmosphere radiates upward to space.
    ds.forward(prop, temf=temf,
               btemp=torch.tensor([float(T_levels_ascending[0])]),
               ttemp=torch.tensor([2.7]),
               temis=torch.zeros((nband, 1)),
               albedo=torch.full((nband, 1), float(albedo)))
    f = ds.gather_flx()[:, 0].numpy()                 # (nband, nlvl, 8), top-down
    net_td = (f[:, :, pydisort.kIRFLDIR] + f[:, :, pydisort.kIFLDN]
              - f[:, :, pydisort.kIFLUP])
    return net_td.sum(axis=0)[::-1].copy()            # band-summed, ascending
