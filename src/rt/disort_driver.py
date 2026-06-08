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


def solve_longwave(tau, ssa, g, T_levels_ascending, wave_lower=10.0,
                   wave_upper=2000.0, albedo=0.0, nstr=8):
    """Thermal DISORT solve with Planck source.

    ``T_levels_ascending`` are level temperatures (surface..top); the surface and
    top temperatures set ``btemp``/``ttemp``.  Returns net downward flux at levels
    (ascending); thermal **cooling** is the divergence of the net *upward* flux.
    """
    ds = _build(tau.size, nstr, "onlyfl,lamber,planck",
                wave_lower=wave_lower, wave_upper=wave_upper)
    prop = _prop_tensor(tau, ssa, g, nstr)
    T_td = np.asarray(T_levels_ascending, float)[::-1]   # top-down
    temf = torch.from_numpy(T_td.copy()).reshape(1, -1)
    ds.forward(prop, temf=temf,
               btemp=torch.tensor([float(T_levels_ascending[0])]),
               ttemp=torch.tensor([float(T_levels_ascending[-1])]),
               temis=torch.tensor([[1.0]]),
               albedo=torch.tensor([[float(albedo)]]))
    return _net_down(ds)[::-1].copy()
