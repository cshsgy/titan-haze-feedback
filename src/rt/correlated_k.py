"""Correlated-k shortwave methane opacity.

Reuses the same CH4 visible k-coefficient table as the reference Fortran model
(``example_bowen_fort/INPUT/DATA/ckc_CH4vis.txt``).  Replaces the gray shortwave
gas placeholder: per spectral band the k-distribution is sampled at Gauss points,
and the layer optical depth is ``tau = k(P,T) * col_abund`` with the column
abundance in km-amagat (matching the Fortran ``get_taukcoeff``).

File layout (per ``setspv``/read loop in the Fortran driver): header; band
centres [cm^-1]; Gauss weights; then for each of ``kP`` pressures a pressure
line, and for each of ``kT`` temperatures a temperature line followed by
``nband`` lines of ``ngauss`` k values.  Grid: P 0.1-1e6 Pa, T 40-400 K.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

R_GAS = 8.314462618
MOLE_TO_AMG = 0.02241397            # amg m^3 mol^-1 (Fortran mole_to_amg)

_ROOT = Path(__file__).resolve().parents[2]
_CKC = _ROOT / "src" / "example_bowen_fort" / "INPUT" / "DATA" / "ckc_CH4vis.txt"
_SOLAR = _ROOT / "src" / "example_bowen_fort" / "INPUT" / "DATA" / "solar_spectrum_houghton.txt"

# shortwave band edges (bwnv): 1e6 / wavelength[nm], identical to the Fortran
_WL_NM = np.array([600, 550, 500, 450, 400, 350, 300,
                   250, 240, 230, 220, 210, 200, 190, 180, 170, 160,
                   150, 145, 140, 135, 130, 125, 120, 115, 110, 105, 100,
                   95, 90, 85, 80, 75, 70, 65, 60, 55, 50, 45, 40, 35, 30, 25],
                  dtype=float)
BWNV = 1.0e6 / _WL_NM               # 43 ascending band edges


def parse_ckc(path=_CKC, kP=20, kT=10):
    lines = Path(path).read_text().split("\n")
    bands = np.array(lines[1].split(), float)
    gw = np.array(lines[2].split(), float)
    nb, ng = bands.size, gw.size
    pres = np.zeros(kP)
    temp = np.zeros(kT)
    k = np.zeros((kP, kT, nb, ng))
    i = 3
    for p in range(kP):
        pres[p] = float(lines[i]); i += 1
        for t in range(kT):
            temp[t] = float(lines[i]); i += 1
            for n in range(nb):
                k[p, t, n, :] = np.array(lines[i].split(), float); i += 1
    return dict(bands=bands, gw=gw, pres=pres, temp=temp, k=k)


def solar_per_band(path=_SOLAR, bwnv=BWNV):
    """Band-integrated TOA solar flux at 1 AU [W/m^2] (integral of solarf dnu)."""
    lines = Path(path).read_text().split("\n")
    n = int(lines[1])
    data = np.array([lines[2 + m].split() for m in range(n)], float)
    wn, sf = data[:, 0], data[:, 1]
    edges = np.empty(n + 1)
    edges[1:-1] = 0.5 * (wn[1:] + wn[:-1])
    edges[0] = wn[0] - 0.5 * (wn[1] - wn[0])
    edges[-1] = wn[-1] + 0.5 * (wn[-1] - wn[-2])
    out = np.zeros(bwnv.size - 1)
    for m in range(bwnv.size - 1):
        lo, hi = bwnv[m], bwnv[m + 1]
        ov = np.clip(np.minimum(edges[1:], hi) - np.maximum(edges[:-1], lo), 0.0, None)
        out[m] = np.sum(sf * ov)
    return out


class CorrelatedKSW:
    """Correlated-k CH4 shortwave optical depths + per-band solar flux."""

    def __init__(self, sma_au: float = 9.58, ckc=_CKC, solar=_SOLAR):
        d = parse_ckc(ckc)
        self.bands = d["bands"]
        self.gw = d["gw"]
        self.pres = d["pres"]
        self.temp = d["temp"]
        self.k = d["k"]                       # (kP, kT, nband, ngauss)
        self.logp = np.log(self.pres)
        self.nband = self.bands.size
        self.ngauss = self.gw.size
        # per-band TOA solar flux scaled to Titan's orbit
        self.solar = solar_per_band(solar, BWNV) / sma_au**2

    def _interp_k(self, P, T):
        """k(P, T) -> (nband, ngauss); log-linear in P, linear in T (cf. Fortran)."""
        lp = np.log(np.clip(P, self.pres[0], self.pres[-1]))
        jp = np.clip(np.searchsorted(self.logp, lp) - 1, 0, self.pres.size - 2)
        fp = (lp - self.logp[jp]) / (self.logp[jp + 1] - self.logp[jp])
        kP = (1 - fp) * self.k[jp] + fp * self.k[jp + 1]      # (kT, nband, ngauss)
        Tc = np.clip(T, self.temp[0], self.temp[-1])
        jt = np.clip(np.searchsorted(self.temp, Tc) - 1, 0, self.temp.size - 2)
        ft = (Tc - self.temp[jt]) / (self.temp[jt + 1] - self.temp[jt])
        return (1 - ft) * kP[jt] + ft * kP[jt + 1]            # (nband, ngauss)

    def layer_tau(self, column, comp=None):
        """CH4 gas optical depth tau[nband, ngauss, nlyr] for the column layers."""
        if comp is None:
            from .cia import Composition
            comp = Composition()
        T, P = column.T_mid, column.P_mid
        # CH4 column abundance per layer [km-amagat]
        col_abund = (comp.x_CH4 * P / (R_GAS * T) * MOLE_TO_AMG
                     * (column.dz / 1e3))
        nlyr = T.size
        tau = np.zeros((self.nband, self.ngauss, nlyr))
        for l in range(nlyr):
            tau[:, :, l] = self._interp_k(P[l], T[l]) * col_abund[l]
        return np.maximum(tau, 0.0)
