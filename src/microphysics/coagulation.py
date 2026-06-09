"""Bimodal coagulation moment tendencies by log-normal quadrature (piece 1b).

The Burgalat & Rannou (2017) two-mode 2-moment scheme writes the coagulation
tendencies of (M0, M3) for a spherical mode S and a fractal mode F.  Rather than
hand-assemble their closed-form ``gamma`` coefficients (which also rely on a
power-law approximation of the free-molecular kernel with numerically-tabulated
``b_k`` factors), we evaluate the moment integrals (B&R Eqs. 6-7) directly by
Gauss-Hermite quadrature over each mode's log-normal, with the FULL Fuchs kernel
(continuum * free-molecular / sum).  Same physics, exact in the free-molecular
limit, conserves monomer volume by construction.

``coag_tendencies`` is vectorized over a leading "node" axis (so a whole BVP grid
is evaluated in one call); scalar inputs return scalars.  It returns the four
coagulation tendencies (dM0S, dM3S, dM0F, dM3F) [per second], including the S->F
inter-mode transfer (an S+S product larger than the monomer goes to F; B&R Eq. 26).
"""

from __future__ import annotations

import numpy as np

from . import constants as C

_PI = np.pi
_FLOOR = 1e-300


def _nodes(M0, M3, sigma, xi, gw):
    """Log-normal GH nodes: bulk radii r[...,nq] and number weights w[...,nq]."""
    s = np.log(sigma)
    r0 = np.cbrt(np.maximum(M3, _FLOOR) / np.maximum(M0, _FLOOR))      # (...,)
    r_g = r0 * np.exp(-1.5 * s * s)
    r = r_g[..., None] * np.exp(np.sqrt(2.0) * s * xi)                 # (...,nq)
    w = (gw / np.sqrt(_PI)) * M0[..., None]
    return r, w


def _apparent(r, Df, r_m):
    return r_m ** ((Df - 3.0) / Df) * r ** (3.0 / Df)


def _beta(ri, rai, rj, raj, T, eta, lam, rho_p, A):
    """Fuchs kernel on a (...,nq,nq) grid (ri,rai end in (...,nq,1); rj,raj in
    (...,1,nq); T,eta,lam in (...,1,1))."""
    Di = C.K_B * T / (6.0 * _PI * eta * rai) * (1.0 + A * lam / rai)
    Dj = C.K_B * T / (6.0 * _PI * eta * raj) * (1.0 + A * lam / raj)
    beta_co = 4.0 * _PI * (rai + raj) * (Di + Dj)
    beta_fm = np.sqrt(6.0 * C.K_B * T / rho_p) * (rai + raj) ** 2 * np.sqrt(ri ** -3 + rj ** -3)
    return beta_co * beta_fm / (beta_co + beta_fm)


def coag_tendencies(M0S, M3S, M0F, M3F, z, atm, p, sigma_s, sigma_f, n=20):
    """(dM0S, dM3S, dM0F, dM3F) [1/s, 1/s].  Vectorized over the shape of z."""
    scalar = np.ndim(z) == 0
    z = np.atleast_1d(np.asarray(z, float))
    M0S, M3S, M0F, M3F = (np.broadcast_to(np.atleast_1d(np.asarray(m, float)), z.shape).copy()
                          for m in (M0S, M3S, M0F, M3F))
    T = np.atleast_1d(atm.temperature(z))[..., None, None]
    eta = np.atleast_1d(atm.viscosity(z))[..., None, None]
    lam = np.atleast_1d(atm.mfp(z))[..., None, None]
    rho, A, r_m, Df = p.rho_p, p.A_slip, p.d_mono, p.D_f
    xi, gw = np.polynomial.hermite.hermgauss(n)

    dM0S = np.zeros_like(z); dM3S = np.zeros_like(z)
    dM0F = np.zeros_like(z); dM3F = np.zeros_like(z)

    def pair(rA, wA, DfA, rB, wB, DfB):
        raA = _apparent(rA, DfA, r_m); raB = _apparent(rB, DfB, r_m)
        ri = rA[..., :, None]; rj = rB[..., None, :]
        b = _beta(ri, raA[..., :, None], rj, raB[..., None, :], T, eta, lam, rho, A)
        ww = wA[..., :, None] * wB[..., None, :]
        vol = ri ** 3 + rj ** 3
        return ri, rj, b, ww, vol

    # ---- S + S (spherical, Df=3) ----
    rS, wS = _nodes(M0S, M3S, sigma_s, xi, gw)
    ri, rj, b, ww, vol = pair(rS, wS, 3.0, rS, wS, 3.0)
    stay = (vol ** (1.0 / 3.0) <= r_m)
    for k in (0, 3):
        g = vol ** (k / 3.0)
        gain = 0.5 * np.sum(ww * b * g, axis=(-2, -1))
        gain_stay = 0.5 * np.sum(ww * b * g * stay, axis=(-2, -1))
        loss = np.sum(ww * b * ri ** k, axis=(-2, -1))
        if k == 0:
            dM0S += gain_stay - loss; dM0F += gain - gain_stay
        else:
            dM3S += gain_stay - loss; dM3F += gain - gain_stay

    # ---- F + F (fractal) ----
    rF, wF = _nodes(M0F, M3F, sigma_f, xi, gw)
    ri, rj, b, ww, vol = pair(rF, wF, Df, rF, wF, Df)
    for k in (0, 3):
        net = (0.5 * np.sum(ww * b * vol ** (k / 3.0), axis=(-2, -1))
               - np.sum(ww * b * ri ** k, axis=(-2, -1)))
        if k == 0:
            dM0F += net
        else:
            dM3F += net

    # ---- S + F (i=S, j=F) -> F ----
    ri, rj, b, ww, vol = pair(rS, wS, 3.0, rF, wF, Df)
    for k in (0, 3):
        lossS = np.sum(ww * b * ri ** k, axis=(-2, -1))
        gainF = np.sum(ww * b * (vol ** (k / 3.0) - rj ** k), axis=(-2, -1))
        if k == 0:
            dM0S += -lossS; dM0F += gainF
        else:
            dM3S += -lossS; dM3F += gainF

    out = (dM0S, dM3S, dM0F, dM3F)
    return tuple(np.asarray(o).reshape(-1)[0].item() for o in out) if scalar else out
