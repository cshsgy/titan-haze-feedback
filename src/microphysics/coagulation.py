"""Bimodal coagulation moment tendencies by log-normal quadrature (piece 1b).

The Burgalat & Rannou (2017) two-mode 2-moment scheme writes the coagulation
tendencies of (M0, M3) for a spherical mode S and a fractal mode F.  Rather than
hand-assemble their closed-form ``gamma`` coefficients (which also rely on a
power-law approximation of the free-molecular kernel with numerically-tabulated
``b_k`` factors), we evaluate the moment integrals (B&R Eqs. 6-7) directly by
Gauss-Hermite quadrature over each mode's log-normal, with the FULL Fuchs kernel
(continuum * free-molecular / sum).  This is the same physics, exact in the
free-molecular limit, and conserves monomer volume by construction.  Affordable
for a 1-D column.

State per layer: (M0S, M3S, M0F, M3F).  Returns the coagulation tendencies
(dM0S, dM3S, dM0F, dM3F) [per second], including the S->F inter-mode transfer
(an S+S product larger than the monomer is delivered to F; B&R Eq. 26).
"""

from __future__ import annotations

import numpy as np

from . import constants as C

_PI = np.pi


def _ln_nodes(M0, M3, sigma, n):
    """Gauss-Hermite nodes for a log-normal mode: bulk radii r[n] and number
    weights w[n] (so that integral g(r) n(r) dr ~= sum w_i g(r_i), sum w = M0)."""
    xi, gw = np.polynomial.hermite.hermgauss(n)        # physicists' weight e^{-x^2}
    s = np.log(sigma)
    r0 = np.cbrt(M3 / max(M0, 1e-300))                 # volume-mean radius
    r_g = r0 * np.exp(-1.5 * s * s)                    # number-median radius
    r = r_g * np.exp(np.sqrt(2.0) * s * xi)
    w = (gw / np.sqrt(_PI)) * M0
    return r, w


def _apparent(r, Df, r_m):
    """Mobility/apparent radius r_a = E r^(3/Df), E = r_m^((Df-3)/Df)."""
    E = r_m ** ((Df - 3.0) / Df)
    return E * r ** (3.0 / Df)


def _beta(ri, rai, rj, raj, T, eta, lam, rho_p, A):
    """Fuchs kernel beta(i,j) [m^3/s] on an outer-product grid.  ri,rai are
    column (ni,1); rj,raj are row (1,nj).  Bulk radii r enter the free-molecular
    mass term (m ~ r^3), apparent radii r_a the collision cross-section."""
    Di = C.K_B * T / (6.0 * _PI * eta * rai) * (1.0 + A * lam / rai)
    Dj = C.K_B * T / (6.0 * _PI * eta * raj) * (1.0 + A * lam / raj)
    beta_co = 4.0 * _PI * (rai + raj) * (Di + Dj)
    vth = np.sqrt(6.0 * C.K_B * T / rho_p)
    beta_fm = vth * (rai + raj) ** 2 * np.sqrt(ri ** -3 + rj ** -3)
    return beta_co * beta_fm / (beta_co + beta_fm)


def coag_tendencies(M0S, M3S, M0F, M3F, z, atm, p, sigma_s, sigma_f, n=20):
    """Coagulation tendencies (dM0S, dM3S, dM0F, dM3F) [1/s, m^3/m^3/s]."""
    T = atm.temperature(z); eta = atm.viscosity(z); lam = atm.mfp(z)
    rho, A, r_m = p.rho_p, p.A_slip, p.d_mono
    Df = p.D_f
    dM0S = dM3S = dM0F = dM3F = 0.0

    # ---- S + S (both spherical, Df=3) ----
    if M0S > 0:
        r, w = _ln_nodes(M0S, M3S, sigma_s, n)
        ra = _apparent(r, 3.0, r_m)
        ri = r[:, None]; rj = r[None, :]
        b = _beta(ri, ra[:, None], rj, ra[None, :], T, eta, lam, rho, A)
        ww = w[:, None] * w[None, :]
        vol = ri ** 3 + rj ** 3                        # coalesced volume
        rp = vol ** (1.0 / 3.0)                        # product radius
        stay = (rp <= r_m)                             # stays spherical
        for k in (0, 3):
            gain = 0.5 * np.sum(ww * b * vol ** (k / 3.0))            # all products
            gain_stay = 0.5 * np.sum(ww * b * vol ** (k / 3.0) * stay)
            loss = np.sum(ww * b * ri ** k)            # = (1/2)int beta(ri^k+rj^k)
            if k == 0:
                dM0S += gain_stay - loss; dM0F += gain - gain_stay
            else:
                dM3S += gain_stay - loss; dM3F += gain - gain_stay

    # ---- F + F (fractal, Df) ----
    if M0F > 0:
        r, w = _ln_nodes(M0F, M3F, sigma_f, n)
        ra = _apparent(r, Df, r_m)
        ri = r[:, None]; rj = r[None, :]
        b = _beta(ri, ra[:, None], rj, ra[None, :], T, eta, lam, rho, A)
        ww = w[:, None] * w[None, :]
        vol = ri ** 3 + rj ** 3
        for k in (0, 3):
            net = 0.5 * np.sum(ww * b * vol ** (k / 3.0)) - np.sum(ww * b * ri ** k)
            if k == 0:
                dM0F += net
            else:
                dM3F += net                            # ~0 by volume conservation

    # ---- S + F (i = spherical, j = fractal) -> F ----
    if M0S > 0 and M0F > 0:
        rs, ws = _ln_nodes(M0S, M3S, sigma_s, n)
        rf, wf = _ln_nodes(M0F, M3F, sigma_f, n)
        ras = _apparent(rs, 3.0, r_m); raf = _apparent(rf, Df, r_m)
        ri = rs[:, None]; rj = rf[None, :]
        b = _beta(ri, ras[:, None], rj, raf[None, :], T, eta, lam, rho, A)
        ww = ws[:, None] * wf[None, :]
        vol = ri ** 3 + rj ** 3
        for k in (0, 3):
            lossS = np.sum(ww * b * ri ** k)           # S loses
            gainF = np.sum(ww * b * (vol ** (k / 3.0) - rj ** k))   # F: product - consumed
            if k == 0:
                dM0S += -lossS; dM0F += gainF
            else:
                dM3S += -lossS; dM3F += gainF

    return dM0S, dM3S, dM0F, dM3F
