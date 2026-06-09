"""K->0 bimodal master ODE (polydisperse piece 2).

Sedimentation-dominated limit of the bimodal 2-moment scheme: with eddy
diffusion off, the downward flux of each moment, ``Phi_k = <w>_k M_k``, evolves
only through coagulation,

    dPhi_k/dz = - (dM_k/dt)_coag ,

integrated downward from the production altitude z0.  At each level the four
moments (M0,M3 for the spherical S and fractal F modes) are recovered from the
four fluxes per mode by a 1-D root find on the mean radius r0 (the moment-averaged
settling velocity depends on the moments only through r0).  Production injects
seeds of radius r_C=r_seed into the S mode at z0; the S->F transfer is handled
inside the coagulation tendencies.  Serves as a standalone profile and as the
initial guess for the full eddy-diffusion BVP (piece 3).
"""

from __future__ import annotations

import contextlib
import os
from dataclasses import dataclass

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq


@contextlib.contextmanager
def _quiet_fortran():
    """Silence the underlying Fortran solver's benign fd-2 warning flood."""
    fd = os.dup(2)
    devnull = os.open(os.devnull, os.O_WRONLY)
    os.dup2(devnull, 2)
    try:
        yield
    finally:
        os.dup2(fd, 2)
        os.close(devnull); os.close(fd)

from .atmosphere import Atmosphere
from .constants import AerosolParams, DEFAULT
from . import moments as mm
from . import coagulation as cg

_FLOOR = 1e-30


@dataclass
class BimodalResult:
    z: np.ndarray
    M0S: np.ndarray; M3S: np.ndarray      # spherical mode moments
    M0F: np.ndarray; M3F: np.ndarray      # fractal mode moments
    r0S: np.ndarray; r0F: np.ndarray      # mean radii [m]
    rho_h: np.ndarray                     # total haze mass density [kg/m^3]
    atm: Atmosphere; params: AerosolParams
    sigma_s: float; sigma_f: float


def to_micro(res):
    """Adapt a bimodal result to the (z, n, Nbar, r_a, params) interface the
    optics expect.  The RDG haze opacity depends only on the monomer-volume column
    ``n*Nbar`` (= total monomer count density ``(M3S+M3F)/d^3``), so the bimodal
    haze plugs straight into ``rt.optics`` / the coupling.  ``Nbar`` and ``r_a``
    are the mass-weighted means (the optically dominant F mode)."""
    from types import SimpleNamespace
    d = res.params.d_mono
    mono = (res.M3S + res.M3F) / d ** 3                 # monomer count density [1/m^3]
    n_tot = np.maximum(res.M0S + res.M0F, _FLOOR)
    Nbar = mono / n_tot
    r = d * np.cbrt(np.maximum(Nbar, 1e-12))
    r_a = d * np.maximum(Nbar, 1e-12) ** (1.0 / res.params.D_f)
    return SimpleNamespace(z=res.z, n=n_tot, Nbar=Nbar, r=r, r_a=r_a,
                           rho_h=res.rho_h, params=res.params)


def _recover(Phi0, Phi3, sigma, Df, z, atm, p, lo=1e-9, hi=2e-4):
    """Moments (M0, M3) from the fluxes (Phi0, Phi3): solve (<w>3/<w>0) r0^3 =
    Phi3/Phi0 for r0, then M0 = Phi0/<w>0.  Settling depends on (M0,M3) only via
    r0, so the unit-moment probe (1, r0^3) evaluates <w>_k at that r0."""
    if Phi0 <= _FLOOR or Phi3 <= _FLOOR:
        return _FLOOR, _FLOOR
    target = Phi3 / Phi0

    def g(lr0):
        r0 = np.exp(lr0)
        w0, w3 = mm.settling_velocities(1.0, r0 ** 3, sigma, Df, z, atm, p)
        return np.log((w3 / w0) * r0 ** 3) - np.log(target)

    a, b = np.log(lo), np.log(hi)
    if g(a) * g(b) > 0:                    # outside bracket: clamp
        r0 = lo if abs(g(a)) < abs(g(b)) else hi
    else:
        r0 = np.exp(brentq(g, a, b, xtol=1e-6))
    w0, _ = mm.settling_velocities(1.0, r0 ** 3, sigma, Df, z, atm, p)
    M0 = Phi0 / w0
    return M0, M0 * r0 ** 3


def solve_bimodal_kzero(atm: Atmosphere | None = None,
                        params: AerosolParams = DEFAULT,
                        sigma_s: float = 1.5, sigma_f: float = 2.0,
                        n_out: int = 200, n_quad: int = 16) -> BimodalResult:
    atm = atm or Atmosphere.titan_reference()
    p = params
    z0 = p.z0
    # production fluxes at z0 (column-integrated): M3 volume-fraction flux and
    # the seed number flux (seeds of radius r_seed).
    P3 = (3.0 / (4.0 * np.pi)) * p.Q_mass / p.rho_p          # [m/s]
    P0 = P3 / p.r_seed ** 3                                  # [1/m^2/s]
    Phi0 = np.array([P3 / p.r_seed ** 3, _FLOOR])            # [S, F] number flux
    Phi3 = np.array([P3, _FLOOR])                            # [S, F] volume flux
    y0 = np.array([Phi0[0], Phi3[0], Phi0[1], Phi3[1]])

    def rhs(z, y):
        P0S, P3S, P0F, P3F = np.maximum(y, _FLOOR)
        M0S, M3S = _recover(P0S, P3S, sigma_s, 3.0, z, atm, p)
        M0F, M3F = _recover(P0F, P3F, sigma_f, p.D_f, z, atm, p)
        d = cg.coag_tendencies(M0S, M3S, M0F, M3F, z, atm, p,
                               sigma_s, sigma_f, n=n_quad)
        return [-d[0], -d[1], -d[2], -d[3]]                 # dPhi/dz = -coag

    z_eval = np.linspace(z0, 0.0, n_out)
    # LSODA is robust here; its benign "t+h=t" warnings at the empty-F start
    # (Fortran fd-2) are suppressed.
    with _quiet_fortran():
        sol = solve_ivp(rhs, (z0, 0.0), y0, t_eval=z_eval, method="LSODA",
                        rtol=1e-4, atol=1e-30)

    M0S = np.empty(n_out); M3S = np.empty(n_out)
    M0F = np.empty(n_out); M3F = np.empty(n_out)
    for i, z in enumerate(sol.t):
        M0S[i], M3S[i] = _recover(*np.maximum(sol.y[0:2, i], _FLOOR), sigma_s, 3.0, z, atm, p)
        M0F[i], M3F[i] = _recover(*np.maximum(sol.y[2:4, i], _FLOOR), sigma_f, p.D_f, z, atm, p)
    rho_h = p.rho_p * (4.0 * np.pi / 3.0) * (M3S + M3F)
    return BimodalResult(z=sol.t, M0S=M0S, M3S=M3S, M0F=M0F, M3F=M3F,
                         r0S=mm.mean_radius(M0S, M3S), r0F=mm.mean_radius(M0F, M3F),
                         rho_h=rho_h, atm=atm, params=p,
                         sigma_s=sigma_s, sigma_f=sigma_f)
