"""Bimodal eddy-diffusion BVP (polydisperse piece 3).

The full steady transport of the bimodal 2-moment state with eddy diffusion on.
For each moment M of mode X carrying downward flux ``Phi = <w>_X M + K dM/dz``:

    dM/dz   = (Phi - <w>_X M) / K
    dPhi/dz = -( (dM/dt)_coag + (dM/dt)_prod )

an 8-component first-order system in ``y = [M0S,Phi0S,M3S,Phi3S, M0F,Phi0F,M3F,Phi3F]``.
Production (a distributed Gaussian) seeds the spherical mode (B&R Eqs. 31-32);
the S->F transfer is inside the coagulation tendencies (piece 1b).  Boundary
conditions: closed top (zero flux, no haze from above) and settling-only
deposition at the surface (zero diffusive flux).  The K->0 bimodal master ODE
(piece 2) supplies the initial guess.

STATUS: the formulation, initial guess, and (vectorized) RHS are correct, but
``scipy.solve_bvp`` (collocation + finite-difference Jacobian) does not converge
in reasonable time for this stiff 8-field system with the quadrature-coagulation
RHS.  The proper solver is a pseudo-time relaxation with operator splitting
(implicit advection-diffusion + explicit/semi-implicit coagulation) -- deferred.
For Titan the eddy-diffusion correction is modest (K->0 gives H~54 km vs the
monodisperse BVP 64 km vs 65 km observed), so downstream work (optics, Step-3
re-coupling) uses the K->0 bimodal profiles (``scaling_law_bimodal``).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.integrate import solve_bvp

from .atmosphere import Atmosphere
from .constants import AerosolParams, DEFAULT
from . import moments as mm
from . import coagulation as cg
from .scaling_law_bimodal import solve_bimodal_kzero

_FLOOR = 1e-30


@dataclass
class BimodalBVPResult:
    z: np.ndarray
    M0S: np.ndarray; M3S: np.ndarray
    M0F: np.ndarray; M3F: np.ndarray
    r0S: np.ndarray; r0F: np.ndarray
    rho_h: np.ndarray
    atm: Atmosphere; params: AerosolParams
    sigma_s: float; sigma_f: float
    success: bool; message: str


def _prod(z, p):
    """Spherical-mode production rates (dM0S, dM3S) [1/m^3/s, 1/s] (Gaussian)."""
    G = np.exp(-((z - p.z0) ** 2) / (2.0 * p.dz ** 2)) / (np.sqrt(2.0 * np.pi) * p.dz)
    dM3 = (3.0 / (4.0 * np.pi)) * (p.Q_mass / p.rho_p) * G
    return dM3 / p.r_seed ** 3, dM3


def solve_bvp_bimodal(atm: Atmosphere | None = None, params: AerosolParams = DEFAULT,
                      sigma_s: float = 1.5, sigma_f: float = 2.0,
                      n_nodes: int = 120, n_quad: int = 10, z_top: float | None = None,
                      tol: float = 1e-3, max_nodes: int = 60000, verbose: int = 0):
    atm = atm or Atmosphere.titan_reference()
    p = params
    z_top = z_top if z_top is not None else p.z0 + 3.0 * p.dz
    z_grid = np.linspace(0.0, z_top, n_nodes)

    # --- initial guess from the K->0 master ODE on [0, z0] ---
    km = solve_bimodal_kzero(atm, p, sigma_s, sigma_f, n_out=max(n_nodes, 200), n_quad=n_quad)
    zk = km.z[::-1]
    def interp_floor(field):
        f = np.interp(z_grid, zk, field[::-1])
        return np.where(z_grid <= p.z0, np.maximum(f, _FLOOR), _FLOOR)
    M0S = interp_floor(km.M0S); M3S = interp_floor(km.M3S)
    M0F = interp_floor(km.M0F); M3F = interp_floor(km.M3F)
    w0S, w3S = mm.settling_velocities(M0S, M3S, sigma_s, 3.0, z_grid, atm, p)
    w0F, w3F = mm.settling_velocities(M0F, M3F, sigma_f, p.D_f, z_grid, atm, p)
    y0 = np.vstack([M0S, w0S * M0S, M3S, w3S * M3S,
                    M0F, w0F * M0F, M3F, w3F * M3F])

    def fun(z, y):
        M0S, F0S, M3S, F3S, M0F, F0F, M3F, F3F = (np.maximum(yi, _FLOOR) for yi in y)
        K = atm.K(z)
        w0S, w3S = mm.settling_velocities(M0S, M3S, sigma_s, 3.0, z, atm, p)
        w0F, w3F = mm.settling_velocities(M0F, M3F, sigma_f, p.D_f, z, atm, p)
        # coagulation per node
        d0S = np.empty_like(z); d3S = np.empty_like(z)
        d0F = np.empty_like(z); d3F = np.empty_like(z)
        for i in range(z.size):
            d0S[i], d3S[i], d0F[i], d3F[i] = cg.coag_tendencies(
                M0S[i], M3S[i], M0F[i], M3F[i], z[i], atm, p, sigma_s, sigma_f, n=n_quad)
        pM0S, pM3S = _prod(z, p)
        return np.vstack([
            (F0S - w0S * M0S) / K, -(d0S + pM0S),
            (F3S - w3S * M3S) / K, -(d3S + pM3S),
            (F0F - w0F * M0F) / K, -d0F,
            (F3F - w3F * M3F) / K, -d3F,
        ])

    def bc(ya, yb):
        a = np.maximum(ya, _FLOOR)
        w0S, w3S = mm.settling_velocities(a[0], a[2], sigma_s, 3.0, 0.0, atm, p)
        w0F, w3F = mm.settling_velocities(a[4], a[6], sigma_f, p.D_f, 0.0, atm, p)
        return np.array([
            ya[1] - w0S * a[0], ya[3] - w3S * a[2],     # surface: settling-only (S)
            ya[5] - w0F * a[4], ya[7] - w3F * a[6],     # surface: settling-only (F)
            yb[1], yb[3], yb[5], yb[7],                 # top: closed (zero flux)
        ])

    sol = solve_bvp(fun, bc, z_grid, y0, tol=tol, max_nodes=max_nodes, verbose=verbose)
    z = sol.x
    M0S, M3S = np.maximum(sol.y[0], _FLOOR), np.maximum(sol.y[2], _FLOOR)
    M0F, M3F = np.maximum(sol.y[4], _FLOOR), np.maximum(sol.y[6], _FLOOR)
    rho_h = p.rho_p * (4.0 * np.pi / 3.0) * (M3S + M3F)
    return BimodalBVPResult(z=z, M0S=M0S, M3S=M3S, M0F=M0F, M3F=M3F,
                            r0S=mm.mean_radius(M0S, M3S), r0F=mm.mean_radius(M0F, M3F),
                            rho_h=rho_h, atm=atm, params=p, sigma_s=sigma_s,
                            sigma_f=sigma_f, success=sol.success, message=sol.message)
