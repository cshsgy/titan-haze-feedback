"""Full eddy-diffusion boundary-value problem (``docs/scaling_law.md`` Sec. 3).

Solves the coupled steady-state continuity equations for the monomer-volume
density ``M = n*Nbar`` and the aggregate number density ``n``:

    d/dz[ omega(Nbar,z) M + K dM/dz ] = -S_M(z)                  (monomers)
    d/dz[ omega(Nbar,z) n + K dn/dz ] = 1/2 beta(Nbar,z) n^2 - S_n(z)   (number)

cast as a first-order system in the state ``y = [M, Phi_M, n, Phi_n]`` where
``Phi_q = omega q + K dq/dz`` is the downward flux of ``q``:

    dM/dz     = (Phi_M - omega M) / K
    dPhi_M/dz = 0
    dn/dz     = (Phi_n - omega n) / K
    dPhi_n/dz = 1/2 beta n^2

The narrow Gaussian production (peak ``z0``, width ``dz``) is concentrated at the
top of the domain and imposed as a flux boundary condition there, rather than as
an interior source.  This is the same idealization the master ODE makes and it
removes the poorly constrained, near-empty region above the source that
otherwise makes the collocation Jacobian singular.  Boundary conditions:

  * top (``z0``): downward monomer flux ``Phi_M = P`` and seed number flux
    ``Phi_n = P / N_seed`` (each incoming seed carries ``N_seed`` monomers);
  * surface (``z=0``): settling-only deposition, zero diffusive flux
    (``dM/dz = dn/dz = 0`` i.e. ``Phi_q = omega q``).

Because there is no interior source, ``Phi_M`` is constant ``= P`` (exact monomer
conservation), and the surface condition then forces ``omega(0) M(0) = P``: all
produced monomers deposit at the surface.  The K -> 0 master ODE
(``scaling_law.solve_scaling_law``) supplies the initial guess on [0, z0].
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.integrate import solve_bvp

from .atmosphere import Atmosphere
from .constants import AerosolParams, DEFAULT
from . import transport as tr
from .scaling_law import solve_scaling_law

_TINY = 1e-30


def production_M(z, p: AerosolParams):
    """Monomer-volume production rate S_M(z) [monomers / m^3 / s] (Gaussian)."""
    z = np.asarray(z, dtype=float)
    G = np.exp(-((z - p.z0) ** 2) / (2.0 * p.dz**2)) / (np.sqrt(2.0 * np.pi) * p.dz)
    return p.P_flux * G  # integrates over z to P_flux


@dataclass
class BVPResult:
    z: np.ndarray
    Nbar: np.ndarray
    n: np.ndarray
    M: np.ndarray
    rho_h: np.ndarray
    r: np.ndarray
    r_a: np.ndarray
    omega: np.ndarray
    flux_M: np.ndarray      # downward monomer flux Phi_M [monomers/m^2/s]
    atm: Atmosphere
    params: AerosolParams
    success: bool
    message: str


def _nbar(M, n):
    return np.maximum(np.maximum(M, _TINY) / np.maximum(n, _TINY), 1e-12)


def solve_bvp_profile(atm: Atmosphere | None = None,
                      params: AerosolParams = DEFAULT,
                      n_nodes: int = 400,
                      z_top: float | None = None,
                      tol: float = 1e-3,
                      max_nodes: int = 200000,
                      verbose: int = 0) -> BVPResult:
    """Solve the eddy-diffusion BVP with the *distributed* Gaussian production.

    Production is an interior source over the domain [0, z_top] (z_top a few
    Gaussian widths above z0), not a delta at the boundary -- this removes the
    sharp top-boundary number spike of the flux-BC formulation.  Boundary
    conditions: a closed top (zero flux, no haze enters from above) and
    settling-only deposition at the surface; together these still enforce global
    mass balance (surface deposition == column production).
    """
    if atm is None:
        atm = Atmosphere.titan_reference()
    p = params
    P = p.P_flux
    if z_top is None:
        z_top = p.z0 + 3.0 * p.dz          # ~3 Gaussian widths above the peak

    z_grid = np.linspace(0.0, z_top, n_nodes)

    # --- initial guess from the K->0 master ODE on [0, z0], extended above ---
    master = solve_scaling_law(atm, p, n_out=max(n_nodes, 400), z_bottom=0.0)
    zm = master.z[::-1]                  # ascending
    Nbar0 = np.where(z_grid <= p.z0,
                     np.interp(z_grid, zm, master.Nbar[::-1]), p.N_seed)
    n_z0 = float(np.interp(p.z0, zm, master.n[::-1]))
    n0 = np.where(z_grid <= p.z0,
                  np.interp(z_grid, zm, master.n[::-1]),
                  n_z0 * np.exp(-(z_grid - p.z0) / p.dz))
    M0 = n0 * Nbar0
    omega0 = tr.settling_velocity(Nbar0, z_grid, atm, p)
    y0 = np.vstack([M0, omega0 * M0, n0, omega0 * n0])

    # --- residual (interior Gaussian source) and BCs ---
    def fun(z, y):
        M, PhiM, n, Phin = y
        Nbar = _nbar(M, n)
        omega = tr.settling_velocity(Nbar, z, atm, p)
        beta = tr.coag_kernel(Nbar, z, atm, p)
        K = atm.K(z)
        SM = production_M(z, p)             # monomer-volume source [1/m^3/s]
        Sn = SM / p.N_seed                  # seed number source
        return np.vstack([
            (PhiM - omega * M) / K,
            -SM,
            (Phin - omega * n) / K,
            0.5 * beta * n * n - Sn,
        ])

    def bc(ya, yb):
        Nbar_s = _nbar(ya[0], ya[2])
        omega_s = tr.settling_velocity(Nbar_s, 0.0, atm, p)
        return np.array([
            ya[1] - omega_s * ya[0],   # surface: zero diffusive flux (M)
            ya[3] - omega_s * ya[2],   # surface: zero diffusive flux (n)
            yb[1],                     # top: closed (zero monomer flux)
            yb[3],                     # top: closed (zero number flux)
        ])

    sol = solve_bvp(fun, bc, z_grid, y0, tol=tol, max_nodes=max_nodes,
                    verbose=verbose)

    z = sol.x
    M, PhiM, n, Phin = sol.y
    n = np.maximum(n, _TINY)
    M = np.maximum(M, _TINY)
    Nbar = _nbar(M, n)
    omega = tr.settling_velocity(Nbar, z, atm, p)
    rho_h = p.m_mono * M
    r = tr.mass_radius(Nbar, p.d_mono)
    r_a = tr.mobility_radius(Nbar, p.d_mono, p.D_f)

    return BVPResult(z=z, Nbar=Nbar, n=n, M=M, rho_h=rho_h, r=r, r_a=r_a,
                     omega=omega, flux_M=PhiM, atm=atm, params=p,
                     success=sol.success, message=sol.message)
