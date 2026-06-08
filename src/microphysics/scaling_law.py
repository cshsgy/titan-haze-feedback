"""The microphysical scaling law (Step 2).

Sedimentation-dominated master ODE (``docs/scaling_law.md`` Eq. 5):

    dNbar/dz = - beta(Nbar, z) * P / (2 omega(Nbar, z)^2)

integrated downward from the production altitude z0 (where Nbar = N_seed) to the
surface, with the two-phase fractal dimension switching at Nbar = 1.  From the
mean size the number density and mass density follow from monomer-flux
conservation (Eq. 4 of the doc):

    n(z)      = P / (omega Nbar)
    rho_h(z)  = Q_p / omega

This is the K -> 0 limit of the full eddy-diffusion BVP and serves both as a
standalone law and as the initial guess for the BVP solver (future work).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.integrate import solve_ivp

from .atmosphere import Atmosphere
from .constants import AerosolParams, DEFAULT
from . import transport as tr


@dataclass
class ScalingResult:
    z: np.ndarray          # altitude grid [m]
    Nbar: np.ndarray       # mean monomer-volume count
    n: np.ndarray          # aggregate number density [1/m^3]
    rho_h: np.ndarray      # haze mass density [kg/m^3]
    r: np.ndarray          # mass radius [m]
    r_a: np.ndarray        # mobility radius [m]
    omega: np.ndarray      # settling velocity [m/s]
    atm: Atmosphere
    params: AerosolParams


def _dlnN_dz(z, lnN, atm, p):
    """RHS for ln(Nbar): d(lnN)/dz = -beta P / (2 omega^2 Nbar)."""
    Nbar = np.exp(lnN[0])
    P = p.P_flux
    omega = tr.settling_velocity(Nbar, z, atm, p)
    beta = tr.coag_kernel(Nbar, z, atm, p)
    dNbar_dz = -beta * P / (2.0 * omega**2)
    return [dNbar_dz / Nbar]


def solve_scaling_law(atm: Atmosphere | None = None,
                      params: AerosolParams = DEFAULT,
                      n_out: int = 300,
                      z_bottom: float = 0.0) -> ScalingResult:
    """Integrate the master ODE downward and assemble the haze profile.

    Parameters
    ----------
    atm : Atmosphere, optional
        Background profiles; defaults to the crude Titan reference column.
    params : AerosolParams
        Aerosol/production parameters.
    n_out : int
        Number of output altitude samples between z0 and z_bottom.
    z_bottom : float
        Lowest altitude to integrate to [m].
    """
    if atm is None:
        atm = Atmosphere.titan_reference()
    p = params

    z0 = p.z0
    lnN0 = np.log(p.N_seed)
    z_eval = np.linspace(z0, z_bottom, n_out)

    sol = solve_ivp(
        _dlnN_dz, (z0, z_bottom), [lnN0],
        t_eval=z_eval, args=(atm, p),
        method="LSODA", rtol=1e-8, atol=1e-12, dense_output=True,
    )
    if not sol.success:
        raise RuntimeError(f"scaling-law integration failed: {sol.message}")

    z = sol.t
    Nbar = np.exp(sol.y[0])

    omega = tr.settling_velocity(Nbar, z, atm, p)
    n = p.P_flux / (omega * Nbar)
    rho_h = p.Q_mass / omega
    r = tr.mass_radius(Nbar, p.d_mono)
    r_a = tr.mobility_radius(Nbar, p.d_mono, p.D_f)

    return ScalingResult(z=z, Nbar=Nbar, n=n, rho_h=rho_h, r=r, r_a=r_a,
                         omega=omega, atm=atm, params=p)
