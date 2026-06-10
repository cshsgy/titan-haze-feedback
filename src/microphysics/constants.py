"""Physical constants and baseline Titan parameters.

Values and sources are documented in ``docs/physics_parameters.md``:
T08 = Tomasko et al. (2008), L23 = Lombardo & Lora (2023),
BR17 = Burgalat & Rannou (2017), dT25 = de Trenquelleon et al. (2025).
"""

from __future__ import annotations

from dataclasses import dataclass
import math

# --- universal constants (SI) ---
K_B = 1.380649e-23          # Boltzmann constant [J/K]
N_A = 6.02214076e23         # Avogadro [1/mol]
R_GAS = 8.314462618         # molar gas constant [J/mol/K]
E_CHARGE = 1.602176634e-19  # elementary charge [C]
K_COULOMB = 8.9875517923e9  # Coulomb constant 1/(4 pi eps0) [N m^2/C^2]

# --- Titan body ---
G_SURF = 1.352              # surface gravity [m/s^2]
R_TITAN = 2.575e6           # radius [m]

# --- background gas (N2-dominated) ---
M_N2 = 0.0280134            # molar mass of N2 [kg/mol]
D_GAS = 3.70e-10            # N2 collision diameter [m] (~3.7 angstrom)


@dataclass(frozen=True)
class AerosolParams:
    """Baseline aerosol / production parameters (dT25 Table 1; see physics_parameters.md)."""

    D_f: float = 2.0            # fractal dimension of aggregates (N >= 1 phase)
    d_mono: float = 50e-9       # monomer radius d [m]
    r_seed: float = 10e-9       # production seed radius r_p [m]
    rho_p: float = 800.0        # aerosol material density [kg/m^3]
    Q_mass: float = 2.1e-13     # mass production rate Q_p [kg/m^2/s]
    z0: float = 415e3           # production altitude [m]
    dz: float = 20e3            # production Gaussian width [m]
    n_e: float = 15.0           # charge density [e- / um] (dT25)
    use_charge: bool = False    # apply the Coulomb coagulation inhibition
                                # (off by default: the validated baselines --
                                # H=64 km, column tau~8 -- were calibrated
                                # without it; flip after re-validating)
    A_slip: float = 1.591       # first-order Cunningham slip constant

    @property
    def m_mono(self) -> float:
        """Mass of one monomer [kg]."""
        return self.rho_p * (4.0 / 3.0) * math.pi * self.d_mono**3

    @property
    def N_seed(self) -> float:
        """Seed size in monomer-volume units: (r_p / d)^3  (sub-monomer if r_p < d)."""
        return (self.r_seed / self.d_mono) ** 3

    @property
    def P_flux(self) -> float:
        """Column monomer production flux P = Q_p / m_mono [monomers / m^2 / s]."""
        return self.Q_mass / self.m_mono


DEFAULT = AerosolParams()
