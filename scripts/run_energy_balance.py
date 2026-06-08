#!/usr/bin/env python3
"""Demo: the Step 1 DISORT energy balance on a Titan column with the Step 2 haze.

Computes shortwave heating and longwave cooling, relaxes to radiative
equilibrium, and plots the heating/cooling rates, flux profiles, and the
initial vs. equilibrium temperature.

Run with the RT venv (which has pydisort):
    .rtenv/bin/python scripts/run_energy_balance.py
Writes writing/figs/energy_balance.png
"""

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from microphysics import Atmosphere, DEFAULT, solve_bvp_profile
from rt.column import Column
from rt.cia import CIABands
from rt.energy_balance import compute_fluxes, radiative_equilibrium, S0_TITAN


def main():
    atm = Atmosphere.titan_reference()
    micro = solve_bvp_profile(atm, DEFAULT, n_nodes=200)
    col = Column.from_atmosphere(atm, nlyr=40, z_top=450e3)

    fx0 = compute_fluxes(col, micro)
    print(f"Solar constant at Titan: {S0_TITAN:.2f} W/m^2")
    inc = fx0.sw_net[-1]
    print(f"SW: TOA net {fx0.sw_net[-1]:.3f}, surface {fx0.sw_net[0]:.3f} W/m^2 "
          f"(absorbed {inc - fx0.sw_net[0]:.3f}, {100*(inc-fx0.sw_net[0])/inc:.0f}% of TOA net)")
    print(f"LW: OLR (TOA net up) {-fx0.lw_net[-1]:.3f} W/m^2, "
          f"surface net down {fx0.lw_net[0]:.3f} W/m^2")

    print("\nRelaxing to radiative equilibrium ...")
    eq, hist = radiative_equilibrium(col, micro, verbose=True)
    print(f"residual {hist[0]:.0f} -> {hist[-50:].mean():.1f} K/Titan-day (bulk)")

    fxe = compute_fluxes(eq, micro)

    # --- figure ---
    fig, ax = plt.subplots(1, 4, figsize=(18, 5.5))
    zk = fx0.z_mid / 1e3

    ax[0].plot(fxe.sw_heating, zk, label="SW heating")
    ax[0].plot(fxe.lw_heating, zk, label="LW cooling")
    ax[0].plot(fxe.net_heating, zk, "k--", lw=1.2, label="net")
    ax[0].axvline(0, color="grey", lw=0.6)
    ax[0].set_xlim(-60, 60)
    ax[0].set_xlabel("rate [K / Titan-day]"); ax[0].set_ylabel("altitude [km]")
    ax[0].set_title("heating / cooling (at equilibrium)"); ax[0].legend()

    ax[1].plot(fxe.sw_net, eq.z / 1e3, label="SW net down")
    ax[1].plot(fxe.lw_net, eq.z / 1e3, label="LW net down")
    ax[1].axvline(0, color="grey", lw=0.6)
    ax[1].set_xlabel("net downward flux [W/m$^2$]")
    ax[1].set_title("radiative fluxes"); ax[1].legend()

    ax[2].plot(col.T, col.z / 1e3, label="initial guess")
    ax[2].plot(eq.T, eq.z / 1e3, "r-", label="radiative eq.")
    ax[2].set_xlabel("temperature [K]")
    ax[2].set_title("temperature profile"); ax[2].legend()

    # CIA column optical depth per pair vs wavenumber
    cia = CIABands()
    bc = 0.5 * (cia.band_lo + cia.band_hi)
    from microphysics.constants import K_B
    from rt.cia import Composition
    comp = Composition()
    n_tot = (eq.P_mid / (K_B * eq.T_mid)) * 1e-6
    dz_cm = eq.dz * 100.0
    dens = {"N2-N2": (comp.x_N2 * n_tot) ** 2,
            "N2-CH4": comp.x_N2 * n_tot * comp.x_CH4 * n_tot,
            "CH4-CH4": (comp.x_CH4 * n_tot) ** 2,
            "N2-H2": comp.x_N2 * n_tot * comp.x_H2 * n_tot}
    for pair in ("N2-N2", "N2-CH4", "N2-H2", "CH4-CH4"):
        k = cia._k_at_T(pair, eq.T_mid)                # (nband, nlyr)
        coltau = (k * dens[pair][None, :] * dz_cm[None, :]).sum(axis=1)
        ax[3].plot(bc, np.maximum(coltau, 1e-6), label=pair)
    ax[3].set_yscale("log"); ax[3].set_ylim(1e-3, 1e3)
    ax[3].set_xlabel("wavenumber [cm$^{-1}$]")
    ax[3].set_ylabel("column CIA optical depth")
    ax[3].set_title("collision-induced absorption"); ax[3].legend()
    ax[3].grid(alpha=0.3)

    for a in ax[:3]:
        a.grid(alpha=0.3); a.set_ylim(0, 450)
    fig.suptitle("Titan 1-D DISORT energy balance (shortwave gray gas + multiband CIA longwave; "
                 "haze from Step 2 microphysics)")
    fig.tight_layout()
    out = ROOT / "writing" / "figs" / "energy_balance.png"
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=130)
    print(f"\nwrote {out.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
