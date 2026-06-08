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
from rt.energy_balance import compute_fluxes, radiative_equilibrium, SolarForcing, S0_TITAN


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
    fig, ax = plt.subplots(1, 3, figsize=(14, 5.5))
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

    for a in ax:
        a.grid(alpha=0.3); a.set_ylim(0, 450)
    fig.suptitle("Titan 1-D DISORT energy balance (two-band, haze from Step 2 microphysics)")
    fig.tight_layout()
    out = ROOT / "writing" / "figs" / "energy_balance.png"
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=130)
    print(f"\nwrote {out.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
