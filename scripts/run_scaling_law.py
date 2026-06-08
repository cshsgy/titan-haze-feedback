#!/usr/bin/env python3
"""Demo: integrate the microphysics scaling law on the Titan reference column
and plot the resulting haze profiles.

Run from the repo root:
    python3 scripts/run_scaling_law.py
Produces writing/figs/scaling_law_profiles.png
"""

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from microphysics import Atmosphere, DEFAULT, solve_scaling_law


def main():
    atm = Atmosphere.titan_reference()
    res = solve_scaling_law(atm, DEFAULT)
    zkm = res.z / 1e3

    # console summary
    print(f"{'z[km]':>7} {'T[K]':>7} {'Nbar':>12} {'r[um]':>9} "
          f"{'n[1/cm3]':>11} {'omega[m/s]':>11}")
    for zk in (415, 300, 200, 100, 40, 0):
        i = int(np.argmin(np.abs(zkm - zk)))
        print(f"{zkm[i]:7.1f} {atm.temperature(res.z[i]):7.1f} "
              f"{res.Nbar[i]:12.3e} {res.r[i]*1e6:9.4f} "
              f"{res.n[i]*1e-6:11.3e} {res.omega[i]:11.3e}")

    # sanity: monomer mass flux rho_h * omega should equal Q_mass everywhere
    flux = res.rho_h * res.omega
    print(f"\nmass-flux check: rho_h*omega / Q_p in "
          f"[{(flux/DEFAULT.Q_mass).min():.4f}, {(flux/DEFAULT.Q_mass).max():.4f}] "
          f"(should be ~1)")

    # plots
    fig, ax = plt.subplots(1, 4, figsize=(15, 5), sharey=True)
    ax[0].plot(res.Nbar, zkm); ax[0].set_xscale("log")
    ax[0].axvline(1.0, ls=":", c="k", lw=0.8)
    ax[0].set_xlabel(r"$\bar N$ (monomers)"); ax[0].set_ylabel("altitude [km]")
    ax[0].set_title("mean aggregate size")

    ax[1].plot(res.r * 1e6, zkm, label="mass radius $r$")
    ax[1].plot(res.r_a * 1e6, zkm, "--", label="mobility radius $r_a$")
    ax[1].set_xscale("log"); ax[1].set_xlabel(r"radius [$\mu$m]")
    ax[1].set_title("particle radius"); ax[1].legend()

    ax[2].plot(res.n * 1e-6, zkm); ax[2].set_xscale("log")
    ax[2].set_xlabel(r"number density [cm$^{-3}$]")
    ax[2].set_title("aggregate number density")

    ax[3].plot(res.rho_h, zkm); ax[3].set_xscale("log")
    ax[3].set_xlabel(r"mass density [kg m$^{-3}$]")
    ax[3].set_title("haze mass density")

    for a in ax:
        a.grid(alpha=0.3)
    fig.suptitle("Titan haze scaling law (sedimentation-dominated limit)")
    fig.tight_layout()

    outdir = ROOT / "writing" / "figs"
    outdir.mkdir(parents=True, exist_ok=True)
    out = outdir / "scaling_law_profiles.png"
    fig.savefig(out, dpi=130)
    print(f"\nwrote {out.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
