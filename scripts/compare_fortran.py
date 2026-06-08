#!/usr/bin/env python3
"""Overlay the reference Fortran model (src/example_bowen_fort) on our DISORT
energy balance.

Runs our Step 1 + Step 2 model, then, if the Fortran model has been built and
run (producing temperatures.txt / sw.txt / lw.txt in its run directory),
overlays its temperature and heating-rate profiles for comparison.

    .rtenv/bin/python scripts/compare_fortran.py [FORTRAN_RUN_DIR]

If the Fortran outputs are not found (the model needs its INPUT/DATA tables,
namelist, and the missing haze.F90 / read_clim modules to run), the comparison
panels are drawn with our model only and a note is printed.
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
from rt.energy_balance import radiative_equilibrium, compute_fluxes


TITAN_DAY = 1.378e6   # s


def load_fortran(run_dir: Path):
    """Load Fortran outputs: temperatures.txt (row0=pressure[Pa], last row=T[K]),
    sw.txt/lw.txt (heating rate [K/s], last row).  Returns dict or None."""
    T_file = run_dir / "temperatures.txt"
    if not T_file.exists():
        return None
    try:
        T = np.loadtxt(T_file)
        if T.ndim != 2 or T.shape[0] < 2:
            return None
        P = T[0]                                   # pressure grid [Pa]
        Tprof = T[-1]                              # last (most evolved) snapshot
        if not np.all(np.isfinite(Tprof)):
            print("warning: Fortran temperature snapshot contains NaN")
            return None
        out = {"P": P, "T": Tprof}
        for key, fname in (("sw", "sw.txt"), ("lw", "lw.txt")):
            f = run_dir / fname
            if f.exists():
                a = np.loadtxt(f)
                out[key] = (a[-1] if a.ndim == 2 else a) * TITAN_DAY  # K/Titan-day
        return out
    except Exception as e:
        print(f"warning: could not parse Fortran outputs: {e}")
        return None


def main():
    run_dir = Path(sys.argv[1]) if len(sys.argv) > 1 else \
        ROOT / "src" / "example_bowen_fort"

    # our model
    atm = Atmosphere.titan_reference()
    micro = solve_bvp_profile(atm, DEFAULT, n_nodes=200)
    col = Column.from_atmosphere(atm, nlyr=40, z_top=430e3)
    eq, _ = radiative_equilibrium(col, micro)
    fx = compute_fluxes(eq, micro)

    fort = load_fortran(run_dir)

    # compare vs PRESSURE (both models share the pressure coordinate)
    fig, ax = plt.subplots(1, 2, figsize=(11, 5.5))
    ax[0].plot(eq.T, eq.P, "r-", lw=1.8, label="this work (DISORT)")
    ax[1].plot(fx.sw_heating, eq.P_mid, "C0", label="SW (this work)")
    ax[1].plot(fx.lw_heating, eq.P_mid, "C1", label="LW (this work)")

    if fort is not None:
        print(f"overlaying Fortran outputs from {run_dir}")
        ax[0].plot(fort["T"], fort["P"], "k--", lw=1.5, label="example_bowen_fort")
        if "sw" in fort:
            ax[1].plot(fort["sw"], fort["P"], "k--", label="SW (Fortran)")
        if "lw" in fort:
            ax[1].plot(fort["lw"], fort["P"], "k:", label="LW (Fortran)")
        title_note = "vs reference Fortran model (TAM-derived)"
    else:
        print("Fortran outputs not found -- plotting this work only.")
        print("Run scripts/run_fortran.sh to enable the overlay.")
        title_note = "(run scripts/run_fortran.sh for the Fortran overlay)"

    for a in ax:
        a.set_yscale("log"); a.set_ylim(1.5e5, 1.0)   # surface at bottom
        a.grid(alpha=0.3)
    ax[0].set_xlabel("temperature [K]"); ax[0].set_ylabel("pressure [Pa]")
    ax[0].set_title("temperature"); ax[0].legend()
    ax[1].axvline(0, color="grey", lw=0.6); ax[1].set_xlim(-80, 80)
    ax[1].set_xlabel("rate [K / Titan-day]"); ax[1].set_title("heating / cooling")
    ax[1].legend()
    fig.suptitle(f"DISORT energy balance {title_note}")
    fig.tight_layout()
    out = ROOT / "writing" / "figs" / "fortran_comparison.png"
    fig.savefig(out, dpi=130)
    print(f"wrote {out.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
