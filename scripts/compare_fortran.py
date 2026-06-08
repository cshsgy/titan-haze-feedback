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


def load_fortran(run_dir: Path):
    """Load Fortran outputs if present: returns dict or None."""
    T_file = run_dir / "temperatures.txt"
    if not T_file.exists():
        return None
    out = {}
    try:
        out["T"] = np.loadtxt(T_file)
        for key, fname in (("sw", "sw.txt"), ("lw", "lw.txt")):
            f = run_dir / fname
            if f.exists():
                out[key] = np.loadtxt(f)
    except Exception as e:
        print(f"warning: could not parse Fortran outputs: {e}")
        return None
    return out


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

    fig, ax = plt.subplots(1, 2, figsize=(11, 5.5))
    ax[0].plot(eq.T, eq.z / 1e3, "r-", label="this work (DISORT)")
    ax[1].plot(fx.sw_heating, fx.z_mid / 1e3, "C0", label="SW (this work)")
    ax[1].plot(fx.lw_heating, fx.z_mid / 1e3, "C1", label="LW (this work)")

    if fort is not None:
        # Fortran grid is pressure-based (nlay layers); plot vs an index proxy if
        # no altitude file -- here we simply overlay on the same axis order.
        print(f"overlaying Fortran outputs from {run_dir}")
        n = fort["T"].size
        yidx = np.linspace(eq.z.min(), eq.z.max(), n) / 1e3
        ax[0].plot(fort["T"], yidx, "k--", label="example_bowen_fort")
        if "sw" in fort:
            ax[1].plot(fort["sw"], np.linspace(0, 430, fort["sw"].size),
                       "k--", label="SW (Fortran)")
        if "lw" in fort:
            ax[1].plot(-np.abs(fort["lw"]), np.linspace(0, 430, fort["lw"].size),
                       "k:", label="LW (Fortran)")
        title_note = "with reference Fortran model"
    else:
        print("Fortran outputs not found -- plotting this work only.")
        print("Build+run src/example_bowen_fort (needs INPUT/DATA, namelist, "
              "and the missing haze.F90 / read_clim modules) to enable overlay.")
        title_note = "(Fortran model not yet run -- see src/example_bowen_fort/README)"

    ax[0].set_xlabel("temperature [K]"); ax[0].set_ylabel("altitude [km]")
    ax[0].set_title("temperature"); ax[0].legend(); ax[0].grid(alpha=0.3)
    ax[1].axvline(0, color="grey", lw=0.6); ax[1].set_xlim(-60, 60)
    ax[1].set_xlabel("rate [K / Titan-day]"); ax[1].set_title("heating / cooling")
    ax[1].legend(); ax[1].grid(alpha=0.3)
    for a in ax:
        a.set_ylim(0, 430)
    fig.suptitle(f"DISORT energy balance vs reference Fortran {title_note}")
    fig.tight_layout()
    out = ROOT / "writing" / "figs" / "fortran_comparison.png"
    fig.savefig(out, dpi=130)
    print(f"wrote {out.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
