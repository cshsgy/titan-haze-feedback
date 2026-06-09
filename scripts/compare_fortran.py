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
from rt.energy_balance import radiative_equilibrium


TITAN_DAY = 1.378e6   # s


def load_fortran(run_dir: Path, navg: int = 20):
    """Load Fortran outputs: temperatures.txt (row0=pressure[Pa], snapshots[K]),
    sw.txt/lw.txt (heating rate [K/s] per snapshot).

    The reference run does NOT reach a steady state in its upper atmosphere
    (P < ~4 Pa): those layers oscillate by tens of K between snapshots about a
    quasi-equilibrium mean (the explicit time-stepper is unstable for the thin,
    radiatively-fast top layers, despite adaptive timestep reduction).  Below
    ~10 Pa the profile is converged (snapshot-to-snapshot std < 2.5 K).  So we
    return the TIME-AVERAGE over the last ``navg`` finite snapshots plus its std,
    which is a far more representative reference than any single snapshot and
    exposes where the reference itself is unconverged.  Returns dict or None.
    """
    T_file = run_dir / "temperatures.txt"
    if not T_file.exists():
        return None
    try:
        T = np.loadtxt(T_file)
        if T.ndim != 2 or T.shape[0] < 2:
            return None
        P = T[0]                                   # pressure grid [Pa]
        valid = [i for i in range(1, T.shape[0]) if np.all(np.isfinite(T[i]))]
        if not valid:
            print("warning: Fortran temperatures.txt has no finite snapshot")
            return None
        sel = valid[-navg:]                         # last navg finite snapshots
        out = {"P": P, "T": T[sel].mean(0), "T_std": T[sel].std(0),
               "n_avg": len(sel), "n_total": len(valid)}
        for key, fname in (("sw", "sw.txt"), ("lw", "lw.txt")):
            f = run_dir / fname
            if f.exists():
                a = np.loadtxt(f)
                if a.ndim == 2:
                    # sw/lw row j-1 corresponds to temperature snapshot row j
                    rows = [k - 1 for k in sel if 0 <= k - 1 < a.shape[0]]
                    out[key] = a[rows].mean(0) * TITAN_DAY       # K/Titan-day
                else:
                    out[key] = a * TITAN_DAY
        return out
    except Exception as e:
        print(f"warning: could not parse Fortran outputs: {e}")
        return None


def main():
    run_dir = Path(sys.argv[1]) if len(sys.argv) > 1 else \
        ROOT / "src" / "example_bowen_fort"

    from rt.optics import OpticsParams
    # our model -- microphysics haze (the coupled mode) and prescribed haze (RT
    # validation: same haze the Fortran uses)
    atm = Atmosphere.titan_reference()
    micro = solve_bvp_profile(atm, DEFAULT, n_nodes=200)
    col = Column.from_atmosphere(atm, nlyr=100, p_top=1.0)
    # use the fluxes the solver converged with (consistent ck + ck_lw); do NOT
    # recompute here -- that previously dropped ck_lw and mis-plotted the heating
    eq, _, fx = radiative_equilibrium(col, micro)

    print("prescribed-haze run (RT validation) ...")
    op_presc = OpticsParams(prescribed_haze=True)
    eqp, _, _ = radiative_equilibrium(col, micro, op=op_presc, n_iter=2500)

    fort = load_fortran(run_dir)

    # compare vs PRESSURE (all share the pressure coordinate)
    fig, ax = plt.subplots(1, 2, figsize=(11, 5.5))
    ax[0].plot(eq.T, eq.P, "r-", lw=1.8, label="DISORT (microphysics haze)")
    ax[0].plot(eqp.T, eqp.P, "C2-", lw=1.8, label="DISORT (prescribed haze)")
    ax[1].plot(fx.sw_heating, eq.P_mid, "C0", label="SW (DISORT)")
    ax[1].plot(fx.lw_heating, eq.P_mid, "C1", label="LW (DISORT)")

    if fort is not None:
        print(f"overlaying Fortran outputs from {run_dir} "
              f"(mean of last {fort['n_avg']}/{fort['n_total']} snapshots)")
        lbl = f"example_bowen_fort (mean of last {fort['n_avg']})"
        ax[0].plot(fort["T"], fort["P"], "k--", lw=1.5, label=lbl)
        if "T_std" in fort:
            # shade +-1 sigma: marks the unconverged (oscillating) upper atmosphere
            ax[0].fill_betweenx(fort["P"], fort["T"] - fort["T_std"],
                                fort["T"] + fort["T_std"], color="k", alpha=0.15,
                                label=r"Fortran $\pm1\sigma$ (unconverged top)")
        if "sw" in fort:
            ax[1].plot(fort["sw"], fort["P"], "k--", label="SW (Fortran)")
        if "lw" in fort:
            ax[1].plot(fort["lw"], fort["P"], "k:", label="LW (Fortran)")
        title_note = "vs reference Fortran model (TAM-derived)"
    else:
        print("Fortran outputs not found -- plotting DISORT only.")
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
