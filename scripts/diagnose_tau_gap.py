#!/usr/bin/env python3
"""Phase 0 of closing the 40 K stratopause gap: WHERE is the haze tau deficit?

The coupled fixed point sits at ~142 K vs ~183 K with the observational
(prescribed) haze, with the SAME RT -- so the difference is the haze.  The mono
model matches the observed column (tau 8.4 vs 8) and mean scale height (64 vs
65 km), so the deficit must be in the vertical placement.  Compare cumulative
visible-SW and LW haze tau(P), observational vs model hazes, and quantify the
deficit factor above the stratopause (P < ~100 Pa).

    .rtenv/bin/python scripts/diagnose_tau_gap.py
Writes writing/figs/tau_gap.png
"""

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams.update({"font.size": 13, "axes.titlesize": 14, "axes.labelsize": 13,
                     "legend.fontsize": 11, "lines.linewidth": 2.2,
                     "figure.titlesize": 15})

from microphysics import Atmosphere, DEFAULT, solve_bvp_profile
from microphysics.scaling_law_bimodal import solve_bimodal_kzero, to_micro
from coupling.presc_haze import observational_haze, microphysics_haze
from rt.optics import OpticsParams

# This script DIAGNOSES the lab-tholin configuration, so the LW haze must use
# the Khare RDG-IR (the discrepancy being shown), not the obs-calibrated default.
OP_KHARE = OpticsParams(haze_mode="rdg", lw_haze="khare")

VIS = 20000.0     # cm^-1, ~0.5 um
LWB = 700.0       # cm^-1, mid-IR


def at_wn(band, wn0):
    """Cumulative tau(P) at the band's wn closest to wn0 -> (pl, tau)."""
    j = int(np.argmin(np.abs(band.wn - wn0)))
    return band.pl, band.tau[:, j], band.wn[j]


def main():
    atm = Atmosphere.titan_reference()
    obs_sw, obs_lw = observational_haze()

    hazes = {"mono BVP": solve_bvp_profile(atm, DEFAULT, n_nodes=200),
             "bimodal sf1.2": to_micro(solve_bimodal_kzero(atm, DEFAULT, 1.5, 1.2,
                                                           n_out=200, n_quad=14)),
             "bimodal sf2.0": to_micro(solve_bimodal_kzero(atm, DEFAULT, 1.5, 2.0,
                                                           n_out=200, n_quad=14))}
    model = {k: microphysics_haze(m, atm, op=OP_KHARE) for k, m in hazes.items()}

    fig, ax = plt.subplots(1, 3, figsize=(15, 5.5), sharey=True)
    for a, (band_i, wn0, label) in zip(
            ax[:2], [(0, VIS, "visible SW"), (1, LWB, "LW")]):
        pl, t_obs, wn = at_wn((obs_sw, obs_lw)[band_i], wn0)
        a.plot(t_obs, pl, "k-", lw=2, label="observational")
        for k, (msw, mlw) in model.items():
            _, t, _ = at_wn((msw, mlw)[band_i], wn0)
            a.plot(t, pl, label=k)
        a.set_xscale("log"); a.set_yscale("log"); a.set_ylim(1.5e5, 1.0)
        a.set_xlim(1e-4, 30)
        a.set_xlabel(f"cumulative tau ({label}, {wn:.0f} cm$^{{-1}}$)")
        a.grid(alpha=0.3); a.legend(fontsize=8)
    ax[0].set_ylabel("pressure [Pa]")

    # deficit factor vs pressure (visible)
    pl, t_obs, _ = at_wn(obs_sw, VIS)
    a = ax[2]
    for k, (msw, _) in model.items():
        _, t, _ = at_wn(msw, VIS)
        a.plot(t_obs / np.maximum(t, 1e-12), pl, label=k)
    a.axvline(1, color="k", lw=0.8)
    a.set_xscale("log"); a.set_yscale("log"); a.set_ylim(1.5e5, 1.0)
    a.set_xlim(0.05, 100)
    a.set_xlabel("deficit factor  tau_obs / tau_model (visible)")
    a.grid(alpha=0.3); a.legend(fontsize=8)
    fig.suptitle("Model (lab-tholin LW) vs observational haze opacity, cumulative from TOA")
    fig.tight_layout()
    out = ROOT / "writing" / "figs" / "tau_gap.png"
    fig.savefig(out, dpi=130)

    print("=== cumulative visible tau from TOA, observational vs models ===")
    rows = [1.0, 10.0, 100.0, 1000.0, 10000.0, float(pl.max())]
    hdr = "P [Pa]      obs   " + "  ".join(f"{k:>13s}" for k in model)
    print(hdr)
    for P0 in rows:
        i = int(np.argmin(np.abs(pl - P0)))
        vals = []
        for k, (msw, _) in model.items():
            _, t, _ = at_wn(msw, VIS)
            vals.append(t[i])
        t_o = t_obs[i]
        print(f"{pl[i]:9.1f}  {t_o:6.3f}  " +
              "  ".join(f"{v:8.3f} (x{t_o/max(v,1e-12):5.1f})" for v in vals))
    print(f"\nwrote {out.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
