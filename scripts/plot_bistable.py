#!/usr/bin/env python3
"""Plot the converged equilibria cached by bistable_states.py.

Two panels: the monodisperse haze settles to two well-separated branches
(bistable, 17 K), the polydisperse (bimodal) haze to two close branches (~8 K,
weakly bistable -- strongly suppressed, not eliminated).  Reads the npz so the
figure can be retuned without re-running DISORT.

    .rtenv/bin/python scripts/plot_bistable.py [sf]      # sf in {1.2, 2.0}; default 2.0
Writes writing/figs/bistable_states.png
"""

import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
SF = sys.argv[1] if len(sys.argv) > 1 else "2.0"
POLY = f"bimodal_sf{SF}"

d = np.load(ROOT / "writing" / "figs" / "bistable_states.npz")


def split(name):
    return abs(d[f"{name}_W_T"].max() - d[f"{name}_C_T"].max())


def panel(a, name, title):
    a.plot(d["base_T"], d["base_P"], color="0.6", ls="--", lw=1.2,
           label="prescribed-haze baseline")
    a.plot(d[f"{name}_W_T"], d[f"{name}_W_P"], color="C3", lw=2,
           label=f"warm start $\\to$ {d[f'{name}_W_T'].max():.0f} K")
    a.plot(d[f"{name}_C_T"], d[f"{name}_C_P"], color="C0", lw=2,
           label=f"cool start $\\to$ {d[f'{name}_C_T'].max():.0f} K")
    dT = split(name)
    verdict = "bistable" if dT > 12 else ("weakly bistable" if dT > 4 else "monostable")
    a.set_title(f"{title}\n{verdict}: stratopause split $\\approx{dT:.0f}$ K")
    a.set_xlabel("temperature [K]")
    a.set_yscale("log"); a.set_ylim(1.5e5, 1.0); a.set_xlim(120, 200)
    a.grid(alpha=0.3); a.legend(loc="lower left", fontsize=8)


fig, ax = plt.subplots(1, 2, figsize=(11, 5.5), sharey=True)
panel(ax[0], "monodisperse", "Monodisperse haze")
panel(ax[1], POLY, fr"Polydisperse haze ($\sigma_F={SF}$)")
ax[0].set_ylabel("pressure [Pa]")
fig.suptitle("Stable equilibria of the coupled haze feedback: polydispersity "
             "strongly suppresses, but does not eliminate, the bistability")
fig.tight_layout()
out = ROOT / "writing" / "figs" / "bistable_states.png"
fig.savefig(out, dpi=130)
print(f"split: monodisperse {split('monodisperse'):.0f} K, "
      f"bimodal_sf1.2 {split('bimodal_sf1.2'):.0f} K, "
      f"bimodal_sf2.0 {split('bimodal_sf2.0'):.0f} K")
print(f"wrote {out.relative_to(ROOT)} (polydisperse panel: sigma_F={SF})")
