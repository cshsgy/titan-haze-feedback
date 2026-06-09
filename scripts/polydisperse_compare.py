#!/usr/bin/env python3
"""Polydisperse vs monodisperse haze: optics + radiative equilibrium (pieces 4-5).

Compares the new bimodal 2-moment (polydisperse) haze against the monodisperse
model, on two things:
  4. the haze optical-depth profile (the polydisperse effect enters via the
     monomer-volume column shaped by gravitational sorting + the S/F bimodality);
  5. the DISORT radiative-convective equilibrium each haze drives, and whether the
     absorbing-haze bistability persists (run from a warm and a cool start).

    .rtenv/bin/python scripts/polydisperse_compare.py
Writes writing/figs/polydisperse_compare.png
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
from microphysics.scaling_law_bimodal import solve_bimodal_kzero, to_micro
from rt.column import Column
from rt.optics import OpticsParams, haze_band_tau, spectral_haze_sw
from rt.correlated_k import CorrelatedKSW
from rt.energy_balance import radiative_equilibrium

VIS = 20000.0     # ~0.5 um


def haze_column_tau(col, micro, ck, op):
    om = spectral_haze_sw(tuple(np.round(ck.bands, 3)))[0]
    tau = haze_band_tau(col, micro, op, ck.bands, "v", om)     # (nband, nlyr)
    bi = int(np.argmin(np.abs(ck.bands - VIS)))
    return tau[bi]                                              # per-layer vis tau


def main():
    atm = Atmosphere.titan_reference()
    mono = solve_bvp_profile(atm, DEFAULT, n_nodes=200)
    bi = to_micro(solve_bimodal_kzero(atm, DEFAULT, 1.5, 2.0, n_out=200, n_quad=14))
    col = Column.from_atmosphere(atm, nlyr=100, p_top=1.0)
    ck = CorrelatedKSW()
    op = OpticsParams(haze_mode="rdg")

    tau_m = haze_column_tau(col, mono, ck, op)
    tau_b = haze_column_tau(col, bi, ck, op)
    zk = col.z_mid / 1e3
    print("=== visible haze optical depth ===")
    print(f"  column tau  monodisperse {tau_m.sum():.2f}   bimodal {tau_b.sum():.2f}")
    # tau-weighted centroid altitude
    cen_m = float(np.sum(tau_m * zk) / tau_m.sum())
    cen_b = float(np.sum(tau_b * zk) / tau_b.sum())
    print(f"  tau centroid altitude  mono {cen_m:.0f} km   bimodal {cen_b:.0f} km")

    print("\n=== DISORT radiative-convective equilibrium ===")
    def eq(micro, label):
        a2 = atm
        c = Column.from_atmosphere(a2, nlyr=100, p_top=1.0)
        e, h, _ = radiative_equilibrium(c, micro, op=op)
        print(f"  {label}: stratopause {e.T.max():.1f} K @ {e.P[np.argmax(e.T)]:.1f} Pa "
              f"(resid {h[-1]:.2f})")
        return e
    eqm = eq(mono, "monodisperse haze")
    eqb = eq(bi, "bimodal haze")

    # bistability check with the bimodal haze: warm vs cool start, same haze
    print("\n=== bistability with the bimodal haze (warm vs cool start) ===")
    def eq_from(scale, label):
        z = atm.z; T = atm.temperature(z); Tc = T.min()
        above = z > z[int(np.argmin(T))]
        Tn = T.copy(); Tn[above] = Tc + scale * (T[above] - Tc)
        a2 = Atmosphere.from_profile(z, Tn, P_surf=atm.P_surf)
        c = Column.from_atmosphere(a2, nlyr=100, p_top=1.0)
        e, h, _ = radiative_equilibrium(c, bi, op=op, n_iter=4000)
        print(f"  {label}: stratopause {e.T.max():.1f} K (resid {h[-1]:.2f})")
        return e.T.max()
    Tw = eq_from(1.35, "warm start"); Tcl = eq_from(0.60, "cool start")
    print(f"  warm vs cool: {Tw:.0f} vs {Tcl:.0f} K -> "
          f"{'BISTABLE (persists)' if abs(Tw-Tcl) > 10 else 'monostable'}")

    # figure
    fig, ax = plt.subplots(1, 3, figsize=(15, 5.5))
    ax[0].plot(tau_m, col.P_mid, label="monodisperse"); ax[0].plot(tau_b, col.P_mid, label="bimodal")
    ax[0].set_yscale("log"); ax[0].set_ylim(1.5e5, 1.0); ax[0].set_xlabel("per-layer vis tau")
    ax[0].set_ylabel("pressure [Pa]"); ax[0].set_title("haze optical depth"); ax[0].legend(); ax[0].grid(alpha=0.3)
    ax[1].plot(mono.Nbar, mono.z / 1e3, label="mono Nbar")
    ax[1].plot(bi.Nbar, bi.z / 1e3, label="bimodal Nbar")
    ax[1].set_xscale("log"); ax[1].set_xlabel(r"$\bar N$"); ax[1].set_ylabel("altitude [km]")
    ax[1].set_title("mean size"); ax[1].legend(); ax[1].grid(alpha=0.3)
    ax[2].plot(eqm.T, eqm.P, label="mono haze"); ax[2].plot(eqb.T, eqb.P, label="bimodal haze")
    ax[2].set_yscale("log"); ax[2].set_ylim(1.5e5, 1.0); ax[2].set_xlabel("T [K]")
    ax[2].set_ylabel("pressure [Pa]"); ax[2].set_title("DISORT equilibrium"); ax[2].legend(); ax[2].grid(alpha=0.3)
    fig.suptitle("Polydisperse (bimodal) vs monodisperse haze")
    fig.tight_layout()
    out = ROOT / "writing" / "figs" / "polydisperse_compare.png"
    fig.savefig(out, dpi=130)
    print(f"\nwrote {out.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
