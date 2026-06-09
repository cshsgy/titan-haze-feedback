#!/usr/bin/env python3
"""Confirm the solver hypothesis: vary the DISORT stream count and see whether
our heating moves toward the reference Fortran's two-stream result.

The opacity breakdown (scripts/opacity_breakdown.py) showed every opacity source
matches the Fortran within a few percent through the mid-atmosphere, so the
residual heating-rate difference must be the radiative-transfer SOLVER: our
eight-stream discrete ordinates vs the Fortran's delta-Eddington two-stream
(shortwave) / hemispheric-mean (longwave).

This script runs the SAME optics (same state, same prescribed haze) through
DISORT at nstr = 2, 4, 8, 16 and compares each to the Fortran heating at the
same snapshot.  If nstr=2 (the closest analogue to the Fortran's two-stream)
matches the Fortran better than nstr=8 -- especially at mid-altitude where the
opacities agree -- the solver hypothesis is confirmed.

    .rtenv/bin/python scripts/solver_streams.py
Writes writing/figs/solver_streams.png
"""

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from microphysics import Atmosphere, DEFAULT, solve_bvp_profile, constants as C
from microphysics.atmosphere import gravity
from rt.column import Column, TITAN_DAY
from rt.optics import OpticsParams
from rt.correlated_k import CorrelatedKSW, CorrelatedKLW
from rt import energy_balance as eb

FORT = ROOT / "src" / "example_bowen_fort"
STREAMS = [2, 4, 8, 16]


def fortran_state():
    """Last finite T snapshot + its sw/lw heating [K/s] per layer."""
    T = np.loadtxt(FORT / "temperatures.txt")
    P = T[0]
    sw = np.loadtxt(FORT / "sw.txt")
    lw = np.loadtxt(FORT / "lw.txt")
    valid = [i for i in range(1, T.shape[0]) if np.all(np.isfinite(T[i]))]
    k = valid[-1]
    return P, T[k], sw[k - 1], lw[k - 1]


def build_column(P, T):
    o = np.argsort(P)[::-1]
    P, T = P[o], T[o]
    z = np.zeros_like(P)
    for i in range(1, P.size):
        Tm = 0.5 * (T[i] + T[i - 1])
        z[i] = z[i - 1] + (C.R_GAS * Tm / (0.028 * gravity(z[i - 1]))) * np.log(P[i - 1] / P[i])
    return Column(z, T, P, gravity(z))


def main():
    P, Tf, sw_f, lw_f = fortran_state()
    col = build_column(P, Tf)
    micro = solve_bvp_profile(Atmosphere.titan_reference(), DEFAULT, n_nodes=200)
    op = OpticsParams(prescribed_haze=True)               # the haze the Fortran uses
    ck, cklw = CorrelatedKSW(), CorrelatedKLW()
    solar = eb.SolarForcing(umu0=1.0 / np.pi)              # match the Fortran insolation

    # our heating [K/Titan-day] at each stream count, on the Fortran snapshot
    runs = {}
    for n in STREAMS:
        fx = eb.compute_fluxes(col, micro, op=op, solar=solar,
                               ck=ck, ck_lw=cklw, nstr=n)
        runs[n] = (fx.sw_heating, fx.lw_heating)
    Pmid = col.P_mid
    o = np.argsort(Pmid)

    # comparison pressures (K/Titan-day); Fortran is K/s -> convert
    Pq = np.array([1., 3., 10., 30., 100., 300., 1000., 3000.])
    of = np.argsort(P)

    def to_q(arr_layer, Pgrid, order):
        return np.interp(Pq, Pgrid[order], arr_layer[order])

    swf_q = to_q(sw_f * TITAN_DAY, P, of)
    lwf_q = to_q(lw_f * TITAN_DAY, P, of)

    print(f"State: Fortran snapshot, prescribed haze, umu0=1/pi.  "
          f"Heating in K/Titan-day.\n")
    print("Does reducing the stream count move our heating toward the Fortran?\n")

    for band, idx, fq in [("SHORTWAVE heating", 0, swf_q),
                          ("LONGWAVE  heating", 1, lwf_q)]:
        print("=" * 72)
        print(f"{band}   (mid-altitude 10-100 Pa flagged *)")
        hdr = "  " + f"{'P[Pa]':>7} |" + "".join(f"{('n='+str(n)):>9}" for n in STREAMS) \
              + f"{'Fortran':>10}"
        print(hdr)
        q = {n: to_q(runs[n][idx], Pmid, o) for n in STREAMS}
        for j, pq in enumerate(Pq):
            mid = " *" if 10 <= pq <= 100 else "  "
            row = f"  {pq:7.0f} |" + "".join(f"{q[n][j]:9.2f}" for n in STREAMS) \
                  + f"{fq[j]:10.2f}{mid}"
            print(row)
        # mid-altitude closeness to Fortran per stream count
        m = (Pq >= 10) & (Pq <= 100)
        print("  mid-altitude mean |ours - Fortran|:")
        for n in STREAMS:
            print(f"    nstr={n:2d}: {np.mean(np.abs(q[n][m] - fq[m])):.2f} K/day")
        best = min(STREAMS, key=lambda n: np.mean(np.abs(q[n][m] - fq[m])))
        print(f"  -> closest to Fortran at mid-altitude: nstr={best}\n")

    # ---- figure: heating profiles per stream count vs Fortran ----
    fig, ax = plt.subplots(1, 2, figsize=(12, 6), sharey=True)
    for a, (title, idx, ff, Pf) in zip(
            ax, [("shortwave", 0, sw_f, P), ("longwave", 1, lw_f, P)]):
        for n in STREAMS:
            a.plot(runs[n][idx][o], Pmid[o], label=f"DISORT nstr={n}")
        a.plot(ff[of] * TITAN_DAY, Pf[of], "k--", lw=1.8, label="Fortran (2-stream)")
        a.set_title(f"{title} heating")
        a.set_yscale("log"); a.set_ylim(1.3e5, 1.0)
        a.set_xlim(-45, 10) if idx == 1 else a.set_xlim(-2, 95)  # focus mid-altitude
        a.axhspan(10, 100, color="orange", alpha=0.1)
        a.axvline(0, color="grey", lw=0.6); a.grid(alpha=0.3); a.legend(fontsize=8)
        a.set_xlabel("rate [K / Titan-day]")
    ax[0].set_ylabel("pressure [Pa]")
    fig.suptitle("Solver test: heating vs DISORT stream count, against the "
                 "Fortran two-stream (shaded = mid-altitude 10-100 Pa)")
    fig.tight_layout()
    out = ROOT / "writing" / "figs" / "solver_streams.png"
    fig.savefig(out, dpi=130)
    print(f"wrote {out.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
