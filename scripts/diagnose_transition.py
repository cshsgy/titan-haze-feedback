#!/usr/bin/env python3
"""Diagnose the coupled loop's ~137 K transition (Step 3).

The coupled loop oscillates between a warm (~183 K) and cool (~131 K) stratopause
and shows a near-discontinuity: a ~1 K change in the input T flips the engine
output ~48 K.  The composite map is T -> haze (microphysics) -> T (RT).  This
script isolates the FIRST half, cheaply and without the Fortran: sweep a family
of input temperature profiles spanning the stratopause range the loop visited and
ask whether the microphysics haze responds SMOOTHLY or with a cliff.

  smooth haze(T)   => the bistability is in the radiative transfer (multiple
                      RT-convective equilibria), not the microphysics.
  cliff in haze(T) => a microphysics regime change drives it (physical or a
                      mapping artifact to chase further).

    python3 scripts/diagnose_transition.py
Writes writing/figs/transition_diagnosis.png
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
from coupling import microphysics_haze

VIS_WN = 20000.0       # ~0.5 um


def scaled_atm(base, scale):
    """titan_reference with the whole stratospheric temperature *excess* above the
    cold point scaled by ``scale`` (troposphere/surface fixed).  Sweeping scale
    spans a weak/cool stratopause (scale<1) to a strong/warm one (scale>1) --
    covering the loop's warm (~183 K) and cool (~131 K) branches and the 137 K
    transition between them."""
    z = base.z
    T = base.temperature(z)
    Tcold = float(T.min())
    above = z > z[int(np.argmin(T))]                # stratosphere (above cold point)
    Tnew = T.copy()
    Tnew[above] = Tcold + scale * (T[above] - Tcold)
    return Atmosphere.from_profile(z, Tnew, P_surf=base.P_surf)


def haze_diagnostics(micro, atm):
    """Column visible tau, the tau-weighted centroid pressure, and microphysics
    moments for one profile."""
    sw, lw = microphysics_haze(micro, atm)
    pl = sw.pl
    bi = int(np.argmin(np.abs(sw.wn - VIS_WN)))
    cum = sw.tau[:, bi]                              # cumulative vis tau on pl
    col = float(cum[-1])
    dtau = np.diff(cum)                              # per-layer
    pmid = np.sqrt(pl[1:] * pl[:-1])
    centroid = float(np.exp(np.sum(dtau * np.log(pmid)) / max(np.sum(dtau), 1e-30)))
    # microphysics moments interpolated to a reference altitude (~1 Pa source)
    o = np.argsort(micro.z)
    zsrc = float(np.interp(np.log(1.0), np.log(atm.pressure(micro.z[o]))[::-1],
                           micro.z[o][::-1]))
    Nbar1 = float(np.interp(zsrc, micro.z[o], micro.Nbar[o]))
    _trapz = getattr(np, "trapezoid", getattr(np, "trapz"))
    rho_col = float(_trapz(micro.rho_h[o], micro.z[o]))
    return col, centroid, Nbar1, rho_col, (pl, cum)


def main():
    base = Atmosphere.titan_reference()
    # scale spans a cool/weak (~128 K) to warm/strong (~225 K) stratopause; above
    # ~230 K the microphysics BVP overflows (a separate numerical limit, not the
    # loop transition), so stop below it.
    amps = np.linspace(0.55, 1.5, 26)               # stratospheric-excess scale
    rows, profiles = [], []
    for a in amps:
        atm = scaled_atm(base, a)
        micro = solve_bvp_profile(atm, DEFAULT, n_nodes=200)
        Tstrat = float(atm.temperature(atm.z).max())
        col, cen, Nbar1, rho_col, prof = haze_diagnostics(micro, atm)
        if not (np.isfinite(col) and col < 1e3):     # skip any BVP blowup
            continue
        rows.append((a, Tstrat, col, cen, Nbar1, rho_col))
        profiles.append((Tstrat, prof))
    rows = np.array(rows)

    print("  amp[K]  stratopause  visTau  centroid[Pa]   Nbar@1Pa   rho_col")
    for r in rows:
        print(f"  {r[0]:6.1f}  {r[1]:10.1f}  {r[2]:6.2f}  {r[3]:11.2f}  {r[4]:9.1f}  {r[5]:.2e}")

    # smoothness: largest fractional jump between adjacent (sorted by stratopause T)
    o = np.argsort(rows[:, 1])
    Ts = rows[o, 1]
    for name, col in (("visTau", 2), ("centroid", 3), ("Nbar@1Pa", 4), ("rho_col", 5)):
        y = rows[o, col]
        dy = np.abs(np.diff(y)) / np.maximum(np.abs(y[:-1]), 1e-30)
        dT = np.diff(Ts)
        j = int(np.argmax(dy))
        print(f"  {name:10s}: max adjacent change {100*dy[j]:.1f}% over dT={dT[j]:.1f} K "
              f"near stratopause {Ts[j]:.0f} K")

    # figure
    fig, ax = plt.subplots(1, 3, figsize=(16, 5))
    ax[0].plot(Ts, rows[o, 2], "o-"); ax[0].set_xlabel("input stratopause T [K]")
    ax[0].set_ylabel("haze column visible tau"); ax[0].set_title("haze amount vs T")
    ax[0].axvspan(135, 139, color="orange", alpha=0.2, label="loop transition ~137 K")
    ax[0].legend(); ax[0].grid(alpha=0.3)
    ax[1].plot(Ts, rows[o, 3], "o-"); ax[1].set_xlabel("input stratopause T [K]")
    ax[1].set_ylabel("haze tau-weighted centroid [Pa]"); ax[1].invert_yaxis()
    ax[1].set_yscale("log"); ax[1].set_title("haze altitude vs T")
    ax[1].axvspan(135, 139, color="orange", alpha=0.2); ax[1].grid(alpha=0.3)
    for Tstrat, (pl, cum) in profiles[::5]:
        ax[2].plot(cum, pl, label=f"{Tstrat:.0f} K")
    ax[2].set_yscale("log"); ax[2].set_ylim(1.5e5, 1.0); ax[2].set_xlabel("cumulative vis tau")
    ax[2].set_ylabel("pressure [Pa]"); ax[2].set_title("haze tau(P) per input T")
    ax[2].legend(fontsize=7); ax[2].grid(alpha=0.3)
    fig.suptitle("Microphysics haze response to input temperature "
                 "(is T->haze smooth across ~137 K?)")
    fig.tight_layout()
    out = ROOT / "writing" / "figs" / "transition_diagnosis.png"
    fig.savefig(out, dpi=130)
    print(f"\nwrote {out.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
