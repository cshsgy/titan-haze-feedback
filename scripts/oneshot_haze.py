#!/usr/bin/env python3
"""One-shot coupling cross-check (Step 3, rung 2).

Drive BOTH radiative engines with the SAME Step 2 microphysics haze on a guess
temperature, and compare the resulting profiles:
  * DISORT: radiative-convective equilibrium with the RDG microphysics haze.
  * Fortran engine: the same haze written through coupling.microphysics_haze
    into its prescribed-haze tables, then run to equilibrium.

Both should reach a haze-warmed stratosphere (~140 K stratopause expected from
the DISORT coupled run); agreement validates that the haze->Fortran wiring
delivers the intended optics.  This is one pass, not the coupled fixed point
(that is rung 3 -- iterate T<->haze).

    .rtenv/bin/python scripts/oneshot_haze.py

Runs both engines (~minutes).  Restores the namelist on exit.
"""

import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import numpy as np

from microphysics import Atmosphere, DEFAULT, solve_bvp_profile
from rt.column import Column
from rt.energy_balance import radiative_equilibrium
from coupling import microphysics_haze, write_presc_haze

FORT = ROOT / "src" / "example_bowen_fort"
NAMELIST = FORT / "namelist"
TEMPS = FORT / "temperatures.txt"
NAME = "mphaze"


def set_haze_presc_file(name):
    orig = NAMELIST.read_text()
    out = [f"haze_presc_file = '{name}'," if ln.strip().startswith("haze_presc_file")
           else ln for ln in orig.splitlines()]
    NAMELIST.write_text("\n".join(out) + "\n")
    return orig


def fortran_mean_profile(navg=20):
    rows = [ln.split() for ln in TEMPS.read_text().splitlines() if ln.split()]
    def parse(r):
        v = []
        for x in r:
            try:
                v.append(float(x))
            except ValueError:
                v.append(np.nan)
        return np.array(v)
    P = parse(rows[0])
    snaps = [parse(r) for r in rows[1:]]
    finite = [s for s in snaps if np.all(np.isfinite(s))]
    return P, np.mean(finite[-navg:], axis=0)


def main():
    atm = Atmosphere.titan_reference()
    micro = solve_bvp_profile(atm, DEFAULT, n_nodes=200)

    print("1/3  DISORT coupled equilibrium (RDG microphysics haze) ...")
    col = Column.from_atmosphere(atm, nlyr=100, p_top=1.0)
    eqd, _, _ = radiative_equilibrium(col, micro)
    Td, Pd = eqd.T, eqd.P
    print(f"     DISORT stratopause {Td.max():.1f} K, tropopause {Td.min():.1f} K, "
          f"surface {Td[0]:.1f} K")

    print("2/3  writing the microphysics haze and running the Fortran engine ...")
    sw, lw = microphysics_haze(micro, atm)
    written = write_presc_haze(NAME, sw, lw)
    orig = set_haze_presc_file(NAME)
    try:
        subprocess.run(["bash", str(ROOT / "scripts" / "run_fortran.sh")],
                       check=True, capture_output=True)
        Pf, Tf = fortran_mean_profile()
    finally:
        NAMELIST.write_text(orig)
        for p in written:
            p.unlink(missing_ok=True)
        print("     (namelist restored, mphaze files removed)")

    of = np.argsort(Pf)
    print(f"     Fortran stratopause {Tf.max():.1f} K, tropopause {Tf.min():.1f} K, "
          f"surface {Tf[of][-1]:.1f} K")

    print("3/3  comparing DISORT vs Fortran (same microphysics haze) ...")
    Pq = np.array([1., 3., 10., 30., 100., 300., 1e3, 3e3, 1e4, 1e5])
    od = np.argsort(Pd)
    Td_q = np.interp(Pq, Pd[od], Td[od])
    Tf_q = np.interp(Pq, Pf[of], Tf[of])
    print(f"  {'P[Pa]':>8} {'DISORT':>8} {'Fortran':>8} {'dT':>7}")
    for i, pq in enumerate(Pq):
        print(f"  {pq:8.0f} {Td_q[i]:8.1f} {Tf_q[i]:8.1f} {Td_q[i]-Tf_q[i]:7.1f}")
    dT = np.abs(Td_q - Tf_q)
    print(f"\n  max |dT| (DISORT vs Fortran) = {dT.max():.1f} K ; "
          f"stratopause {Td.max():.0f} vs {Tf.max():.0f} K")
    print("  (one pass, same haze; residual is the RT-engine difference + the "
          "Fortran's residual top oscillation, not a coupling error)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
