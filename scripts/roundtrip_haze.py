#!/usr/bin/env python3
"""End-to-end round-trip of the prescribed-haze writer (Step 3, rung 1).

Writes the OBSERVATIONAL haze back out through coupling.write_presc_haze under a
new name, points the Fortran engine at it, runs, and checks it reproduces the
stock prescribed-haze T(z).  This validates the writer + file format against the
actual Fortran reader -- the prerequisite for swapping in the microphysics haze.

    .rtenv/bin/python scripts/roundtrip_haze.py   # (system python also works)

Runs the engine twice (~minutes).  Restores the namelist on exit.
"""

import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import numpy as np

from coupling import observational_haze, write_presc_haze

FORT = ROOT / "src" / "example_bowen_fort"
NAMELIST = FORT / "namelist"
TEMPS = FORT / "temperatures.txt"
COUPLED_NAME = "rthaze"          # writes rthaze{v,i}_{tau,ssa,g}.txt


def run_fortran():
    subprocess.run(["bash", str(ROOT / "scripts" / "run_fortran.sh")],
                   check=True, capture_output=True)


def last_finite_T():
    rows = [ln.split() for ln in TEMPS.read_text().splitlines() if ln.split()]
    def parse(r):
        out = []
        for x in r:
            try:
                out.append(float(x))
            except ValueError:
                out.append(np.nan)
        return np.array(out)
    P = parse(rows[0])
    snaps = [parse(r) for r in rows[1:]]
    finite = [s for s in snaps if np.all(np.isfinite(s))]
    if not finite:
        raise RuntimeError("no finite temperature snapshot")
    return P, finite[-1]


def set_haze_presc_file(name):
    """Rewrite the namelist's haze_presc_file entry; return the original text."""
    orig = NAMELIST.read_text()
    out = []
    for ln in orig.splitlines():
        if ln.strip().startswith("haze_presc_file"):
            out.append(f"haze_presc_file = '{name}',")
        else:
            out.append(ln)
    NAMELIST.write_text("\n".join(out) + "\n")
    return orig


def main():
    print("1/3  baseline run with the stock 'preschaze' haze ...")
    orig_namelist = set_haze_presc_file("preschaze")
    try:
        run_fortran()
        P0, T0 = last_finite_T()

        print("2/3  writing the observational haze back out as "
              f"'{COUPLED_NAME}' and running with it ...")
        sw, lw = observational_haze()
        written = write_presc_haze(COUPLED_NAME, sw, lw)
        set_haze_presc_file(COUPLED_NAME)
        run_fortran()
        P1, T1 = last_finite_T()

        print("3/3  comparing ...")
        assert np.allclose(P0, P1), "pressure grids differ"
        dT = np.abs(T1 - T0)
        o = np.argsort(P0)
        print(f"  max |dT| = {dT.max():.4f} K  (at P={P0[np.argmax(dT)]:.2f} Pa)")
        print(f"  RMS |dT| = {np.sqrt((dT**2).mean()):.4f} K")
        print("  sample (P[Pa], T_presc, T_rthaze):")
        for i in o[::20]:
            print(f"    {P0[i]:9.2f}  {T0[i]:7.2f}  {T1[i]:7.2f}")
        ok = dT.max() < 0.5
        print(f"\n{'PASS' if ok else 'FAIL'}: round-trip reproduces the stock haze "
              f"T(z) to {dT.max():.3f} K (writer + format validated)")
        # clean up the round-trip haze files
        for p in written:
            p.unlink(missing_ok=True)
        return 0 if ok else 1
    finally:
        NAMELIST.write_text(orig_namelist)
        print("  (namelist restored)")


if __name__ == "__main__":
    sys.exit(main())
