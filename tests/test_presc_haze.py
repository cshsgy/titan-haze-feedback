#!/usr/bin/env python3
"""Round-trip test for the prescribed-haze writer (Step 3 coupling, rung 1).

Writing the observational haze back out and re-reading it must reproduce the
tables to the writer's precision, and the written files must be in the exact
format the Fortran reader expects.  (The end-to-end check that the Fortran reads
them and yields the same T(z) is scripts/roundtrip_haze.py, which needs the
built engine.)

    python3 tests/test_presc_haze.py
"""

import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import numpy as np

from coupling import observational_haze, write_presc_haze, parse_presc
from coupling.presc_haze import HazeBand

_PASS = 0
_FAIL = 0


def check(name, cond, info=""):
    global _PASS, _FAIL
    if cond:
        _PASS += 1
        print(f"PASS  {name}")
    else:
        _FAIL += 1
        print(f"FAIL  {name}  {info}")


def test_roundtrip_observational():
    """Observational tables -> write -> read == original (to writer precision)."""
    sw, lw = observational_haze()
    with tempfile.TemporaryDirectory() as d:
        write_presc_haze("coupledhaze", sw, lw, data_dir=d)
        for band, hb in (("v", sw), ("i", lw)):
            for field, ref in (("tau", hb.tau), ("ssa", hb.ssa), ("g", hb.g)):
                wn, pl, data = parse_presc("coupledhaze", band, field, data_dir=d)
                check(f"roundtrip {band}_{field} wn", np.allclose(wn, hb.wn, rtol=1e-6),
                      f"max dwn {np.abs(wn-hb.wn).max():.2e}")
                check(f"roundtrip {band}_{field} pl", np.allclose(pl, hb.pl, rtol=1e-6))
                rel = np.abs(data - ref) / np.maximum(np.abs(ref), 1e-30)
                check(f"roundtrip {band}_{field} data",
                      np.allclose(data, ref, rtol=1e-5, atol=1e-12),
                      f"max rel {rel.max():.2e}")


def test_format_matches_fortran_reader():
    """Written files must have the layout the Fortran read_presc expects:
    header, nwn, npl, wn-line, pl-line, then npl rows of nwn values."""
    sw, _ = observational_haze()
    with tempfile.TemporaryDirectory() as d:
        write_presc_haze("coupledhaze", sw, sw, data_dir=d)   # sw for both bands here
        lines = (Path(d) / "coupledhazev_tau.txt").read_text().splitlines()
        nwn, npl = int(lines[1]), int(lines[2])
        check("header line present", "nwns" in lines[0] or lines[0] != "")
        check("nwn matches", nwn == sw.wn.size, f"{nwn} vs {sw.wn.size}")
        check("npl matches", npl == sw.pl.size, f"{npl} vs {sw.pl.size}")
        check("wn line count", len(lines[3].split()) == nwn)
        check("pl line count", len(lines[4].split()) == npl)
        check("total rows", len(lines) >= 5 + npl)
        check("a data row has nwn values", len(lines[5].split()) == nwn)


def test_cumulative_tau_monotone():
    """Sanity: the observational tau we round-trip is cumulative (non-decreasing
    down the column) and ssa/g are in range -- guards against axis mixups."""
    sw, lw = observational_haze()
    for nm, hb in (("sw", sw), ("lw", lw)):
        # cumulative top->surface: tau increases with pressure (axis 0 ascending)
        dtau = np.diff(hb.tau, axis=0)
        check(f"{nm} tau cumulative (non-decreasing)", (dtau >= -1e-9).all(),
              f"min dtau {dtau.min():.2e}")
        check(f"{nm} ssa in [0,1]", (hb.ssa >= 0).all() and (hb.ssa <= 1.0001).all())
        check(f"{nm} g in [0,1)", (hb.g >= -1e-9).all() and (hb.g < 1.0).all())
    # longwave is a pure absorber
    check("lw ssa == 0", np.allclose(lw.ssa, 0.0))
    check("lw g == 0", np.allclose(lw.g, 0.0))


def test_microphysics_haze_matches_disort():
    """The microphysics->preschaze mapping must reproduce the DISORT coupled
    run's haze optical depth (same RDG optics), with a valid cumulative table."""
    from microphysics import Atmosphere, DEFAULT, solve_bvp_profile
    from coupling import microphysics_haze, observational_haze
    from coupling.presc_haze import _column_on_pl
    from rt.optics import haze_band_tau, OpticsParams, spectral_haze_sw

    atm = Atmosphere.titan_reference()
    micro = solve_bvp_profile(atm, DEFAULT, n_nodes=200)
    sw, lw = microphysics_haze(micro, atm)
    obs_sw, _ = observational_haze()

    check("mphaze sw shape", sw.tau.shape == obs_sw.tau.shape)
    check("mphaze sw cumulative", (np.diff(sw.tau, axis=0) >= -1e-9).all())
    check("mphaze sw top row ~0", np.allclose(sw.tau[0], 0.0))
    check("mphaze lw pure absorber", np.allclose(lw.ssa, 0) and np.allclose(lw.g, 0))
    # column-total (surface row) == direct DISORT haze_band_tau column sum
    col = _column_on_pl(atm, obs_sw.pl)
    om = spectral_haze_sw(tuple(np.round(obs_sw.wn, 3)))[0]
    direct = haze_band_tau(col, micro, OpticsParams(haze_mode="rdg"),
                           obs_sw.wn, "v", om).sum(axis=1)
    check("mphaze sw total == DISORT haze sum",
          np.allclose(sw.tau[-1], direct, rtol=1e-9),
          f"max rel {np.abs(sw.tau[-1]-direct).max():.2e}")
    # ssa/g keep the observational shape
    check("mphaze sw ssa == observational", np.allclose(sw.ssa, obs_sw.ssa))


if __name__ == "__main__":
    test_roundtrip_observational()
    test_format_matches_fortran_reader()
    test_cumulative_tau_monotone()
    test_microphysics_haze_matches_disort()
    print(f"\n{_PASS}/{_PASS + _FAIL} passed")
    sys.exit(1 if _FAIL else 0)
