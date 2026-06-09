"""Write the reference Fortran model's prescribed-haze input tables.

The Fortran (``haze_data='presc'``) reads, for the shortwave band ``v`` and the
longwave band ``i``, three files ``<name><band>_<field>.txt`` with
``field in {tau, g, ssa}``, each in the format

    <header line>
    <nwn>                          # number of wavenumbers
    <npl>                          # number of pressure levels
    <wn[1..nwn]>                   # wavenumbers [cm^-1]
    <pl[1..npl]>                   # pressures [Pa], ascending
    <row 1 .. row npl>             # npl rows of nwn values

where ``tau`` is the *cumulative* haze optical depth (top-down), and ``ssa``,
``g`` are the *local* single-scattering albedo and asymmetry at each ``(pl, wn)``.
The Fortran bilinearly interpolates ``(log p, wn)`` onto its own model grid, so
the ``wn``/``pl`` grids are arbitrary -- we reuse the observational grids verbatim.

This module is the Step 3 interface: regenerate these tables each coupled
iteration (from the microphysics; see :func:`write_presc_haze`) and re-run the
engine.  :func:`observational_haze` returns the observational tables unchanged,
for the round-trip validation that the writer + format are correct.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

_FORT = Path(__file__).resolve().parents[2] / "src" / "example_bowen_fort"
_DATA = _FORT / "INPUT" / "DATA"
_HEADER = "statclim -- nwns, npls, wn, pl, data at each pressure level"


@dataclass
class HazeBand:
    """A haze table for one band: cumulative tau + local ssa, g on (pl, wn)."""
    wn: np.ndarray       # wavenumbers [cm^-1]            (nwn,)
    pl: np.ndarray       # pressures [Pa], ascending       (npl,)
    tau: np.ndarray      # cumulative optical depth        (npl, nwn)
    ssa: np.ndarray      # single-scattering albedo        (npl, nwn)
    g: np.ndarray        # asymmetry parameter             (npl, nwn)

    def __post_init__(self):
        npl, nwn = self.tau.shape
        for nm, a, shape in (("ssa", self.ssa, (npl, nwn)),
                             ("g", self.g, (npl, nwn))):
            if a.shape != shape:
                raise ValueError(f"{nm} shape {a.shape} != tau shape {(npl, nwn)}")
        if self.wn.size != nwn or self.pl.size != npl:
            raise ValueError("wn/pl sizes inconsistent with tau")


def parse_presc(name: str, band: str, field: str, data_dir=_DATA) -> tuple:
    """Parse ``<name><band>_<field>.txt`` -> (wn[nwn], pl[npl], data[npl, nwn])."""
    lines = (Path(data_dir) / f"{name}{band}_{field}.txt").read_text().split("\n")
    nwn, npl = int(lines[1]), int(lines[2])
    wn = np.array(lines[3].split(), float)
    pl = np.array(lines[4].split(), float)
    data = np.array([lines[5 + i].split() for i in range(npl)], float)
    return wn, pl, data


def _write_field(path: Path, wn, pl, data, header=_HEADER):
    """Write one preschaze-format file. ``data`` is (npl, nwn)."""
    npl, nwn = data.shape
    if wn.size != nwn or pl.size != npl:
        raise ValueError("wn/pl sizes inconsistent with data")
    with open(path, "w") as f:
        f.write(header + "\n")
        f.write(f"{nwn}\n{npl}\n")
        f.write(" ".join(f"{x:.8g}" for x in wn) + "\n")
        f.write(" ".join(f"{x:.8g}" for x in pl) + "\n")
        for row in data:
            f.write(" ".join(f"{x:.7e}" for x in row) + "\n")


def write_presc_haze(name: str, sw: HazeBand, lw: HazeBand, data_dir=_DATA):
    """Write the 6 prescribed-haze files for the shortwave (``sw``) and longwave
    (``lw``) bands, named ``<name>{v,i}_{tau,ssa,g}.txt`` in ``data_dir``.

    Point the Fortran at them with ``haze_presc_file = '<name>'`` in the namelist
    (use a name other than ``preschaze`` to leave the observational tables intact).
    Returns the list of written paths.
    """
    out = []
    for band, hb in (("v", sw), ("i", lw)):
        for field, arr in (("tau", hb.tau), ("ssa", hb.ssa), ("g", hb.g)):
            p = Path(data_dir) / f"{name}{band}_{field}.txt"
            _write_field(p, hb.wn, hb.pl, arr)
            out.append(p)
    return out


def _column_on_pl(atm, pl):
    """Build a :class:`rt.column.Column` whose levels are the pressures ``pl``
    (ascending in pressure = top->surface), ordered ascending in altitude."""
    from rt.column import Column
    P_asc = pl[::-1]                                  # descending P = ascending z
    zg = atm.z
    Pg = atm.pressure(zg)
    order = np.argsort(np.log(Pg))
    z_asc = np.interp(np.log(P_asc), np.log(Pg)[order], zg[order])
    return Column(z_asc, atm.temperature(z_asc), P_asc, atm.gravity(z_asc))


def microphysics_haze(micro, atm, op=None, data_dir=_DATA) -> tuple[HazeBand, HazeBand]:
    """Prescribed-haze tables driven by the Step 2 microphysics ``micro``.

    Mirrors the DISORT coupled run's haze exactly: the optical-depth *amount* is
    the mean-field (RDG) aggregate cross-section times the monomer column
    ``n*Nbar`` (``rt.optics.haze_band_tau`` in 'rdg' mode, same ``haze_abs_scale``
    and monomer radius), evaluated on the observational ``wn``/``pl`` grids and
    cumulated top-down.  The single-scattering albedo and asymmetry keep the
    observational spectral shape (a composition property the microphysics does not
    change); the longwave is a pure absorber (ssa=g=0).  Returns (sw, lw) bands
    ready for :func:`write_presc_haze`.
    """
    from rt.optics import haze_band_tau, OpticsParams, spectral_haze_sw
    op = op or OpticsParams(haze_mode="rdg")
    obs_sw, obs_lw = observational_haze(data_dir)
    col = _column_on_pl(atm, obs_sw.pl)               # pl shared by both bands

    def cumulate(tau_layer, nwn):
        # tau_layer (nwn, nlyr), ascending z (surface->top).  Cumulative-from-top
        # at each pl level (top->surface): tau_pl[0]=0 (top), tau_pl[-1]=column total.
        revcum = np.cumsum(tau_layer[:, ::-1], axis=1)[:, ::-1]        # sum i..top
        cum_full = np.concatenate([revcum, np.zeros((nwn, 1))], axis=1)
        return cum_full[:, ::-1].T                                    # (npl, nwn)

    om = spectral_haze_sw(tuple(np.round(obs_sw.wn, 3)))[0]
    tau_v = cumulate(haze_band_tau(col, micro, op, obs_sw.wn, "v", om), obs_sw.wn.size)
    tau_i = cumulate(haze_band_tau(col, micro, op, obs_lw.wn, "i"), obs_lw.wn.size)
    sw = HazeBand(wn=obs_sw.wn, pl=obs_sw.pl, tau=tau_v, ssa=obs_sw.ssa, g=obs_sw.g)
    lw = HazeBand(wn=obs_lw.wn, pl=obs_lw.pl, tau=tau_i, ssa=obs_lw.ssa, g=obs_lw.g)
    return sw, lw


def observational_haze(data_dir=_DATA) -> tuple[HazeBand, HazeBand]:
    """The observational ``preschaze`` tables as (sw, lw) :class:`HazeBand`s.

    Used for the round-trip test: writing these back out and running the engine
    must reproduce the stock prescribed-haze result, proving the writer/format.
    """
    bands = {}
    for b in ("v", "i"):
        wn, pl, tau = parse_presc("preschaze", b, "tau", data_dir)
        _, _, ssa = parse_presc("preschaze", b, "ssa", data_dir)
        _, _, g = parse_presc("preschaze", b, "g", data_dir)
        bands[b] = HazeBand(wn=wn, pl=pl, tau=tau, ssa=ssa, g=g)
    return bands["v"], bands["i"]
