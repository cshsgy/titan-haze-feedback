#!/usr/bin/env python3
"""Per-source opacity breakdown: our DISORT optics vs the reference Fortran.

At the SAME state (the Fortran's converged P/T snapshot + the prescribed haze
both models use), compute our per-source CUMULATIVE optical depth (top-down) by
band and overlay the Fortran's own opacity diagnostics, so we can see WHICH
opacity source diverges and WHERE -- the question behind the mid-altitude
heating-rate discrepancy.

Sources compared (cumulative tau vs pressure, per representative band):
  shortwave : haze (prescribed), CH4 correlated-k gas, Rayleigh
  longwave  : haze (IR), gas lines (CH4/C2H2/C2H6/C2H4/HCN), CIA

The Fortran writes per-band cumulative tau for haze (tau_hazev/tau_hazei.txt),
CIA (tau_ciai.txt) and the Rayleigh coefficient (model_interp_ray.txt) on its
203-level half-layer grid; gas tau is NOT dumped by default (see
write_gas_tau_diagnostics in the driver to enable tau_gasv/tau_gasi.txt).

    .rtenv/bin/python scripts/opacity_breakdown.py
Writes writing/figs/opacity_breakdown.png
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
from rt.column import Column
from rt.optics import (OpticsParams, haze_band_tau, rayleigh_band_tau,
                       spectral_haze_sw)
from rt.correlated_k import CorrelatedKSW, CorrelatedKLW
from rt.cia import CIABands, Composition

FORT = ROOT / "src" / "example_bowen_fort"


# ----------------------------------------------------------------------------- Fortran state + diagnostics
def fortran_snapshot():
    """Fortran converged (P[Pa], T[K]) layer profile (last finite snapshot)."""
    T = np.loadtxt(FORT / "temperatures.txt")
    P = T[0]
    valid = [i for i in range(1, T.shape[0]) if np.all(np.isfinite(T[i]))]
    return P, T[valid[-1]]


def _read_rows(path):
    """Read a Fortran diagnostic file as a list of 1-D rows (ragged-safe)."""
    rows = []
    with open(path) as fh:
        for line in fh:
            v = line.split()
            if not v:
                continue
            try:
                rows.append(np.array([float(x) for x in v]))
            except ValueError:
                continue                       # skip trailing 'Closed'/text lines
    return rows


def fortran_haze_v():
    """Cumulative SW haze tau: returns (pressure[Pa], tau[nband, nlev])."""
    r = _read_rows(FORT / "tau_hazev.txt")
    P = r[0]                                   # row0 = level pressures
    tau = np.array(r[1:])                      # rows 1.. = bands
    return P, tau


def fortran_haze_i():
    """Cumulative LW haze tau: (pressure, wno[cm^-1], tau[nband, nlev])."""
    r = _read_rows(FORT / "tau_hazei.txt")
    P, wno, tau = r[0], r[1], np.array(r[2:2 + len(r[1])])
    return P, wno, tau


def fortran_cia_i():
    """Cumulative LW CIA tau: (pressure, wno[cm^-1], tau[nband, nlev])."""
    r = _read_rows(FORT / "tau_ciai.txt")
    wno, P = r[0], r[1]
    tau = np.array([row for row in r[2:] if row.size == P.size])[:len(wno)]
    return P, wno, tau


def fortran_rayleigh():
    """Rayleigh coefficient: (wno[cm^-1], coeff[nband])."""
    r = _read_rows(FORT / "model_interp_ray.txt")
    return r[0], r[1]


def fortran_gas(band):
    """Gauss-weighted per-LAYER gas tau: (pressure[Pa], tau[nband, nlev]).

    band='v' shortwave (CH4), band='i' longwave (CH4/C2H2/C2H6/C2H4/HCN).
    Layout (row0 = plev, then one row per band) matches the tau_gas*.txt dump
    added to the Fortran driver.
    """
    r = _read_rows(FORT / f"tau_gas{band}.txt")
    P = r[0]
    tau = np.array([row for row in r[1:] if row.size == P.size])
    return P, tau


# ----------------------------------------------------------------------------- our column on the Fortran state
def build_column(P, T):
    o = np.argsort(P)[::-1]                     # surface -> top
    P, T = P[o], T[o]
    z = np.zeros_like(P)
    for i in range(1, P.size):
        Tm = 0.5 * (T[i] + T[i - 1])
        g = gravity(z[i - 1])
        z[i] = z[i - 1] + (C.R_GAS * Tm / (0.028 * g)) * np.log(P[i - 1] / P[i])
    return Column(z, T, P, gravity(z))          # ascending altitude


def cum_from_top(tau_layer):
    """Per-layer tau[nband, nlyr] (ascending) -> cumulative-from-top at the LOWER
    interface of each layer [nband, nlyr], i.e. tau above that layer's base."""
    return np.cumsum(tau_layer[:, ::-1], axis=1)[:, ::-1]


def gauss_mean(tau_bgl, gw):
    """Collapse correlated-k tau[nband, ngauss, nlyr] to a per-band representative
    via the Gauss-weighted mean (a clean 'typical' opacity for the breakdown)."""
    w = np.asarray(gw, float)
    return (tau_bgl * w[None, :, None]).sum(axis=1) / w.sum()


def main():
    P, T = fortran_snapshot()
    col = build_column(P, T)
    micro = solve_bvp_profile(Atmosphere.titan_reference(), DEFAULT, n_nodes=200)
    op = OpticsParams(prescribed_haze=True)
    comp = Composition()
    ck, cklw, cia = CorrelatedKSW(), CorrelatedKLW(), CIABands()

    Pmid = col.P_mid                                          # ascending? -> surface..top
    # ---- our per-source per-band per-layer tau, then cumulative-from-top ----
    o0, _ = spectral_haze_sw(tuple(np.round(ck.bands, 3)))
    sw = {
        "haze":     haze_band_tau(col, micro, op, ck.bands, "v", o0),     # (42, nlyr)
        "gas CH4":  gauss_mean(ck.layer_tau(col, comp), ck.gw),           # (42, nlyr)
        "Rayleigh": rayleigh_band_tau(col, ck.bands),                     # (42, nlyr)
    }
    lw_wno = 0.5 * (cia.band_lo + cia.band_hi)                            # (40,)
    lw = {
        "haze":     haze_band_tau(col, micro, op, lw_wno, "i"),           # (40, nlyr)
        "gas line": gauss_mean(cklw.layer_tau(col), cklw.gw),             # (40, nlyr)
        "CIA":      cia.optical_depth(col, comp),                         # (40, nlyr)
    }
    sw_cum = {k: cum_from_top(v) for k, v in sw.items()}
    lw_cum = {k: cum_from_top(v) for k, v in lw.items()}

    # ---- Fortran cumulative tau (interp onto our comparison pressures) ----
    fPv, fhv = fortran_haze_v()
    fPi, fwi, fhi = fortran_haze_i()
    fPc, fwc, fci = fortran_cia_i()
    rwn, rco = fortran_rayleigh()
    # Fortran gas is per-layer (top->surface order); cumulate from the top.
    fPgv, fgv = fortran_gas("v"); fgv_cum = np.cumsum(fgv, axis=1)
    fPgi, fgi = fortran_gas("i"); fgi_cum = np.cumsum(fgi, axis=1)

    def finterp(fP, ftau_band, bidx, Pq):
        """Fortran cumulative tau in band bidx, interpolated to pressures Pq."""
        s = np.argsort(fP)
        return np.interp(Pq, fP[s], ftau_band[bidx][s])

    # representative bands
    sw_bands = [("0.5 um (vis)", int(np.argmin(np.abs(ck.bands - 20000.0)))),
                ("0.9 um (CH4)", int(np.argmin(np.abs(ck.bands - 11100.0)))),
                ("2.0 um (CH4)", int(np.argmin(np.abs(ck.bands - 5000.0))))]
    lw_bands = [("100 cm-1", int(np.argmin(np.abs(lw_wno - 100.0)))),
                ("300 cm-1", int(np.argmin(np.abs(lw_wno - 300.0)))),
                ("600 cm-1", int(np.argmin(np.abs(lw_wno - 600.0))))]

    Pq = np.array([1., 3., 10., 30., 100., 300., 1000.])     # incl. mid-altitude
    o = np.argsort(Pmid)

    def our_cum(cum, bidx):
        return np.interp(Pq, Pmid[o], cum[bidx][o])

    print(f"State: Fortran snapshot, P {P.min():.2f}-{P.max():.1e} Pa, "
          f"T {T.min():.0f}-{T.max():.0f} K\n")
    print("Cumulative optical depth from TOA (ours vs Fortran).  "
          "Mid-altitude band = 10-100 Pa.\n")

    # ---------- SHORTWAVE ----------
    print("=" * 78)
    print("SHORTWAVE")
    for name, bi in sw_bands:
        fhv_b = finterp(fPv, fhv, bi, Pq)
        fgv_b = finterp(fPgv, fgv_cum, bi, Pq)
        print(f"\n  band {name}  (nu={ck.bands[bi]:.0f} cm-1)")
        print(f"  {'P[Pa]':>7} | {'haze(us)':>9} {'haze(F)':>9} | "
              f"{'gasCH4(us)':>10} {'gasCH4(F)':>10} | {'Ray(us)':>8} {'Ray(F)':>8}")
        ray_cum = our_cum(sw_cum["Rayleigh"], bi)
        haze_cum = our_cum(sw_cum["haze"], bi)
        gas_cum = our_cum(sw_cum["gas CH4"], bi)
        # Fortran Rayleigh cumulative: coeff * cumulative dp[mbar] from top
        s = np.argsort(fPv); fp = fPv[s]
        dp_mbar = np.diff(np.concatenate([[fp[0]], fp])) / 100.0
        fray_cum_full = rco[np.argmin(np.abs(rwn - ck.bands[bi]))] * np.cumsum(dp_mbar)
        fray_cum = np.interp(Pq, fp, fray_cum_full)
        for j, pq in enumerate(Pq):
            mid = " *" if 10 <= pq <= 100 else "  "
            print(f"  {pq:7.0f} | {haze_cum[j]:9.3f} {fhv_b[j]:9.3f} | "
                  f"{gas_cum[j]:10.3f} {fgv_b[j]:10.3f} | "
                  f"{ray_cum[j]:8.4f} {fray_cum[j]:8.4f}{mid}")

    # ---------- LONGWAVE ----------
    print("\n" + "=" * 78)
    print("LONGWAVE")
    for name, bi in lw_bands:
        # Fortran LW band index (its wno grid)
        bfi_h = int(np.argmin(np.abs(fwi - lw_wno[bi])))
        bfi_c = int(np.argmin(np.abs(fwc - lw_wno[bi])))
        bfi_g = int(np.argmin(np.abs(fwi - lw_wno[bi])))   # gas dump uses same wno order
        fhi_b = finterp(fPi, fhi, bfi_h, Pq)
        fci_b = finterp(fPc, fci, bfi_c, Pq)
        fgi_b = finterp(fPgi, fgi_cum, bfi_g, Pq)
        haze_cum = our_cum(lw_cum["haze"], bi)
        gas_cum = our_cum(lw_cum["gas line"], bi)
        cia_cum = our_cum(lw_cum["CIA"], bi)
        print(f"\n  band {name}  (nu={lw_wno[bi]:.0f} cm-1)")
        print(f"  {'P[Pa]':>7} | {'haze(us)':>9} {'haze(F)':>9} | "
              f"{'gas(us)':>8} {'gas(F)':>8} | {'CIA(us)':>9} {'CIA(F)':>9}")
        for j, pq in enumerate(Pq):
            mid = " *" if 10 <= pq <= 100 else "  "
            print(f"  {pq:7.0f} | {haze_cum[j]:9.3f} {fhi_b[j]:9.3f} | "
                  f"{gas_cum[j]:8.3f} {fgi_b[j]:8.3f} | "
                  f"{cia_cum[j]:9.3f} {fci_b[j]:9.3f}{mid}")

    # ---------- VERDICT: mid-altitude (10-100 Pa) relative discrepancy ----------
    print("\n" + "=" * 78)
    print("VERDICT  --  mid-altitude (10-100 Pa) cumulative-tau ratio ours/Fortran")
    print("(averaged over the sampled bands where the source is non-negligible)")
    midmask = (Pq >= 10) & (Pq <= 100)

    def ratio(our_src_cum, fP, ftau, bands, fbands):
        rs = []
        for _, bi in bands:
            fbi = bi if fbands is None else int(np.argmin(np.abs(fbands - lw_wno[bi])))
            ours = our_cum(our_src_cum, bi)[midmask]
            fort = finterp(fP, ftau, fbi, Pq)[midmask]
            m = fort > 1e-4
            if m.any():
                rs.extend((ours[m] / fort[m]).tolist())
        return (np.mean(rs), np.min(rs), np.max(rs)) if rs else (np.nan,)*3

    for label, args in [
        ("SW haze    ", (sw_cum["haze"], fPv, fhv, sw_bands, None)),
        ("SW gas CH4 ", (sw_cum["gas CH4"], fPgv, fgv_cum, sw_bands, None)),
        ("LW haze    ", (lw_cum["haze"], fPi, fhi, lw_bands, fwi)),
        ("LW gas line", (lw_cum["gas line"], fPgi, fgi_cum, lw_bands, fwi)),
        ("LW CIA     ", (lw_cum["CIA"], fPc, fci, lw_bands, fwc)),
    ]:
        mean, lo, hi = ratio(*args)
        if mean != mean:
            print(f"  {label}: (negligible in sampled bands)")
            continue
        tag = "OK" if 0.8 <= mean <= 1.25 else "<-- check"
        print(f"  {label}: ours/F = {mean:5.2f}  (range {lo:4.2f}-{hi:4.2f})  {tag}")
    print("\n  -> if all sources are ~1.0 here, the mid-altitude HEATING gap is NOT")
    print("     opacity; it is the RT solver (8-stream DISORT vs 2-stream/hemispheric).")

    # ---------- figure: cumulative tau by source vs pressure ----------
    fig, axes = plt.subplots(2, 3, figsize=(15, 9), sharey=True)
    Pp = Pmid[o]
    for ax, (name, bi) in zip(axes[0], sw_bands):
        ax.plot(our_cum(sw_cum["haze"], bi)*0 + np.nan, Pq)  # keep colour cycle
    # SW row
    for ax, (name, bi) in zip(axes[0], sw_bands):
        ax.plot(sw_cum["haze"][bi][o], Pp, "C0-", label="haze (us)")
        ax.plot(sw_cum["gas CH4"][bi][o], Pp, "C1-", label="gas CH4 (us)")
        ax.plot(sw_cum["Rayleigh"][bi][o], Pp, "C2-", label="Rayleigh (us)")
        s, sg = np.argsort(fPv), np.argsort(fPgv)
        ax.plot(fhv[bi][s], fPv[s], "C0--", label="haze (Fortran)")
        ax.plot(fgv_cum[bi][sg], fPgv[sg], "C1--", label="gas CH4 (Fortran)")
        ax.set_title(f"SW {name}")
        ax.set_xscale("symlog", linthresh=1e-3)
    # LW row
    for ax, (name, bi) in zip(axes[1], lw_bands):
        ax.plot(lw_cum["haze"][bi][o], Pp, "C0-", label="haze (us)")
        ax.plot(lw_cum["gas line"][bi][o], Pp, "C1-", label="gas line (us)")
        ax.plot(lw_cum["CIA"][bi][o], Pp, "C3-", label="CIA (us)")
        bfi_h = int(np.argmin(np.abs(fwi - lw_wno[bi])))
        bfi_c = int(np.argmin(np.abs(fwc - lw_wno[bi])))
        sh, sc, sgi = np.argsort(fPi), np.argsort(fPc), np.argsort(fPgi)
        ax.plot(fhi[bfi_h][sh], fPi[sh], "C0--", label="haze (Fortran)")
        ax.plot(fgi_cum[bfi_h][sgi], fPgi[sgi], "C1--", label="gas line (Fortran)")
        ax.plot(fci[bfi_c][sc], fPc[sc], "C3--", label="CIA (Fortran)")
        ax.set_title(f"LW {name}")
        ax.set_xscale("symlog", linthresh=1e-3)
    for ax in axes.flat:
        ax.set_yscale("log"); ax.set_ylim(1.3e5, 1.0)
        ax.grid(alpha=0.3); ax.legend(fontsize=7)
        ax.axhspan(10, 100, color="orange", alpha=0.08)      # mid-altitude band
    for ax in axes[:, 0]:
        ax.set_ylabel("pressure [Pa]")
    for ax in axes[1]:
        ax.set_xlabel("cumulative tau from TOA")
    fig.suptitle("Per-source cumulative optical depth: DISORT (solid) vs Fortran "
                 "(dashed).  Shaded = mid-altitude (10-100 Pa)")
    fig.tight_layout()
    out = ROOT / "writing" / "figs" / "opacity_breakdown.png"
    fig.savefig(out, dpi=130)
    print(f"\nwrote {out.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
