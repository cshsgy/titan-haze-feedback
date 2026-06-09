"""Collision-induced absorption (CIA) for Titan's thermal-IR opacity.

Parses HITRAN CIA files (Karman et al. 2019; hitran.org/cia), band-averages them
onto the longwave band grid, and computes per-layer CIA optical depth from the
pair number-density products.  On Titan the far-IR continuum is dominated by
N2-N2, with N2-CH4, CH4-CH4, and N2-H2 contributing.

HITRAN CIA convention: the tabulated coefficient k(nu, T) has units
cm^5 molecule^-2, so the absorption coefficient is

    alpha(nu) [cm^-1] = sum_pairs  k_pair(nu, T) * n_A * n_B          (n in cm^-3)

and the layer optical depth is alpha * dz.  The compact band-averaged table is
cached in data/cia/titan_cia_bands.npz; run scripts/fetch_cia.py to (re)build it
from the raw HITRAN files.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

_trapz = getattr(np, "trapezoid", getattr(np, "trapz", None))  # numpy 2.x rename

_ROOT = Path(__file__).resolve().parents[2]
_RAW = _ROOT / "data" / "cia" / "raw"
_NPZ = _ROOT / "data" / "cia" / "titan_cia_bands.npz"

# raw HITRAN files per pair (see scripts/fetch_cia.py)
_RAW_FILES = {
    "N2-N2": "N2-N2_2021.cia",
    "N2-CH4": "N2-CH4_2024.cia",
    "N2-H2": "N2-H2_2024.cia",
    "CH4-CH4": "CH4-CH4_2011.cia",
}

# default longwave band edges [cm^-1]: match the IR correlated-k gas bands
# (bwni: 0-2000 cm^-1 in 50 cm^-1 bands, centres 25..1975) so CIA and the gas
# k-tables share one band grid.
DEFAULT_LW_EDGES = np.arange(0.0, 2001.0, 50.0)

# Titan-relevant temperature grid for the table [K]
DEFAULT_T_GRID = np.arange(70.0, 230.0, 10.0)


@dataclass
class Composition:
    """Uniform mole fractions (simplification; CH4 actually rises near surface)."""
    x_CH4: float = 0.014
    x_H2: float = 1.0e-3

    @property
    def x_N2(self) -> float:
        return 1.0 - self.x_CH4 - self.x_H2


# ----------------------------------------------------------------------------
# parsing + band-averaging (build step)
# ----------------------------------------------------------------------------
def _parse_blocks(path: Path, pair: str):
    """Yield (T, nu_array, k_array) blocks from a HITRAN .cia file."""
    with open(path) as f:
        lines = f.readlines()
    i, n = 0, len(lines)
    while i < n:
        tok = lines[i].split()
        if tok and tok[0] == pair:
            T = float(tok[4])
            npts = int(tok[3])
            data = np.array([row.split() for row in lines[i + 1:i + 1 + npts]],
                            dtype=float)
            i += 1 + npts
            yield T, data[:, 0], data[:, 1]
        else:
            i += 1


def _band_average(nu, k, edges):
    """Average k over each band [edges[j], edges[j+1]] (k=0 outside data range)."""
    out = np.zeros(edges.size - 1)
    for j in range(edges.size - 1):
        lo, hi = edges[j], edges[j + 1]
        m = (nu >= lo) & (nu <= hi)
        if m.sum() >= 2:
            out[j] = _trapz(k[m], nu[m]) / (hi - lo)
        elif m.sum() == 1:
            out[j] = k[m][0]
    return out


def build_band_table(edges=None, T_grid=None, out=_NPZ):
    """Build and cache the band-averaged CIA table k_pair[nband, nT]."""
    edges = DEFAULT_LW_EDGES if edges is None else np.asarray(edges, float)
    T_grid = DEFAULT_T_GRID if T_grid is None else np.asarray(T_grid, float)
    table = {"edges": edges, "T_grid": T_grid}
    for pair, fname in _RAW_FILES.items():
        path = _RAW / fname
        if not path.exists():
            print(f"  [skip] {pair}: {path} missing")
            continue
        # keep the first block at each rounded temperature (main dataset)
        byT = {}
        for T, nu, k in _parse_blocks(path, pair):
            key = round(T)
            if key not in byT:
                byT[key] = (nu, k)
        kband = np.zeros((edges.size - 1, T_grid.size))
        for it, T in enumerate(T_grid):
            key = int(round(T))
            # nearest available temperature block
            avail = np.array(sorted(byT))
            kkey = int(avail[np.argmin(np.abs(avail - key))])
            nu, k = byT[kkey]
            kband[:, it] = _band_average(nu, k, edges)
        table[pair] = kband
        print(f"  {pair}: table {kband.shape}, peak {kband.max():.3e} cm^5/molec^2")
    out.parent.mkdir(parents=True, exist_ok=True)
    np.savez(out, **table)
    print(f"  wrote {out}")
    return out


# ----------------------------------------------------------------------------
# runtime use
# ----------------------------------------------------------------------------
class CIABands:
    """Band-averaged CIA table with per-layer optical-depth evaluation."""

    def __init__(self, npz=_NPZ):
        if not Path(npz).exists():
            raise FileNotFoundError(
                f"{npz} not found; run scripts/fetch_cia.py to build it.")
        d = np.load(npz)
        self.edges = d["edges"]
        self.band_lo = self.edges[:-1]
        self.band_hi = self.edges[1:]
        self.T_grid = d["T_grid"]
        self.pairs = {k: d[k] for k in d.files if k not in ("edges", "T_grid")}
        self.nband = self.band_lo.size

    def _k_at_T(self, pair, T):
        """Interpolate band-averaged k[nband] to temperatures T[nlyr]."""
        kb = self.pairs[pair]                       # (nband, nT)
        T = np.clip(T, self.T_grid[0], self.T_grid[-1])
        return np.array([np.interp(T, self.T_grid, kb[b]) for b in range(self.nband)])

    def optical_depth(self, column, comp: Composition | None = None):
        """CIA optical depth tau[nband, nlyr] for the column layers."""
        from microphysics.constants import K_B
        comp = comp or Composition()
        # number density per layer [molecules / cm^3]
        n_tot = (column.P_mid / (K_B * column.T_mid)) * 1e-6
        n_N2 = comp.x_N2 * n_tot
        n_CH4 = comp.x_CH4 * n_tot
        n_H2 = comp.x_H2 * n_tot
        prod = {
            "N2-N2": n_N2 * n_N2,
            "N2-CH4": n_N2 * n_CH4,
            "CH4-CH4": n_CH4 * n_CH4,
            "N2-H2": n_N2 * n_H2,
        }
        dz_cm = column.dz * 100.0                   # [cm]
        T = column.T_mid
        alpha = np.zeros((self.nband, column.nlyr))  # [cm^-1]
        for pair, kpair in self.pairs.items():
            if pair in prod:
                alpha += self._k_at_T(pair, T) * prod[pair][None, :]
        return alpha * dz_cm[None, :]


# ----------------------------------------------------------------------------
# Fortran-matching exponential-sum CIA
# ----------------------------------------------------------------------------
# The reference model does not band-average HITRAN cross-sections; it evaluates a
# 3-term exponential-sum fit to the band *transmission* (INPUT/DATA/trans_*.txt):
#
#     T(band) = sum_i a_i(band,Temp) * exp(-k_i(band,Temp) * u),   tau = -ln T
#
# with the absorber amount u = n1 * n2 * dist * 1e-7 in the model's own scaled
# units: number density n = vmr * (P/(T*Runiv)) * avo with avo = 6.02214086 (the
# 1e23 dropped, absorbed into the fit's k_i), and dist the layer thickness in km.
# The fit captures sub-band absorber-amount nonlinearity a single band-averaged
# cross-section cannot, and uses different underlying data; porting it here makes
# our CIA identical to the reference's by construction.
_TRANS_DIR = _ROOT / "src" / "example_bowen_fort" / "INPUT" / "DATA"
_TRANS_FILES = {
    "N2-N2": ("N2", "N2"),
    "N2-CH4": ("N2", "CH4"),
    "CH4-CH4": ("CH4", "CH4"),
    "N2-H2": ("N2", "H2"),
}
_RUNIV = 8.3144598          # universal gas constant [J/mol/K] (Fortran Runiv)
_AVO_SCALED = 6.02214086    # Avogadro / 1e23  (Fortran avo; the 1e23 is in the fit)


def _parse_trans(path):
    """Parse a trans_<pair>.txt exp-sum file.

    Returns (wno[nband], temps[ntemp], coeff[ntemp, nband, 2*nfits]); the first
    nfits coefficients are the weights a_i, the next nfits the exponents k_i.
    """
    lines = [ln for ln in Path(path).read_text().splitlines() if ln.strip()]
    wno = np.array([float(x) for x in lines[1].split()])     # line 0 is a header
    nb = wno.size
    temps, blocks, i = [], [], 2
    while i < len(lines):
        temps.append(float(lines[i].split()[0])); i += 1
        blocks.append(np.array([[float(x) for x in lines[i + b].split()]
                                for b in range(nb)]))
        i += nb
    return wno, np.array(temps), np.array(blocks)


class CIAExpSum:
    """Exp-sum-transmission CIA, matching the reference Fortran's ``get_tauCIA``.

    Same ``optical_depth(column, comp) -> tau[nband, nlyr]`` interface and band
    grid as :class:`CIABands`, so it is a drop-in replacement.
    """

    def __init__(self, trans_dir=_TRANS_DIR, pairs=_TRANS_FILES):
        self.tables, wno = {}, None
        for pair, species in pairs.items():
            w, temps, coeff = _parse_trans(Path(trans_dir) / f"trans_{pair}.txt")
            self.tables[pair] = (species, temps, coeff)
            wno = w
        self.wno = wno
        self.band_lo = wno - 25.0                 # 50 cm^-1 bands, centres 25..1975
        self.band_hi = wno + 25.0
        self.nband = wno.size
        self.nfits = next(iter(self.tables.values()))[2].shape[2] // 2

    def _pair_tau(self, temps, coeff, T, u):
        """Per-pair tau[nband, nlyr] by linear-in-T interpolation of the band
        transmission (exactly as the Fortran: interpolate T, not the coeffs)."""
        nf = self.nfits
        a, k = coeff[:, :, :nf], coeff[:, :, nf:]      # (ntemp, nband, nfits)

        def trans_at(it, uu):                          # transmission at temp index it
            return (a[it] * np.exp(-k[it] * uu)).sum(axis=1)   # (nband,)

        out = np.zeros((self.nband, T.size))
        for l in range(T.size):
            Tl, ul = float(T[l]), float(u[l])
            if Tl <= temps[0]:
                tr = trans_at(0, ul)
            elif Tl >= temps[-1]:
                tr = trans_at(-1, ul)
            else:
                it = int(np.searchsorted(temps, Tl))   # temps[it-1] < Tl <= temps[it]
                tlo, thi = trans_at(it - 1, ul), trans_at(it, ul)
                tr = thi - (temps[it] - Tl) / (temps[it] - temps[it - 1]) * (thi - tlo)
            tr = np.clip(tr, 1e-307, 1.0)
            out[:, l] = -np.log(tr)
        return out

    def optical_depth(self, column, comp: Composition | None = None):
        comp = comp or Composition()
        rho_moles = column.P_mid / (column.T_mid * _RUNIV)     # mol/m^3
        nd = {s: x * rho_moles * _AVO_SCALED for s, x in
              (("N2", comp.x_N2), ("CH4", comp.x_CH4), ("H2", comp.x_H2))}
        dist_km = column.dz / 1000.0
        tau = np.zeros((self.nband, column.nlyr))
        for pair, (species, temps, coeff) in self.tables.items():
            s1, s2 = species
            u = nd[s1] * nd[s2] * dist_km * 1e-7               # col_abund (scaled units)
            tau += self._pair_tau(temps, coeff, column.T_mid, u)
        return tau
