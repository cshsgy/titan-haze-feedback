#!/usr/bin/env python3
"""sigma_F sweep: how the aggregate distribution width controls haze opacity and
the absorbing-haze bistability (polydisperse follow-up).

Pieces 4-5 found that at sigma_F=2.0 the bimodal haze is ~4.4x optically thinner
than monodisperse (gravitational sorting) and the radiative-convective
equilibrium becomes monostable.  Both effects are sigma_F-dependent.  This sweep:
  (a) cheap optics scan -- visible column tau and centroid vs sigma_F, to find the
      width that recovers the observed column tau ~ 8;
  (b) DISORT bistability (warm vs cool start, same haze) at a few sigma_F, to map
      where the two-equilibria behaviour appears.

    .rtenv/bin/python scripts/sigma_f_sweep.py
Writes writing/figs/sigma_f_sweep.png
"""

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from microphysics import Atmosphere, DEFAULT
from microphysics.scaling_law_bimodal import solve_bimodal_kzero, to_micro
from rt.column import Column
from rt.optics import OpticsParams, haze_band_tau, spectral_haze_sw
from rt.correlated_k import CorrelatedKSW
from rt.energy_balance import radiative_equilibrium

VIS = 20000.0
TAU_OBS = 8.0


def main():
    atm = Atmosphere.titan_reference()
    col = Column.from_atmosphere(atm, nlyr=100, p_top=1.0)
    ck = CorrelatedKSW()
    op = OpticsParams(haze_mode="rdg")
    om = spectral_haze_sw(tuple(np.round(ck.bands, 3)))[0]
    bj = int(np.argmin(np.abs(ck.bands - VIS)))
    zk = col.z_mid / 1e3

    def haze(sf):
        return to_micro(solve_bimodal_kzero(atm, DEFAULT, 1.5, sf, n_out=200, n_quad=14))

    # ---- (a) optics sweep (cheap) ----
    sig = np.linspace(1.2, 2.2, 11)
    tau, cen = [], []
    for sf in sig:
        t = haze_band_tau(col, haze(sf), op, ck.bands, "v", om)[bj]
        tau.append(float(t.sum())); cen.append(float(np.sum(t * zk) / t.sum()))
    tau = np.array(tau); cen = np.array(cen)
    sf_obs = float(np.interp(TAU_OBS, tau[::-1], sig[::-1]))   # sigma_F giving tau~8
    print("=== optics sweep: visible column tau vs sigma_F ===")
    for s, t, c in zip(sig, tau, cen):
        print(f"  sigma_F={s:.2f}  tau={t:5.2f}  centroid={c:5.0f} km")
    print(f"  -> tau~{TAU_OBS:.0f} (observed) at sigma_F ~ {sf_obs:.2f}")

    # ---- (b) DISORT bistability at bracketing sigma_F ----
    def split(sf):
        mic = haze(sf)
        def eq(scale):
            z = atm.z; T = atm.temperature(z); Tc = T.min()
            above = z > z[int(np.argmin(T))]
            Tn = T.copy(); Tn[above] = Tc + scale * (T[above] - Tc)
            a2 = Atmosphere.from_profile(z, Tn, P_surf=atm.P_surf)
            c = Column.from_atmosphere(a2, nlyr=70, p_top=1.0)
            e, h, _ = radiative_equilibrium(c, mic, op=op, n_iter=2000)
            return e.T.max()
        Tw, Tcl = eq(1.35), eq(0.60)
        return Tw, Tcl

    sig_b = [1.2, 1.6, 2.0]                        # thick -> thin haze bracket
    print("\n=== bistability (warm vs cool start, same haze) vs sigma_F ===")
    res = []
    for sf in sig_b:
        Tw, Tcl = split(sf)
        d = abs(Tw - Tcl)
        res.append((sf, Tw, Tcl, d))
        print(f"  sigma_F={sf:.2f}: warm {Tw:.0f} K, cool {Tcl:.0f} K -> "
              f"split {d:.0f} K  {'BISTABLE' if d > 10 else 'monostable'}")
    res = np.array(res)

    # ---- figure ----
    fig, ax = plt.subplots(1, 2, figsize=(12, 5))
    ax[0].plot(sig, tau, "o-"); ax[0].axhline(TAU_OBS, color="k", ls="--", label=r"observed $\tau\approx8$")
    ax[0].axvline(sf_obs, color="C1", ls=":", label=fr"$\sigma_F\approx{sf_obs:.2f}$")
    ax[0].set_xlabel(r"$\sigma_F$"); ax[0].set_ylabel("visible column $\\tau$")
    ax[0].set_title("haze opacity vs distribution width"); ax[0].legend(); ax[0].grid(alpha=0.3)
    ax[1].plot(res[:, 0], res[:, 3], "s-")
    ax[1].axhline(10, color="grey", ls=":", label="bistability threshold")
    ax[1].set_xlabel(r"$\sigma_F$"); ax[1].set_ylabel("warm$-$cool stratopause split [K]")
    ax[1].set_title("bistability vs distribution width"); ax[1].legend(); ax[1].grid(alpha=0.3)
    fig.suptitle(r"$\sigma_F$ controls haze opacity and the absorbing-haze bistability")
    fig.tight_layout()
    out = ROOT / "writing" / "figs" / "sigma_f_sweep.png"
    fig.savefig(out, dpi=130)
    print(f"\nwrote {out.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
