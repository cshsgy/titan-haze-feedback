#!/usr/bin/env python3
"""Step 3 rung 3: the coupled radiative<->microphysical fixed-point loop.

Iterate the feedback to self-consistency using the reference Fortran model as the
Step 1 RT engine:

    T  --(Atmosphere.from_profile)-->  microphysics (Step 2 BVP)
       --(coupling.microphysics_haze)-->  prescribed-haze tables
       --(Fortran engine)-->  new T   ... repeat (under-relaxed) to a fixed point.

Each Fortran pass is warm-started from the current T (init_t=0,
initial_temperatures.txt) so it re-equilibrates quickly.  Reports the feedback:
the coupled fixed point vs. the prescribed-haze (no-feedback) baseline, and the
stratopause trajectory.  Writes writing/figs/coupled_feedback.png and
docs/coupled_history.txt.

    .rtenv/bin/python scripts/run_coupled.py [max_iter]

Restores the namelist on exit.  This is the scientific payoff of the project;
expect a handful of Fortran runs (minutes each).
"""

import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from microphysics import Atmosphere, DEFAULT, solve_bvp_profile
from coupling import microphysics_haze, write_presc_haze

FORT = ROOT / "src" / "example_bowen_fort"
NAMELIST = FORT / "namelist"
TEMPS = FORT / "temperatures.txt"
NAME = "cplhaze"
M_AIR, G_TITAN, R = 0.028, 1.352, 8.314


def set_namelist(updates: dict):
    """Set namelist keys (preserving the rest); return the original text."""
    orig = NAMELIST.read_text()
    lines = orig.splitlines()
    for i, ln in enumerate(lines):
        key = ln.split("=")[0].strip()
        if key in updates:
            lines[i] = f"{key} = {updates[key]},"
    NAMELIST.write_text("\n".join(lines) + "\n")
    return orig


def run_engine():
    """Full equilibration from the seed (cold start) for the current haze.

    Cold-starting each pass (rather than warm-starting from the previous T) means
    every Fortran run reaches the TRUE equilibrium for its haze, so the outer-loop
    residual measures the feedback fixed point -- not a half-relaxed engine that
    barely moved from a warm start (which falsely reads as converged)."""
    subprocess.run(["bash", str(ROOT / "scripts" / "run_fortran.sh")],
                   check=True, capture_output=True)


def mean_T(navg=20):
    """Mean of the last navg finite snapshots: (P[100], T[100]) top->surface."""
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


def atm_from(P, T):
    """Atmosphere from the Fortran (P, T) grid (top->surface), extended a little
    above the top isothermally so the microphysics BVP has headroom over the
    ~1 Pa haze source."""
    o = np.argsort(P)[::-1]                         # surface -> top
    Ps, Ts = P[o], T[o]
    ext = Ps[-1] * np.array([0.3, 0.1, 0.03])
    Ps = np.concatenate([Ps, ext])
    Ts = np.concatenate([Ts, np.full(3, Ts[-1])])
    z = np.zeros_like(Ps)
    for i in range(1, Ps.size):
        Tm = 0.5 * (Ts[i] + Ts[i - 1])
        z[i] = z[i - 1] + (R * Tm / (M_AIR * G_TITAN)) * np.log(Ps[i - 1] / Ps[i])
    return Atmosphere.from_profile(z, Ts, P_surf=float(Ps[0]))


def stratopause(P, T):
    return float(T.max()), float(P[np.argmax(T)])


def main():
    max_iter = int(sys.argv[1]) if len(sys.argv) > 1 else 8
    # The feedback is strong and near-bistable (a warm ~183 K and a cool ~131 K
    # stratopause attractor), so naive iteration oscillates; under-relax hard.
    relax, tol = float(sys.argv[2]) if len(sys.argv) > 2 else 0.3, 1.5
    # full equilibration (cold start) each pass; keep the namelist's num_i
    orig = set_namelist({"haze_presc_file": "'preschaze'", "init_t": "0."})
    history = []
    try:
        # --- baseline: prescribed (observational) haze, no feedback ---
        print("baseline: prescribed-haze equilibrium (no feedback) ...")
        subprocess.run(["bash", str(ROOT / "scripts" / "run_fortran.sh")],
                       check=True, capture_output=True)
        P, T = mean_T()
        Tbase = T.copy()
        sp, spp = stratopause(P, T)
        print(f"  baseline stratopause {sp:.1f} K @ {spp:.1f} Pa, surface {T[-1]:.1f} K")

        # --- coupled iteration, warm-started ---
        set_namelist({"haze_presc_file": f"'{NAME}'"})
        profiles = [("baseline (prescribed haze)", P.copy(), Tbase.copy())]
        T_new = Tbase.copy()
        for it in range(max_iter):
            atm = atm_from(P, T)
            micro = solve_bvp_profile(atm, DEFAULT, n_nodes=200)
            sw, lw = microphysics_haze(micro, atm)
            write_presc_haze(NAME, sw, lw)
            run_engine()                             # cold start, full equilibration
            P, T_new = mean_T()
            resid = float(np.max(np.abs(T_new - T)))
            sp, spp = stratopause(P, T_new)
            history.append((it, sp, spp, float(T_new[-1]), resid))
            profiles.append((f"iter {it}", P.copy(), T_new.copy()))
            print(f"  iter {it}: stratopause {sp:.1f} K @ {spp:.1f} Pa, "
                  f"surface {T_new[-1]:.1f} K, max|dT| {resid:.2f} K")
            T = T + relax * (T_new - T)              # under-relax
            if resid < tol:
                print(f"  converged (max|dT| {resid:.2f} < {tol} K)")
                break
    finally:
        NAMELIST.write_text(orig)
        for b in ("v", "i"):
            for fld in ("tau", "ssa", "g"):
                (FORT / "INPUT" / "DATA" / f"{NAME}{b}_{fld}.txt").unlink(missing_ok=True)
        print("  (namelist + haze files restored)")

    # --- report ---
    spb, _ = stratopause(profiles[0][1], profiles[0][2])
    spc, _ = stratopause(P, T_new)
    print("\n=== feedback ===")
    print(f"  prescribed-haze stratopause : {spb:.1f} K")
    print(f"  coupled fixed-point stratopause : {spc:.1f} K")
    print(f"  feedback shift : {spc - spb:+.1f} K at the stratopause")
    print(f"  max |coupled - prescribed| over the column : {np.max(np.abs(T_new - Tbase)):.1f} K")

    # --- figure ---
    hist = np.array(history)
    fig, ax = plt.subplots(1, 2, figsize=(12, 5.5))
    for label, Pp, Tp in profiles:
        style = "k--" if label.startswith("baseline") else "-"
        lw = 2.0 if label.startswith("baseline") else 1.2
        ax[0].plot(Tp, Pp, style, lw=lw, label=label)
    ax[0].set_yscale("log"); ax[0].set_ylim(1.5e5, 1.0)
    ax[0].set_xlabel("temperature [K]"); ax[0].set_ylabel("pressure [Pa]")
    ax[0].set_title("coupled feedback: T(p) per iteration"); ax[0].legend(fontsize=7)
    ax[0].grid(alpha=0.3)
    if hist.size:
        ax[1].plot(hist[:, 0], hist[:, 1], "o-", label="stratopause T")
        ax[1].axhline(spb, color="k", ls="--", lw=1, label="prescribed (no feedback)")
        ax[1].set_xlabel("iteration"); ax[1].set_ylabel("stratopause T [K]")
        ax[1].set_title("convergence"); ax[1].legend(fontsize=8); ax[1].grid(alpha=0.3)
    fig.tight_layout()
    out = ROOT / "writing" / "figs" / "coupled_feedback.png"
    fig.savefig(out, dpi=130)
    print(f"\nwrote {out.relative_to(ROOT)}")

    hpath = ROOT / "docs" / "coupled_history.txt"
    hpath.write_text("# iter  stratopause_K  stratopause_Pa  surface_K  max_dT_K\n" +
                     "\n".join(f"{int(h[0])} {h[1]:.2f} {h[2]:.2f} {h[3]:.2f} {h[4]:.3f}"
                              for h in history) + "\n")
    print(f"wrote {hpath.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
