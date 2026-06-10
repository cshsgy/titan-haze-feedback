#!/usr/bin/env python3
"""Summarize/plot the damped branch-tracking (continuation) solves.

Reads writing/figs/continuation_<tag>.npz written by continuation_solve.py and
shows, per closure, the coupled fixed point each branch converges to: two
distinct fixed points = the coupled system is bistable; same = monostable.

    .rtenv/bin/python scripts/plot_continuation.py
Writes writing/figs/continuation_branches.png
"""

from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
FIGS = ROOT / "writing" / "figs"

CASES = [("mono", "Monodisperse"),
         ("bi1.2", r"Polydisperse $\sigma_F=1.2$"),
         ("bi2", r"Polydisperse $\sigma_F=2.0$")]


def load(tag):
    f = FIGS / f"continuation_{tag}.npz"
    return np.load(f) if f.exists() else None


fig, ax = plt.subplots(2, 3, figsize=(15, 9))
print(f"{'case':22s} {'warm fp':>9s} {'cool fp':>9s} {'split':>7s}  verdict")
for j, (pre, title) in enumerate(CASES):
    dW, dC = load(f"{pre}_warm"), load(f"{pre}_cool")
    a, b = ax[0, j], ax[1, j]
    if dW is None or dC is None:
        a.set_title(f"{title}\n(missing)")
        continue
    sW, sC = float(dW["T"].max()), float(dC["T"].max())
    split = abs(sW - sC)
    okW, okC = bool(dW["converged"]), bool(dC["converged"])
    verdict = ("BISTABLE" if split > 5 else
               "weakly bistable" if split > 2 else "MONOSTABLE")
    conv = "" if (okW and okC) else "  [UNCONVERGED CHAIN]"
    print(f"{title:22s} {sW:7.1f} K {sC:7.1f} K {split:5.1f} K  {verdict}{conv}")

    # converged profiles
    a.plot(dW["T"], dW["P"], "C3", lw=2, label=f"warm-tracked fp: {sW:.0f} K"
           + ("" if okW else " (unconv)"))
    a.plot(dC["T"], dC["P"], "C0", lw=2, label=f"cool-tracked fp: {sC:.0f} K"
           + ("" if okC else " (unconv)"))
    a.set_yscale("log"); a.set_ylim(1.5e5, 1.0); a.set_xlim(110, 200)
    a.set_title(f"{title}\n{verdict}: split {split:.1f} K")
    a.set_xlabel("T [K]")
    a.grid(alpha=0.3); a.legend(fontsize=8, loc="lower left")

    # stratopause history of the damped iteration
    for d, c, lab in ((dW, "C3", "warm"), (dC, "C0", "cool")):
        h = d["hist"]
        b.plot(h[:, 0], h[:, 1], f"{c}o-", ms=3, label=f"{lab} start: T_k")
        b.plot(h[:, 0], h[:, 2], f"{c}s--", ms=3, alpha=0.5, label=f"{lab}: F(T_k)")
    b.set_xlabel("outer iteration"); b.set_ylabel("stratopause T [K]")
    b.grid(alpha=0.3); b.legend(fontsize=7)
ax[0, 0].set_ylabel("pressure [Pa]")
fig.suptitle("Damped branch-tracking solve of the coupled fixed point "
             r"$T=F(\mathrm{haze}(T))$: converged branches (top) and iteration histories (bottom)")
fig.tight_layout()
out = FIGS / "continuation_branches.png"
fig.savefig(out, dpi=130)
print(f"wrote {out.relative_to(ROOT)}")
