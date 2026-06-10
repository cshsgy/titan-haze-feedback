# Step 3 — Coupling scope: Step 2 haze → Fortran RT

How to drive the Step 1 engine (`src/example_bowen_fort/`) with the Step 2
microphysics haze, and iterate the two to a self-consistent fixed point.

## Architecture: Python orchestration, file-based coupling

The Fortran reads its haze once at init from the **prescribed-haze path**
(`haze_data='presc'`), so we couple **without touching the Fortran source**: a
Python loop regenerates the prescribed-haze files from the microphysics each
outer iteration and re-runs the model.

```
T0  (guess, e.g. titan_reference)
repeat:
    atm   = Atmosphere.from_profile(z, T)         # T from the Fortran (or T0 on iter 0)
    micro = solve_bvp_profile(atm)                # Step 2 -> n(z), Nbar(z), ...
    write_presc_haze(micro, atm)                  # -> coupledhaze{v,i}_{tau,g,ssa}.txt
    run_fortran(seed_T=T)                          # Step 1 -> new T(z), P(z)
    T = read_fortran_T()
until  max|ΔT| < tol
```

The microphysics rates depend on T (`β ∝ T/η(T)`, settling, free-molecular
`T^{1/2}`), so this is a genuine feedback fixed point, not a one-shot.

## The interface the Fortran expects

`haze_init` → `read_presc` reads, per band `v` (shortwave) and `i` (longwave),
three files `INPUT/DATA/<haze_presc_file><band>_<field>.txt` with
`field ∈ {tau, g, ssa}`. Format (see `preschazev_tau.txt`):

```
<header line>
<nwn>                      # 43 for v, 41 for i
<npl>                      # 600
<wn[1..nwn]>               # wavenumbers [cm^-1]
<pl[1..npl]>               # pressures [Pa], ascending (0.0177 .. 1.46e5)
<row 1: nwn values @ pl[1]>
...
<row npl: nwn values @ pl[npl]>
```

- **`tau`** = *cumulative* haze optical depth, top-down (the Fortran differences
  it to per-layer `dtau`). **`ssa`, `g`** = *local* values at each `(pl, wn)`.
- The grid is independent of the model layers — the Fortran bilinearly
  interpolates `(log p, wn)` onto its own `wnov/wnoi × plev`. **Reuse the existing
  `pl` grid and `wn` grids verbatim** so the interpolation behaves identically.
- LW haze is a pure absorber: write `ssa=0, g=0` for band `i`.

Set `haze_presc_file = 'coupledhaze'` in the namelist so we write
`coupledhaze*.txt` and leave the observational `preschaze*.txt` untouched.

## What we already have (reuse)

| need | exists as |
|---|---|
| Atmosphere from arbitrary `T(z)` | `Atmosphere.from_profile(z, T, P_surf)` |
| microphysics haze profile | `solve_bvp_profile(atm) -> BVPResult(z, n, Nbar, rho_h, r, r_a, omega)` |
| per-band haze τ from microphysics (RDG) | `aggregate_optics.aggregate_haze_layer_tau(column, micro, wn_cm, omega0_band, d_mono, pure_absorber)` |
| spectral ω₀(λ), g(λ) | `optics.spectral_haze_sw(...)` (observational; composition-only, amount set by microphysics) |
| monomer radius `d` | `microphysics.constants.AerosolParams.d_mono` (50 nm) |
| read the Fortran T snapshot | `compare_fortran.load_fortran` / `rt_diagnostics.fortran_state` (last finite snapshot) |

So the haze **amount** (τ) comes from the microphysics RDG cross-section × `n·Nbar`;
the **scattering** (ω₀, g) stays the observational spectral shape — exactly the
mapping our DISORT coupled run already uses.

## New code (small)

1. **`src/coupling/presc_haze.py`** — `write_presc_haze(micro, atm, out_dir, name='coupledhaze')`:
   - reuse the existing `wn` grids (43 SW / 41 LW) and `pl` grid (600) read from
     the observational files;
   - build a `Column` on `pl` (hydrostatic z from `atm`);
   - SW: `aggregate_haze_layer_tau(col, micro, wn_v, omega0_v, d_mono)` + `spectral_haze_sw` → τ, ssa, g;
   - LW: `aggregate_haze_layer_tau(col, micro, wn_i, None, d_mono, pure_absorber=True)`, ssa=g=0;
   - cumulate τ top-down; write the 6 files in the exact format.
2. **`scripts/run_coupled.py`** — the outer loop above; under-relax `T` between
   iterations; log `max|ΔT|` and the stratopause; save the fixed-point profile +
   a convergence plot.
3. **`scripts/run_fortran.sh`** tweak (or `run_fortran_coupled.sh`) — seed
   `initial_temperatures.txt` from the previous iteration's T (`init_t=1`) for a
   warm start, and optionally a smaller `num_i` once warm.
4. **namelist**: `haze_presc_file='coupledhaze'`, `init_t=1` (warm start).

## Validation ladder

1. **Round-trip** — DONE. `coupling.presc_haze.write_presc_haze` writes the
   prescribed-haze tables; `tests/test_presc_haze.py` (33 checks) validates the
   format, and `scripts/roundtrip_haze.py` feeds the *observational* haze back
   through the writer → Fortran → reproduces the stock `preschaze` `T(z)` to
   **0.000 K**. Writer + format validated against the real reader.
2. **One-shot** — DONE. `coupling.microphysics_haze(micro, atm)` maps the Step 2
   profile to the tables (RDG τ from `n·N̄` on the SW/LW wn grids, cumulated;
   ω₀/g keep the observational spectral shape; LW pure absorber) and reproduces
   the DISORT coupled run's haze τ to `rtol=1e-9` (`tests/test_presc_haze.py`).
   `scripts/oneshot_haze.py` drives BOTH engines with that haze: DISORT
   stratopause 142 K, Fortran 131 K (tropopause 72 K both); max |ΔT| 15.5 K in
   the lower stratosphere — the RT-engine difference (8-stream vs two-stream) +
   the Fortran's residual top oscillation, not a coupling error. The wiring
   delivers the intended optics.
3. **Iterate** — DONE (`scripts/run_coupled.py`). Closes the loop: prescribed-haze
   baseline, then `Fortran T → Atmosphere.from_profile → solve_bvp_profile →
   microphysics_haze → rewrite tables → cold-start re-run`, under-relaxed, with a
   feedback report + `writing/figs/coupled_feedback.png` + `docs/coupled_history.txt`.
   Each pass cold-starts (full equilibration) so the residual reflects the
   feedback, not a half-relaxed engine (an early warm-start version with reduced
   `num_i` falsely "converged" at the baseline).

## Step 3 result: a strong, near-bistable feedback

The coupled iteration does **not** settle to a single fixed point under naive
under-relaxation. It oscillates between two attractors:

- **warm**: stratopause $\sim\!183\,\mathrm{K}$ at $8\,\mathrm{Pa}$ (close to the
  prescribed-haze baseline, $184.5\,\mathrm{K}$);
- **cool**: stratopause $\sim\!131$–$135\,\mathrm{K}$ at $1.3$–$1.7\,\mathrm{Pa}$
  (the microphysics-haze state; matches the one-shot's $131\,\mathrm{K}$).

They are **anti-correlated** — a warm $T$ makes a haze that radiates to a cool
equilibrium and vice-versa, i.e. period-2 oscillation. With `relax=0.5` the loop
damped toward $\sim\!133\,\mathrm{K}$ (iters 2–5, `max|dT|` $60\to6\,\mathrm{K}$)
but then jumped back to the warm state — a near-discontinuity around
$T\approx137\,\mathrm{K}$ (a $1.3\,\mathrm{K}$ input change flipped the output
$48\,\mathrm{K}$). So the radiative$\leftrightarrow$microphysical feedback is
**strong and destabilising** (the stratopause-warming haze feedback has gain
$>1$ near the cool branch), with a sharp transition between branches.

### Diagnosis of the ~137 K transition

Is the jump a microphysics regime change, a mapping artifact, or a real
radiative feedback? Two cheap experiments (no Fortran) settle it:

1. **Microphysics is smooth** (`scripts/diagnose_transition.py`,
   `writing/figs/transition_diagnosis.png`). Sweeping a family of input
   temperature profiles (stratopause $128$–$225\,\mathrm{K}$) and computing the
   microphysics haze for each, the column $\tau$, the $\tau$-weighted centroid,
   and `rho_h` all vary **smoothly** — largest adjacent change $\le3.5\%$, with
   **no cliff at $137\,\mathrm{K}$**. (The BVP only overflows above
   $\sim\!230\,\mathrm{K}$, a separate numerical limit.) The map $T\to$ haze is
   continuous; the bistability is not a microphysics artifact. The haze centroid
   rises smoothly with $T$ ($1714\to770\,\mathrm{Pa}$) — a positive feedback
   (warmer $\Rightarrow$ haze higher $\Rightarrow$ absorbs higher $\Rightarrow$ warmer).

2. **The radiative transfer is genuinely bistable** (`scripts/rt_multiplicity.py`).
   Holding the haze **fixed** (the transition-state microphysics haze) and
   relaxing DISORT to radiative--convective equilibrium from a warm vs.\ a cool
   start lands on **two different equilibria** — stratopause
   $158\,\mathrm{K}$ vs.\ $141\,\mathrm{K}$ (max $|\Delta T|\approx31\,\mathrm{K}$).
   The gap is stable as the residual falls ($1.84\to0.73$, $2000\to4000$
   iterations), so it is real multiplicity, not incomplete convergence. DISORT is
   an independent engine from the Fortran, so this is engine-independent.

**Conclusion.** The $\sim\!137\,\mathrm{K}$ transition is a **genuine radiative
feedback**: the strongly-absorbing haze, smoothly mobile in altitude with
temperature, gives the 1-D radiative--convective system two stable states
(warm/high-stratopause and cool/low-stratopause), with gain $>1$ between them. It
is not a microphysics regime change or a coupling-mapping artifact.

**Next steps.** A converged coupled branch needs harder damping (`relax<=0.25`,
now the default) or a continuation / damped-Newton solve of
$T=F(\mathrm{haze}(T))$ that follows one branch; the two branches and the
absorbing-haze multiple-equilibria mechanism are the headline Step-3 science
result (Titan's haze feedback is strong enough to be bistable in this model).

**Polydisperse caveat (see `docs/polydisperse_scheme.md`).** The wide bistability
above is for the **monodisperse** closure. Repeating the fixed-haze test with the
bimodal (polydisperse) haze gives a converged warm−cool split of only ~8 K
(stratopause; 17 K monodisperse) at the transition-state haze and ~3–4 K at the
nominal haze (`scripts/bistable_states.py`, paper Fig. 6): the **frozen-haze**
multiplicity is **strongly suppressed but not eliminated** — largely a
single-size artifact.

## RESOLVED: the coupled system is monostable (continuation solve)

The deferred branch-tracking solve is done (`scripts/continuation_solve.py`,
`scripts/plot_continuation.py`, paper sec:coupling-cont +
`writing/figs/continuation_branches.png`). Damped iteration
`T_{k+1} = 0.7 T_k + 0.3 F(haze(T_k))` with the haze recomputed every pass and
the inner RT **initialized from the current iterate** (branch-continuous), fp
tolerance max|F(T)−T| < 1 K, DISORT, 100 layers. Six chains
({mono, bi σ_F=1.2, bi σ_F=2.0} × {warm, cool}) ALL converge to the **same**
coupled fixed point: stratopause 141.7–143.1 K, warm−cool ≤1 K per closure —
**monostable even for the monodisperse haze**. Mechanism: the haze a warm
stratosphere produces admits no warm RT equilibrium (F of the warm start lands
at ~143 K even warm-initialized) — the haze response is net stabilizing and
removes the warm branch; the frozen-haze bistability does not survive coupling.
The naive loop's oscillation = overshoot of a steep (|slope|≫1) single-valued
map, not two coupled attractors.

Caveats: (i) the radiatively-stiff 50–200 km band keeps 10–18 K of
initial-condition memory at tol 0.3 K/day; a tightened-tol (0.05) relax drifts
it monotonically toward the cool-tracked state with no barrier
(`scripts/drift_test.py`: mono 90-km gap 7.5→0.9 K) — convergence memory, not a
second attractor. (ii) The fixed point (~142 K) is well below the observed-like
~183 K — the Step-2 haze-vertical-distribution gap, not the solve.

## Risks / decisions to make

1. **Unconverged Fortran top — RESOLVED.** The haze source sits at ~1 Pa, where
   the Fortran used to oscillate by ~24 K. Fixed via option (a): a **per-layer
   step cap** (`dT_cap=1.0 K` in `run_planetary_radiation.F90`) bounds the top
   oscillation to ~2.5 K @1 Pa with the converged profile unchanged (see
   `docs/rt_discrepancies.md`). The coupled feedback at the haze source is now on
   stable ground; the residual ~2.5 K is averaged out by the last-20-snapshot
   mean.
2. **Cost.** Each outer iteration is a full Fortran run (num_i steps). Warm-start
   from the previous T and cut `num_i` once warm; expect a handful of outer
   iterations.
3. **Fixed-point stability.** Under-relax `T` (and/or the haze τ) between
   iterations; the τ ∝ 1/(1−ω₀) amplification can overshoot.
4. **z↔P consistency.** Use the Fortran's hydrostatic z (or `atm`'s) consistently
   when mapping microphysics(z) onto the `pl` grid.

## Effort

- Writer + round-trip test: small (reuses existing optics).
- Orchestration + warm-start: small–moderate.
- Fixed-point stability + the unconverged top: the real work.
