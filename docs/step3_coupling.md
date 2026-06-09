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

1. **Round-trip**: feed the *observational* haze through `write_presc_haze`
   (bypass microphysics, write the observational τ/ssa/g) → Fortran → confirm the
   same `T(z)` as the stock `preschaze` run. Proves the writer/format is correct.
2. **One-shot**: microphysics haze on the `titan_reference` T → one Fortran run →
   compare to the DISORT coupled run (~140 K stratopause expected).
3. **Iterate**: run the loop to a fixed point; report the feedback (dT vs
   d(opacity), gain, stability).

## Risks / decisions to make

1. **Unconverged Fortran top.** The haze source sits at ~1 Pa (Gaussian at
   ~400 km), exactly where the Fortran oscillates (≳ tens of K, P ≲ 4 Pa). The
   coupled feedback there will inherit that noise. Options: (a) fix the Fortran's
   explicit top stepper / sub-step the thin top layers; (b) cap the production
   region / freeze the top haze; (c) couple only below ~4 Pa and prescribe above.
   **This is the main blocker and should be decided first.**
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
