# Discrepancies between the two radiative-transfer models

Our DISORT model (`src/rt/`) vs. the reference Fortran TAM-lineage model
(`src/example_bowen_fort/`). Found by (a) reading the Fortran RT routines
(`setspv/setspi`, `get_taukcoeff`, `get_tauCIA`, `optc`, `sfluxv/sfluxi`) and
(b) running both engines on the **same** state — the Fortran's converged
pressure/temperature profile and the prescribed (observational) haze — and
comparing per-layer heating rates and fluxes (`scripts/rt_diagnostics.py`).

## Reference top convergence (fixed)

The reference's explicit forward-Euler time-stepper did not reach a steady state
above ~4 Pa: the thin, radiatively- and conductively-fast top layers oscillated
by **up to ~24 K** (snapshot-to-snapshot std), and its global step limiter
(reduce the whole column to 0.1·dt only when *some* layer's step exceeds 80 K)
still let those layers move ~20 K/step. Replacing it with a **per-layer step cap**
(`dT_cap`, run_planetary_radiation.F90 — cap each layer's `|ΔT|` per step
independently) bounds the oscillation without touching the deep layers (small
tendency, never capped) and is exact at steady state. At `dT_cap = 1.0 K` the run
is stable (40/40 finite snapshots) and the top oscillation drops to **~2.5 K @1 Pa
(≤3.6 K for P<10 Pa)**; the converged profile is unchanged (stratopause 195 K,
tropopause 68 K, surface 92 K). Smaller caps (0.5 K) destabilise via the slower
spin-up's coupling to convection/surface, so ~1 K is the practical floor. This
makes the haze-source region (~1 Pa) usable for the Step 3 coupling.

## Same-state quantitative comparison

Fed identical T(p) and haze, the two engines **agree to ~10–20 % through the bulk
of the atmosphere** but **diverge sharply at the very top** (≲ 3 Pa, the
mesosphere):

| pressure | SW ours | SW Fortran | LW ours | LW Fortran |
|---|---|---|---|---|
| 1 Pa   | 59  | 82  | **−102** | **−57** |
| 3 Pa   | 45  | 47  | −25 | −22 |
| 25 Pa  | 26  | 24  | −13 | −15 |
| 200 Pa | 5.0 | 5.0 | −2.7 | −3.0 |
| 1.8 kPa| 0.29| 0.30| −0.16| −0.14 |

(K / Titan-day.) Below ~3 Pa the SW heating and LW cooling track closely; the
big differences are in the topmost layers.

## Ranked discrepancies

**1. Longwave top boundary condition — largest single difference.**
Ours imposes **cold space** (top thermal emission 2.7 K, `temis = 0`): no
downwelling longwave enters from above. The Fortran injects a **warm downwelling
`btop`** ∝ the top-level Planck function (scaled by the padded top optical
depth). Consequently our top layer cools ~2× as fast (−102 vs −57 K/day at 1 Pa).
Ours is the more physically correct boundary for the *model top as the top of the
atmosphere*; the Fortran's warm top partly stands in for the real atmosphere
above its truncation (≳ 415 km). Either way it dominates the upper-atmosphere
disagreement.

**2. Solver: multi-stream DISORT (ours) vs. two-stream (Fortran).**
The Fortran shortwave is **delta-Eddington two-stream** (`gfluxv`); the longwave
is **hemispheric-mean** with a single diffusivity μ̄ = 0.5 (`gfluxi`). Ours is
8-stream discrete ordinates with the full Planck/scattering angular integration.
This is the likely cause of the top SW-heating gap (59 vs 82 K/day) and small
bulk differences, especially where haze scattering matters.

**3. Rayleigh scattering — present in the Fortran, absent in ours.**
The Fortran adds a pressure-scaled, pure-scattering (ω₀=1, g=0) Rayleigh
component (`Rayleigh.txt`, coefficient ∝ ν⁴, steeply rising into the blue/UV).
Our shortwave omits Rayleigh entirely. For Titan the column Rayleigh optical
depth is modest (~0.1–0.2 in the visible), consistent with the close bulk SW
agreement, but it backscatters in the UV/blue and should be added for fidelity.

**4. CIA: HITRAN cross-sections (ours) vs. exp-sum transmission fits (Fortran).**
Ours band-averages HITRAN CIA coefficients and forms τ = Σ k·n_A·n_B·Δz. The
Fortran evaluates a **3-term exponential-sum fit to transmission**
(`trans_*.txt`): T = Σᵢ aᵢ exp(−kᵢ u), τ = −ln T, interpolating *transmission*
linearly in temperature. The exp-sum captures sub-band absorber-amount
nonlinearity (a baked-in mini-k-distribution) that a single per-band cross-section
does not. Different data and different nonlinearity; the close bulk LW agreement
shows the *net* effect is similar, but they are not identical.

**5. Surface longwave boundary.**
Ours: emissivity 1, no reflection (`albedo = 0`). Fortran: emissivity
1 − α_i = 0.95 and the surface **reflects** downwelling longwave (α_i = 0.05).
Small but systematic in the dense lower atmosphere.

**6. Delta-scaling of the haze.**
The Fortran delta-scales the haze in the **shortwave only** (g → g/(1+g),
ω → ω(1−g²)/(1−ωg²), Δτ → Δτ(1−ωg²)) and not in the longwave. Ours passes the
Henyey-Greenstein moments to DISORT directly (no delta-scaling; the 8 streams
resolve the forward peak). Different treatment of the forward-scattering peak.

**7. Surface energy balance.**
The Fortran has an **interactive slab surface** (radiative + bulk sensible/latent
fluxes, time-stepped T_surf). Ours **fixes** the surface temperature. This changes
the lower boundary and the troposphere, though not the radiation per se.

**8. Grid + Planck source.**
The Fortran uses a **doubled half-layer grid** (`L_LEVELS = 2·nlay+3`) with a
linear-in-τ Planck source per half-layer and a precomputed band-integrated Planck
table; ours uses a standard level/layer grid with DISORT's internal Planck. Minor
but affects exactly where sources are evaluated.

## Solar input (verified, then aligned)

The shortwave **solar spectrum is identical**: both read
`INPUT/DATA/solar_spectrum_houghton.txt` and band-integrate it onto the same
`BWNV` edges. Two *scalar* inputs differed and have been aligned in the
diagnostic so the SW comparison is apples-to-apples:

- **Orbit distance.** Fortran `sma = 9.5` AU (namelist), `sol = solarf·rrsun/sma²`;
  ours had `CorrelatedKSW(sma_au=9.58)`. ~1.7 % flux difference (∝ 1/sma²).
- **Insolation cosine.** Fortran is **diurnally averaged** (`diurnal=.false.`) at
  `testing_lat=0`, `testing_ls=0` → declination 0 → `cosz = 1/π ≈ 0.318`; ours
  used `umu0 = 0.35`. ~9 % geometry difference.

With `sma=9.5` and `umu0=1/π` the mid-atmosphere SW heating matches almost
exactly (8 Pa: 36.6 vs 36.6; 25 Pa: 23.9 vs 23.7 K/day). (The coupled model
itself still uses 9.58 AU — Titan's true semi-major axis — and a representative
`umu0=0.35`; the *diagnostic* uses the Fortran's exact values for the one-to-one
comparison.) `albv=0.15`, `albi=0.05` already matched the namelist.

## What is NOT a discrepancy (verified the same)

- **Correlated-k gas overlap:** both sum k at matched Gauss points under the
  perfectly-correlated assumption (no random overlap / resorting).
- **k-coefficient interpolation:** both log-pressure / linear-in-k, linear-in-T.
- **Band grids, per-band solar flux, Gauss weights, gas/CIA pair lists** — shared
  (ours reads the same data files).

## Priority fixes (radiation fidelity)

1. Add **Rayleigh scattering** to the shortwave (cheap, the Fortran data is in
   `INPUT/DATA/Rayleigh.txt`).
2. Decide the **longwave top BC**: keep cold space (physically correct for a
   true top of atmosphere) or add an above-model contribution to mimic the real
   atmosphere above ~415 km.
3. Optionally match the **CIA** treatment (exp-sum fits) and add **surface LW
   reflection** + **delta-scaling** for closer one-to-one agreement.

## Fixes applied (maximizing the match)

The solvable discrepancies above were implemented behind `OpticsParams.match_fortran`
(default on), bringing the two engines into close agreement through the *entire*
column, not just the bulk:

| pressure | SW ours | SW Fortran | LW ours | LW Fortran |
|---|---|---|---|---|
| 1 Pa   | 61  | 82  | **−45** | **−57** |
| 3 Pa   | 46  | 47  | −25 | −22 |
| 8 Pa   | 37  | 37  | −23 | −28 |
| 25 Pa  | 25  | 24  | −13 | −15 |
| 72 Pa  | 13  | 13  | −6.0| −6.6 |
| 200 Pa | 4.6 | 5.0 | −2.7| −3.0 |
| 1.8 kPa| 0.27| 0.30| −0.16| −0.14|

What changed:
- **Longwave top BC (#1).** Replaced cold space with a **warm downwelling top**:
  the truncated atmosphere above the model top emits at the top-level temperature
  with effective emissivity ε = 1 − exp(−τ_above/μ̄), τ_above = `LW_TOP_SCALE`·τ_top·
  (P_top/ΔP_top), μ̄ = 0.5. `LW_TOP_SCALE = 0.02` (tuned) brings the 1 Pa cooling
  from −102 → −45 K/day (Fortran −57). Passed via `P_levels_ascending` to
  `solve_longwave_spectral`.
- **Rayleigh scattering (#3).** Added a pure-scattering (ω₀=1, g=0) molecular
  component to the shortwave optics (`rayleigh_band_tau`, τ_ray = Δp[mbar]·coef(ν),
  `Rayleigh.txt`); the per-wave asymmetry is scattering-weighted so the haze
  forward peak is correctly diluted.
- **Surface LW reflection (#5).** Longwave surface albedo 0.05 (ε = 0.95).
- **Surface SW albedo.** Visible albedo set to 0.15 (`albv`) to match.

Remaining gaps are the **near-massless top layer** (1 Pa: SW 61 vs 82, LW −45 vs
−57), where multi-stream DISORT (#2) and the Fortran two-stream inherently differ,
and the unmatched **CIA exp-sum fit (#4)** and **delta-scaling (#6)** — left as-is
since the net bulk effect already agrees. Set `match_fortran=False` to recover the
physically-pure boundaries (cold space, no Rayleigh, black surface).

## SW per-mass check (closing the +11 K overshoot, phase 3)

Is the model's 2.6x visible-tau excess aloft (vs the observational haze, 100 Pa,
`scripts/diagnose_tau_gap.py`) an optics artifact of the fixed observational
omega0 ~ 0.9 (C_ext = C_abs/(1-omega0))? **No.** For Df=2 the RDG-FA scattering
saturates: C_sca = N^2 C_sca,mono S(kR_g) with S -> (kR_g)^-2 and R_g^2 ∝ N gives
C_sca ∝ N — the same mass scaling as absorption. Per-mass extinction is
size-independent within 5% for N >= 100 (model aloft has Nbar = 1.7e3-3e4,
kR_g = 26-109), and the RDG-FA asymptotic omega0 = 0.944 matches the
observational ~0.92 used. The aloft excess is therefore the **monodisperse
closure's mass distribution**; the bimodal haze already brings tau(P) to
0.96-1.3x of observational below 100 Pa (Phase-0 table). No optics change
warranted. Residual: the 1-10 Pa production zone (model ~1.5-2x obs, tau<0.2).
