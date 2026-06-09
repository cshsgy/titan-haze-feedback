# Discrepancies between the two radiative-transfer models

Our DISORT model (`src/rt/`) vs. the reference Fortran TAM-lineage model
(`src/example_bowen_fort/`). Found by (a) reading the Fortran RT routines
(`setspv/setspi`, `get_taukcoeff`, `get_tauCIA`, `optc`, `sfluxv/sfluxi`) and
(b) running both engines on the **same** state — the Fortran's converged
pressure/temperature profile and the prescribed (observational) haze — and
comparing per-layer heating rates and fluxes (`scripts/rt_diagnostics.py`).

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
