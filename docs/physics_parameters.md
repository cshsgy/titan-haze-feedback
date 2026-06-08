# Baseline physical parameters & profiles

Extracted from the four reference papers. Values to initialize and validate the
model. Sources: Tomasko 2008 (T08), Lombardo & Lora 2023 (L23), Burgalat &
Rannou 2017 (BR17), de Trenquelléon 2025 (dT25).

## Aerosol / microphysics parameters

| Quantity | Symbol | Value | Source |
|---|---|---|---|
| Fractal dimension | D (D_f) | **2.0** (free param, 2–3; D=3 ⇒ spheres) | dT25 Table 1 (Cabane 1993); BR17 |
| Monomer radius | d (r_m) | **50 nm** | dT25 Table 1 (Tomasko 2008) |
| Production seed radius | r_p | **10 nm** (BR17 uses 20 nm, r_C=2×10⁻⁸ m) | dT25 Table 1 (Rannou 2004) |
| Aerosol material density | ρ_p | **800 kg m⁻³** | dT25 Table 1 (Trainer 2006) |
| Mass production rate | 𝒬_p | **2.1×10⁻¹³ kg m⁻² s⁻¹** (BR17: 3.5×10⁻¹³ kg m⁻³ s⁻¹) | dT25 (Tomasko 2008) |
| Production altitude | z₀ | **415 km ≈ 1 Pa** | dT25 (Rannou 2004) |
| Production Gaussian width | Δz | **20 km** | dT25 |
| Electric charge density | n_e | **15 e⁻ µm⁻¹** (range 15–30; BR17 tests 0/15/30/45) | dT25 (Lavvas 2010) |
| Slip-flow (Cunningham) constant | A | **1.591** (valid Kn∈[0,5]) | BR17, dT25 |
| Distribution char. radius | r_c | **4.582×10⁻⁷ m** | dT25 |
| Distribution normalization | C₀ | **2.122×10⁻³** | dT25 |

**Fractal geometry (BR17 §2.2.1, the load-bearing relations):**
- Mass / bulk-equivalent radius: `r = d · N^(1/3)`  ⇒ `N = (r/d)³`
- Apparent / mobility radius:     `r_a = d · N^(1/D)`
- Conversion (BR17 Eq. 8):        `r_a = E⁻¹ r^(3/D)`, with `E = d^((D−3)/D)`
- Aggregate volume:               `V = (4/3)π d³ N`  (mass ∝ N ∝ r³ = M₃ moment)
- Defining scaling (dT25):        mass within radius 𝓡 scales as `𝓜 ∝ 𝓡^D`

## Atmospheric structure (Titan)

| Quantity | Value | Source |
|---|---|---|
| Surface temperature | **93.5 K** (Huygens; dT25 93.65 K) | T08 |
| Surface pressure | ~1.465×10⁵ Pa (1.47 bar) | L23 SVRS grid |
| Tropopause | ~68.75 K at ~1.4×10⁴ Pa (~40 km region) | dT25 |
| Stratopause / mesopause | ~30 Pa / ~10 Pa | dT25 |
| T(z) source | Vinatier 2007b (120–500 km) merged to HASI (<120 km) | T08 |
| Stratospheric CH₄ | **1.4 %** mole fraction | T08 |
| CO | 47 ppm (uniform) | T08 |
| H₂ | 0.1 % (uniform) | T08 |
| D/H | 1.32×10⁻⁴ | T08 |
| Gravity | g(z) altitude-dependent (~1.35 m s⁻² surface) | dT25 |

Reported simulated stratospheric T: winter polar stratopause >195 K near 10 Pa;
mid-latitudes 175–185 K near 50 Pa; equinox stratopause ~180 K (L23).

## Radiative transfer setup

| Quantity | Value | Source |
|---|---|---|
| RT grid | 0–540 km (T08); 50–55 layers to ~500 km / 10⁻²–10⁻³ Pa top | T08, L23, dT25 |
| Spectral split | shortwave/longwave boundary at **5 µm** | L23 |
| Shortwave band | 2000–40000 cm⁻¹ (<5 µm); RT integ. 300–3000 nm | L23, T08 |
| Longwave band | 1–2000 cm⁻¹ (>5 µm); thermal 10–1400 cm⁻¹ LBL | L23, T08 |
| Gas opacity method | correlated-k from HITRAN (L23); DISR + Irwin 2006 CH₄ k (T08) | L23, T08 |
| Line absorbers | CH₄, C₂H₆, C₂H₂, C₂H₄, C₃H₄, HCN | L23 |
| CIA pairs (IR only) | N₂–N₂, N₂–CH₄, N₂–H₂, CH₄–CH₄ (negligible >600 cm⁻¹) | L23 |
| CIA data (implemented) | HITRAN, band-averaged on ~33 far-IR bands → `src/rt/cia.py` | Karman 2019 |
| N₂–CH₄ CIA scaling | ×1.50 (below 125 K) | T08 |

**Haze optical structure / coupling:**
- Vertical: extinction ~constant 0–80 km, then **scale height 65 km** above 80 km. | T08, L23
- IR: aerosol treated as **pure absorber** (ω₀ = 0, g = 0). | T08, L23
- Solar optics (ω₀, g, Qext): Doose et al. 2016 (DISR); fractal-aggregate Mie/
  effective-medium table — **needs Tomasko 2008a** for full spectra. | L23, T08
- Second far-IR "Haze A": scale height 50 km, opacity 90–190 cm⁻¹. | T08
- dT25 optics table: 3 particle types × 46 bands × 33 modal radii (10⁻⁹–10⁻⁵ m).

**Per-layer optics from moments (dT25 Eqs. 16–18) — the Step-3 coupling map:**
- Characteristic radius: `r_i = (M₃_i / (M₀_i · α(3)))^(1/3)`
- Layer optical depth:   `Δτ(λ) = Σ_i σ_i(λ, r_i) · M₀_i`
- Mean SSA:              `ω̄ = Σ ω̄_i Δτ_i / Σ Δτ_i`
- Mean asymmetry:        `ḡ = Σ ḡ_i ω̄_i Δτ_i / Σ ω̄_i Δτ_i`

## Energy-balance validation targets (Tomasko 2008, Tables 1–2)

Solar net flux (positive down) / fraction of incident, probe latitude:

| Alt (km) | Net solar (W/m²) | frac | Heating (K/Titan-day) |
|---|---|---|---|
| 404 | 4.129 | 0.808 | 37.4 (@351) |
| 208 | 3.504 | 0.686 | 6.7 |
| 120 | 2.677 | 0.524 | 1.32 |
| 80  | 2.067 | 0.404 | 0.20 |
| 0   | 0.574 | 0.112 | 0.006 |

Thermal net flux (positive up) / cooling rate near Huygens latitude:

| Alt (km) | Net thermal (W/m²) | Cooling (K/Titan-day) |
|---|---|---|
| 489 | 3.010 | 14.1 |
| 342 | 2.972 | 21.5 (peak) |
| 202 | 2.533 | 6.9 |
| 118 | 1.845 | 1.05 |
| 0   | 0.266 | 0.031 |

Energy budget: ~80 % of TOA sunlight absorbed (60 % below 150 km, 40 % below
80 km, ~11 % at surface). Net heating exceeds cooling by ≤0.5 K/Titan-day near
120 km; the residual is exported by circulation. Haze IR emission warms the
troposphere (cuts cooling ~30 % below 50 km). Radiative time constant in the
lower atmosphere exceeds a Saturn year (⇒ a 1D radiative-equilibrium column is a
reasonable leading-order target, with seasonal lag a known correction).

## Gas microphysical properties (needed by the scaling law)

These appear as inputs in BR17/dT25 rather than as hardcoded formulas — supply
Titan-appropriate expressions:
- Dynamic viscosity `η(T)` of N₂ (e.g. Sutherland / power law).
- Mean free path `λ(T,P) = k_B T / (√2 π σ_g² P)` with N₂ collision diameter σ_g.
- Knudsen number `Kn = λ / r_a`.
- Continuum kernel prefactor `K_CO = 2 k_B T / (3η)`.
- Free-molecular kernel prefactor `K_FM = (6 k_B T / ρ_p)^(1/2)`.

## Notes / open items
- DISORT is **not** used in any of these papers (T08 = doubling-adding + LBL;
  L23 = two-stream; dT25 = two-stream source function). This project's DISORT
  multi-stream solver is an upgrade — validate against T08's flux tables.
- The α(k) inter-moment shape function (BR17 Eq. 3 / Burgalat 2014) encodes the
  fixed distribution shape; needed if using a two-moment closure rather than the
  full N-resolved C(z,N).
- Obtain **Tomasko 2008a** (haze optics + aggregate geometry) and **2008b**
  (CH₄ k-coeffs) to complete Step 1/Step 3.
