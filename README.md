# titan-haze-feedback

A coupled radiative–microphysical model of Titan's organic haze. The goal is to
close the loop between the **radiative energy balance** (which sets the
temperature profile) and the **aerosol microphysics** (coagulation +
sedimentation of fractal aggregates, which sets the haze density/size profile,
which in turn sets the opacity that drives the radiation).

## Scientific motivation

On Titan the haze is not a passive tracer: it absorbs sunlight (heating the
stratosphere), emits in the thermal IR (cooling), and its vertical/size
structure is controlled by microphysics that is itself temperature dependent
(gas viscosity η(T), mean free path λ(T,P), thermal velocity, Brownian kernel
∝ T/η). This is a genuine **feedback loop**, and most existing models either
prescribe the haze (Tomasko 2008, Lombardo 2023) or run full 3D microphysics in
a GCM (de Trenquelléon 2025). This project isolates the 1D radiative ↔
microphysics feedback with a fast, analytic **scaling law** for the
microphysics so the loop can be iterated cheaply to a self-consistent steady
state.

## Project plan

### Step 1 — Radiative-transfer energy balance (DISORT)
Build a 1D plane-parallel column for Titan and solve the shortwave + longwave
radiative transfer with **DISORT** (discrete ordinates, multi-stream).
- Shortwave: solar absorption by CH₄ + haze; Lombardo splits at 5 µm
  (2000–40000 cm⁻¹). Haze scattering with single-scattering albedo ω₀ and
  asymmetry g.
- Longwave: thermal cooling, gas line opacity + N₂–N₂/N₂–CH₄/CH₄–CH₄/N₂–H₂ CIA;
  haze treated as a (near-)pure absorber in the IR.
- Gas opacity via correlated-k (HITRAN) or DISR/Irwin CH₄ coefficients.
- Iterate heating = cooling to **radiative(-convective) equilibrium** → T(z).
- **Validation target:** Tomasko (2008) net-flux and heating/cooling tables
  (solar net flux 0.574 W/m² at surface, ~4.1 W/m² at 400 km; net heating peak
  ~0.5 K/Titan-day near 120 km; cooling peak ~21 K/Titan-day near 340 km).

### Step 2 — Microphysics scaling law
Derive a closed analytic relation for the aerosol population — coagulation of
fractal aggregates + gravitational sedimentation → steady-state density/size
profile. Target form:

> **C(z,N) dz dN = f(z, N, T, d, D)**

where `C(z,N)` is the concentration of aggregates containing `N` monomers at
altitude `z`, `T` = temperature, `d` = monomer radius, `D` = fractal dimension.
The full derivation (settling velocity, Brownian kernel, Smoluchowski balance,
all in (z,N) coordinates with explicit T/d/D dependence) lives in
[`docs/scaling_law.md`](docs/scaling_law.md). The physical building blocks are
taken from Burgalat & Rannou (2017) and de Trenquelléon (2025).

### Step 3 — Coupled radiative–microphysical iteration
Close the loop:
1. Start from a guess T(z).
2. Microphysics scaling law (Step 2) → C(z,N) → mass/number/size profiles.
3. Map (M0, M3 → characteristic radius per layer) to optical properties
   (Qext, ω₀, g) via a precomputed aggregate optics table → layer optical depths.
4. DISORT (Step 1) → new heating/cooling → new T(z).
5. Repeat to a fixed point. Diagnose the **feedback** sign/strength
   (e.g. dT/d(opacity), gain factor) and stability.

### Step 4 — Photochemistry in the loop
Add photochemical haze production (currently the source term Q(z) is imposed as
a Gaussian at ~1 Pa / 400 km). Couple production rate / monomer properties to a
photochemistry module (CH₄/N₂ photolysis) so the haze mass flux responds to the
radiation field and composition, completing the radiation ↔ chemistry ↔
microphysics feedback.

## Repository layout (planned)

```
papers/                 source PDFs (reference literature)
docs/
  physics_parameters.md baseline constants & profiles extracted from the papers
  scaling_law.md        derivation of C(z,N) dz dN = f(z,N,T,d,D)
  rt_discrepancies.md   ranked differences between our DISORT RT and the Fortran RT
src/microphysics/       Step 2 scaling-law solver
  constants.py          physical constants + AerosolParams
  atmosphere.py         background T(z), P(z), g(z), eta(T), lambda(T,P), K(z)
  transport.py          fractal geometry, settling omega, coagulation kernel beta
  scaling_law.py        K->0 master-ODE integrator -> n(z), Nbar(z), rho_h(z)
  bvp.py                full eddy-diffusion BVP (master ODE as initial guess)
src/rt/                 Step 1 DISORT energy balance (needs pydisort/.rtenv)
  column.py             layered column from an Atmosphere; heating-rate helper
  cia.py                HITRAN collision-induced absorption (N2-N2, N2-CH4, …)
  optics.py             per-layer (tau, ssa, g): haze (Step 2) + CIA + gray SW gas
  disort_driver.py      pydisort wrapper: shortwave beam + multiband longwave Planck
  energy_balance.py     SW heating, LW cooling, radiative-equilibrium relaxation
data/cia/               band-averaged CIA table (raw HITRAN files gitignored)
scripts/fetch_cia.py    download HITRAN CIA + (re)build the band table
scripts/run_scaling_law.py   demo: integrate + plot the haze profiles
scripts/cross_validate.py    validate against Tomasko/dT25 constraints
scripts/run_energy_balance.py  demo: DISORT heating/cooling + equilibrium T(z)
tests/test_scaling_law.py    master-ODE sanity checks
tests/test_bvp.py            BVP sanity + cross-validation checks
tests/test_rt.py             DISORT energy-balance structural checks (.rtenv)
writing/                paper-style LaTeX writeup + built PDF
```

### Environments

Two local envs (both gitignored, see `.gitignore`):
- **Microphysics (Step 2)** — runs on the system Python (numpy/scipy/matplotlib).
- **RT (Step 1)** — needs `pydisort` (torch-based); installed in `.rtenv`. Run RT
  scripts/tests with `.rtenv/bin/python`.

### Running Step 2

```bash
python3 tests/test_scaling_law.py tests/test_bvp.py   # 12 checks
python3 scripts/run_scaling_law.py     # K->0 profiles -> writing/figs/
python3 scripts/cross_validate.py      # BVP + validation -> writing/figs/
```

### Running Step 1

```bash
.rtenv/bin/python tests/test_rt.py            # structural checks
.rtenv/bin/python scripts/run_energy_balance.py  # heating/cooling + T(z)
```

### Reference Fortran model (`src/example_bowen_fort/`)

A TAM-derived (Lora/Lombardo lineage) Fortran radiation model is included for
cross-comparison. It builds with stock gfortran and runs on prescribed data:

```bash
bash scripts/build_fortran.sh && bash scripts/run_fortran.sh
.rtenv/bin/python scripts/compare_fortran.py   # overlay -> writing/figs/fortran_comparison.png
```

Our DISORT model and the full TAM model agree closely on T(z) (stratopause
~185–190 K, tropopause ~70–75 K, surface ~90 K). Two of its modules
(`haze_mod`, `read_clim_mod`) were missing from the upload and are **reconstructed**
in `reconstructed_stubs.F90` (notably an approximate `saturate`); see that dir's
README for provenance.

**Cross-validation (vs. published Titan constraints).** The eddy-diffusion BVP
reproduces the haze extinction scale height of Tomasko et al. (2008),
**H = 64 km vs. 65 km observed** (the K→0 settling limit alone gives 54 km), and
a main-haze characteristic radius of **~0.3–0.5 µm**, matching de Trenquelléon
et al. (2025) (r_c ≈ 0.46 µm). Monomer mass flux is conserved to machine
precision and the BVP agrees with the master ODE to ~9% in the lower haze.

## Key references

| Paper | Role in this project |
|---|---|
| **Tomasko et al. (2008)**, *Planet. Space Sci.* 56, 648 — *Heat balance in Titan's atmosphere* | Energy-balance targets, gas composition, T(z) anchoring, haze opacity structure, validation tables (Step 1). |
| **Lombardo & Lora (2023)**, *Icarus* 390, 115291 — *Seasonal radiative … Titan* | RT methodology (bands, correlated-k, CIA pairs, haze optics inputs, layering) for the RT solver (Step 1). |
| **Burgalat & Rannou (2017)**, *J. Aerosol Sci.* 105, 151 — *Brownian coagulation of … fractal aerosols* | **Core of the scaling law**: Smoluchowski moment equations, Brownian kernels (continuum/free-molecular/transition), fractal radius laws, settling velocity (Step 2). |
| **de Batz de Trenquelléon et al. (2025)**, *PSJ* 6, 79 — *Titan PCM II: Haze & cloud cycles* | Integrated coupling recipe (moments → optics → RT → T → microphysics) and baseline parameter values (Steps 2–4). |

> Note: the companion papers **Tomasko 2008a** (haze optical constants ω₀, g,
> Qₑₓₜ; fractal-aggregate geometry) and **Tomasko 2008b** (CH₄ k-coefficients)
> are referenced but not in `papers/` — they should be obtained to fully
> populate the optics table and CH₄ opacity in Step 1/Step 3.

## Status

- [x] Literature extracted → `docs/physics_parameters.md`, `docs/scaling_law.md`
- [x] Step 1 — DISORT energy balance (`src/rt/`) via pydisort, **physics complete
      and validated**: **correlated-k CH₄ shortwave** + **spectral haze ω₀(λ)/g(λ)**
      + **mean-field (RDG) aggregate haze cross-section** (`aggregate_optics.py`,
      Khare tholin index; absorption ∝ N monomers, not r_a²) + **longwave = HITRAN
      CIA + correlated-k gas lines (CH₄/C₂H₂/C₂H₆/C₂H₄/HCN)** + convective adjustment.
      **Validation:** on the *prescribed* (observational) haze the engine reaches a
      ~200 K stratopause, matching the reference Fortran (`example_bowen_fort`, ~195 K).
      The RDG cross-section corrected the haze optical-depth magnitude (visible column
      τ ~8, matching Doose, vs ~70 with the old gray mobility-radius cross-section) and
      warmed the coupled stratopause ~115→145 K. The residual ~50 K gap is the haze
      **vertical distribution** (set by the Step 2 production/settling, ρ_h ∝ 1/ω),
      not the radiation. Haze model selectable via `OpticsParams.haze_mode`
      ('rdg' default / 'gray' / 'prescribed').
- [x] Step 2 — scaling-law implementation (`src/microphysics/`): K→0 master ODE
      **and** full eddy-diffusion BVP, cross-validated against Tomasko/dT25
- [ ] Step 3 — coupled iteration
- [ ] Step 4 — photochemistry coupling
