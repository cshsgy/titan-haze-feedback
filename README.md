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

### Step 1 — Radiative-transfer energy balance
Solve the shortwave + longwave radiative transfer for a 1D plane-parallel Titan
column. **The Step 1 engine is the reference TAM-lineage Fortran model**
(`src/example_bowen_fort/`, correlated-$k$ two-stream; see below); our
multi-stream **DISORT** implementation (`src/rt/`) is a validated cross-check.
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

## Repository layout

```
papers/                 source PDFs (reference literature)
docs/
  physics_parameters.md baseline constants & profiles extracted from the papers
  scaling_law.md        derivation of C(z,N) dz dN = f(z,N,T,d,D)
  rt_discrepancies.md   ranked differences between our DISORT RT and the Fortran RT
  step3_coupling.md     Step 3 wiring plan + bistability result/diagnosis
  polydisperse_scheme.md bimodal 2-moment log-normal scheme (B&R 2017) + results
src/microphysics/       Step 2 scaling-law solver
  constants.py          physical constants + AerosolParams
  atmosphere.py         background T(z), P(z), g(z), eta(T), lambda(T,P), K(z)
  transport.py          fractal geometry, settling omega, coagulation kernel beta
  scaling_law.py        K->0 master-ODE integrator -> n(z), Nbar(z), rho_h(z)
  bvp.py                full eddy-diffusion BVP (master ODE as initial guess)
  moments.py            polydisperse: log-normal alpha(k), moment-avg settling
  coagulation.py        polydisperse: bimodal Fuchs-kernel coag (GH quadrature)
  scaling_law_bimodal.py polydisperse: K->0 bimodal (S+F) master ODE + to_micro
  bvp_bimodal.py        polydisperse 8-field BVP (relaxation implemented but
                        impractically stiff; needs operator splitting)
src/coupling/           Step 3: presc_haze.py writes Fortran prescribed-haze tables
src/example_bowen_fort/ Step 1 RT engine: reference TAM-lineage Fortran model
src/rt/                 DISORT cross-check of Step 1 (needs pydisort/.rtenv)
  column.py             layered column from an Atmosphere; heating-rate helper
  cia.py                CIA: HITRAN band-avg (CIABands) + exp-sum fits (CIAExpSum)
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

### Step 1 RT engine — the reference Fortran model (`src/example_bowen_fort/`)

**Decision (Jun 2026):** the **TAM-derived (Lora/Lombardo lineage) Fortran
radiation model is the Step 1 RT engine.** It is a mature, complete
correlated-$k$ model (CH₄ shortwave; CH₄/C₂H₂/C₂H₆/C₂H₄/HCN gas lines;
N₂–N₂/N₂–CH₄/CH₄–CH₄/N₂–H₂ CIA via exp-sum transmission fits; prescribed haze;
delta-Eddington two-stream SW + hemispheric-mean LW) with surface energy balance
and convective adjustment. It builds with stock gfortran and runs on prescribed
data:

```bash
bash scripts/build_fortran.sh && bash scripts/run_fortran.sh
.rtenv/bin/python scripts/compare_fortran.py   # overlay -> writing/figs/fortran_comparison.png
```

Step 3 will drive **this** RT with the Step 2 microphysics haze. Two of its
modules (`haze_mod`, `read_clim_mod`) were missing from the upload and are
**reconstructed** in `reconstructed_stubs.F90` (notably an approximate
`saturate` and the prescribed-haze interpolation); see that dir's README for
provenance.

**Top convergence (fixed).** The explicit stepper used to oscillate by up to
~24 K above ~4 Pa; a **per-layer step cap** (`dT_cap=1.0 K` in
`run_planetary_radiation.F90`) now bounds it to ~2.5 K @1 Pa (≤3.6 K for
P<10 Pa) with the converged profile unchanged — the haze-source region is now
usable for Step 3. The comparison still uses the time-mean of the last 20
snapshots. See `docs/rt_discrepancies.md`.

### DISORT cross-check (`src/rt/`)

Our `pydisort`/DISORT implementation (multi-stream, Python) is retained as an
**independent validation** of the Step 1 engine, not the production RT. It was
matched to the reference source-by-source: with the same insolation
(`sma=9.5` AU, diurnally-averaged `cosz=1/π`), Rayleigh, surface albedos, the
reference's **exp-sum CIA** (`CIAExpSum`, reading the same `trans_*.txt` tables),
and the **same 100-layer log-pressure grid** (top 1 Pa), **every opacity source
agrees to within a few percent** (`scripts/opacity_breakdown.py`); on the
prescribed haze both reach a ~200 K stratopause. The residual same-state
heating-rate differences are the solver (8-stream DISORT vs the reference's
two-stream) and boundary details, not the optics. Run
`scripts/rt_diagnostics.py` and `scripts/opacity_breakdown.py` for the
breakdown.

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
- [x] Step 1 — **engine chosen: the reference TAM-lineage Fortran model**
      (`src/example_bowen_fort/`), a mature correlated-$k$ two-stream RT with
      surface energy balance + convective adjustment. Step 3 drives it with the
      Step 2 haze.
- [x] Step 1 cross-check — **DISORT (`src/rt/`) validated against the engine and
      now matches it source-by-source**: correlated-k CH₄ shortwave, spectral haze
      ω₀(λ)/g(λ), mean-field (RDG) aggregate cross-section, **exp-sum CIA matching
      the reference's `trans_*.txt`** + gas lines, same insolation, same 100-layer
      log-pressure grid. **Every opacity source agrees to within a few percent**;
      on the prescribed haze both reach a ~200 K stratopause. Residual same-state
      heating differences are the solver (8-stream vs two-stream), not the optics.
      (DISORT's coupled microphysics-haze run gives a ~140 K stratopause; the
      ~50 K gap to observed is the haze **vertical distribution** from Step 2, not
      the radiation.)
- [x] Step 2 — scaling-law implementation (`src/microphysics/`): K→0 master ODE
      **and** full eddy-diffusion BVP, cross-validated against Tomasko/dT25
- [~] Step 3 — coupled iteration **built and run** (`scripts/run_coupled.py`;
      `src/coupling/`): a Python loop regenerates the Fortran's prescribed-haze
      tables from the Step 2 microphysics each pass and cold-starts the engine —
      **no Fortran source changes** (`haze_data='presc'`). Round-trip + one-shot
      validated (`tests/test_presc_haze.py`, 39 checks). **Finding:** the feedback
      is **strong and near-bistable** — the loop oscillates between a warm
      (~183 K) and a cool (~131 K) stratopause and does not settle under naive
      iteration. **Diagnosed** (`scripts/diagnose_transition.py`,
      `scripts/rt_multiplicity.py`): the microphysics T→haze map is smooth (no
      cliff), but the RT is **genuinely bistable** for a fixed haze (two
      radiative–convective equilibria, 17 K apart at the stratopause / up to 31 K
      through the profile, engine-independent) — a real absorbing-haze radiative
      feedback (haze rises with T → absorbs higher → warmer), not a microphysics
      or mapping artifact (see `docs/step3_coupling.md`,
      `writing/figs/transition_diagnosis.png`). The earlier
      top-oscillation blocker is resolved (per-layer step cap, ~24→2.5 K).
      **Resolved by the continuation solve** (`scripts/continuation_solve.py`,
      paper §coupling-cont): the **coupled** system is **monostable** — damped
      branch-tracking converges warm- and cool-tracked chains to the same fixed
      point for every closure.
      **Headline (supersedes the K-numbers above): the haze IR opacity sets the
      feedback regime.** The oscillation/bistability/142 K fixed point all
      belong to the *lab-tholin* LW configuration: Khare (1984) k over-absorbs
      the 300–900 cm⁻¹ thermal window 20–90× vs the observational haze
      (`scripts/diagnose_tau_gap.py`, `figs/tau_gap.png`), over-cooling the
      stratopause ~40 K. Calibrating the LW spectral absorptivity per unit
      visible extinction to the observational tables (`lw_haze='obs'`, default)
      gives +53 K one-shot (141.5→194.2 K) and a **benign** coupling: flat
      composite map (slope ~0.1), all six chains → **~194 K stratopause** (mono
      split 0.1 K), frozen-haze splits ≤6 K — in agreement with the
      prescribed-haze references (Fortran engine 195 K, same-engine DISORT
      ~200 K; never compare stratopauses across engines). Lab-tholin
      IR constants manufacture a dramatic feedback that is not there.
- [x] Step 3 follow-up — **polydisperse microphysics**
      (`src/microphysics/{moments,coagulation,scaling_law_bimodal}.py`,
      `docs/polydisperse_scheme.md`): bimodal 2-moment log-normal scheme
      (Burgalat & Rannou 2017), spheres + fractal aggregates with moment-averaged
      settling and full-Fuchs-kernel coagulation. **Findings:** (i) gravitational
      sorting thins the haze steeply with the aggregate width σ_F (visible column
      τ 8.4 monodisperse → 6.2 at σ_F=1.2 → 1.9 at σ_F=2.0; observed τ≈8 needs
      σ_F≲1.2 — a constraint on Titan's aggregate spread;
      `scripts/sigma_f_sweep.py` — LW-independent, stands). (ii) In the
      *lab-tholin* configuration the frozen-haze bistability is strongly
      suppressed by polydispersity (~8 K vs 17 K mono;
      `scripts/{bistable_states,plot_bistable,check_bimodal_converge}.py`); with
      the obs-calibrated LW there is no robust multiplicity for either closure.
      The 8-field eddy-diffusion BVP: a MOL+banded-LSODA pseudo-time relaxation
      is implemented (`solve_bimodal_relax`) but impractically stiff even at
      60 nodes; operator splitting (implicit transport + sub-cycled coag) is
      the identified path. Bimodal results use K→0 profiles (eddy correction
      bounded by the mono case, H 54→64 km).
- [ ] Step 4 — photochemistry coupling
