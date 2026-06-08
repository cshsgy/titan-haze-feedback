# example_bowen_fort — TAM-derived planetary radiation model (reference)

A standalone Fortran radiative-transfer model (J. M. Lora, based on the TAM
radiation code; some routines derived from the legacy Ames Mars GCM), included
as a **reference / comparison** for our Step 1 DISORT energy balance. Same
lineage as Lombardo & Lora (2023).

- `planetary_radiation.F90` — RT routines (`setspv/setspi`, `get_taukcoeff`,
  `get_tauCIA`, `optc`, `sfluxv/sfluxi`): correlated-k gas opacity, CIA, haze,
  two-stream-style solar/IR fluxes.
- `planetary_radiation_driver.F90` — bands, reads k-coefficients / CIA / haze,
  drives the column RT, convection, diagnostics.
- `run_planetary_radiation.F90` — `PROGRAM`: pressure grid (`nlay=100`), iterates
  to radiative-convective equilibrium, writes profiles.
- `reconstructed_stubs.F90` — **reconstructed** `haze_mod` + `read_clim_mod`
  (see below).

## Build & run

```bash
bash ../../scripts/build_fortran.sh   # or: make FC=/usr/bin/gfortran
bash ../../scripts/run_fortran.sh     # stages the seed restart, runs the binary
```

Outputs `temperatures.txt` (row 0 = pressure [Pa], later rows = T snapshots),
`sw.txt`/`lw.txt` (heating rate [K/s]), `tsurf.txt`, …
`../../scripts/compare_fortran.py` overlays the result on our DISORT figure
(`writing/figs/fortran_comparison.png`).

## What had to be reconstructed / fixed (provenance)

The upload + the data archive (`~/archive_20260515_cleanup.zip`) did **not**
include the `haze_mod` / `read_clim_mod` source (only Mac-compiled artifacts) and
the driver had non-portable constructs. To get a clean build+run with stock
gfortran we:

- **`reconstructed_stubs.F90`** — minimal `haze_mod` (`haze_init` reads &
  bilinearly interpolates the prescribed `INPUT/DATA/preschaze*` files; climatology
  routines stub out) and `read_clim_mod` (`saturate` = an approximate
  Clausius-Clapeyron cap; `readclim_MOD` stubs out). Used in the prescribed-data
  configuration (`haze_data='presc'`, `clim_gas=''`) of `namelist`. **These are
  reimplementations, not the original physics** — `saturate` in particular is an
  approximation; drop in the original `haze.F90`/`read_clim.F90` for full fidelity.
- Driver fixes (commented in-source): non-standard array-range constructors
  `(/600:300:-50/)` → explicit literals; `error_mesg` stub signature
  `character,dimension(300)` → `character(len=*)`; `if (1)`/`if (0)` →
  `if (.true.)`/`.false.`; commented the dead `use netcdf` (no `nf90_` calls are
  reached, so netCDF is not linked).

## Data & gotchas

- `INPUT/DATA/` holds the k-coefficients, CIA transmission fits, gas profiles,
  prescribed haze, solar spectrum, Rayleigh (from the archive). The 23 MB
  `haze.nc` is **not used** (prescribed haze) and is gitignored.
- `initial_temperatures.txt` is BOTH the restart input AND the checkpoint the run
  overwrites — a diverged run poisons the next one. `run_fortran.sh` re-stages the
  pristine `initial_temperatures_seed.txt` before every run.
