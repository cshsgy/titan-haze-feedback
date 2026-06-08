# example_bowen_fort — TAM-derived planetary radiation model (reference)

A standalone Fortran radiative-transfer model (J. M. Lora, based on the TAM
radiation code; some routines derived from the legacy Ames Mars GCM), uploaded
to the repo as a **reference / comparison** implementation for our Step 1 DISORT
energy balance. It is the same lineage as Lombardo & Lora (2023).

- `planetary_radiation.F90` — RT routines (`setspv/setspi`, `get_taukcoeff`,
  `get_tauCIA`, `optc`, `sfluxv/sfluxi`): correlated-k gas opacity, CIA, haze,
  two-stream-style solar/IR fluxes.
- `planetary_radiation_driver.F90` — sets up bands, reads k-coefficients / CIA /
  haze, drives the column RT, applies convection, writes diagnostics.
- `run_planetary_radiation.F90` — `PROGRAM`: builds a pressure grid (`nlay=100`),
  iterates to (radiative-convective) equilibrium, writes profiles.

## Build

Needs a Fortran compiler + netCDF-Fortran. We install both in `../../.fortenv`:

```bash
bash ../../scripts/build_fortran.sh     # creates .fortenv (conda) and runs make
# or, with your own toolchain:
make FC=gfortran NFCONFIG=nf-config
```

## Run — required input data (NOT in this repo)

The program reads an `INPUT/DATA/` directory and a few run-dir files that are
**physical data tables and were not uploaded**; it cannot run without them:

| File(s) | Contents |
|---|---|
| `namelist` | run + radiation namelists (`run_planetary_radiation_nml`, `planetary_radiation_nml`) |
| `initial_temperature.txt` | starting T profile (`nlay` values) |
| `INPUT/DATA/solar_spectrum_houghton.txt` | solar spectrum |
| `INPUT/DATA/Rayleigh.txt` | Rayleigh optical depth |
| `INPUT/DATA/ckc_<gas>vis.txt`, `ckc_<gas>ir.txt` | correlated-k coefficients per gas (CH₄, …) |
| `INPUT/DATA/trans_<gas>.txt` | transmission fits |
| `INPUT/DATA/profile_<gas>.txt` | gas VMR profiles |
| CIA + haze files (`haze.nc` or `*_hazev/i.txt`) | aerosol optical properties |

**To run it, provide the `INPUT/DATA/` directory + `namelist` + initial
temperatures** (e.g. from the TAM distribution this was extracted from). Place
them in a run directory alongside the `run_planetary_radiation` binary.

## Outputs (used by the comparison plot)

On a successful run it writes `temperatures.txt`, `sw.txt`/`lw.txt` (heating
rates), `swspec.txt`/`lwspec.txt`, `tsurf.txt`, `surface_fluxes.txt`, …
`scripts/compare_fortran.py` overlays `temperatures.txt` and the SW/LW heating
rates on our DISORT energy-balance figure.
