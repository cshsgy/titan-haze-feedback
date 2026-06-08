#!/usr/bin/env bash
# Build the reference Fortran radiation model (src/example_bowen_fort).
# Installs a conda Fortran toolchain + netCDF into .fortenv if missing.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"

if [ ! -x .fortenv/bin/gfortran ]; then
  echo "Installing Fortran toolchain into .fortenv ..."
  conda create -y -p ./.fortenv -c conda-forge gfortran netcdf-fortran make
fi

export PATH="$ROOT/.fortenv/bin:$PATH"
cd src/example_bowen_fort
make clean || true
make FC="$ROOT/.fortenv/bin/gfortran" NFCONFIG="$ROOT/.fortenv/bin/nf-config"
echo "Built: src/example_bowen_fort/run_planetary_radiation"
