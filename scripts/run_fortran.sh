#!/usr/bin/env bash
# Build (if needed) and run the reference Fortran model, producing
# temperatures.txt / sw.txt / lw.txt in src/example_bowen_fort.
#
# NOTE: the program both READS and OVERWRITES initial_temperatures.txt (it is the
# restart checkpoint).  We stage the pristine seed restart before each run so a
# previous diverged/!= run can't poison the next one.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
D="$ROOT/src/example_bowen_fort"
cd "$D"

[ -x run_planetary_radiation ] || make FC=/usr/bin/gfortran

if [ ! -f initial_temperatures_seed.txt ]; then
  echo "ERROR: initial_temperatures_seed.txt (pristine restart) missing." >&2
  exit 1
fi
cp initial_temperatures_seed.txt initial_temperatures.txt
./run_planetary_radiation
echo "Done. Outputs: $D/{temperatures,sw,lw,tsurf}.txt"
