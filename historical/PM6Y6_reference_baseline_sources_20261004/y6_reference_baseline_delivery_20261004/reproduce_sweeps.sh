#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
export OPENBLAS_NUM_THREADS=1
mkdir -p logs
# Under223K this runner explicitly selects the unverified mobility extrapolation.
for T in 300 250 200 150 100; do python code/run.py "$T" 321 > "logs/T${T}_N321.log" 2>&1; done
for T in 300 200 100; do
 for I in .25 .08; do ONLY_LIGHT=1 python code/run.py "$T" 321 "$I" blocking 0 0 .05; done
done
