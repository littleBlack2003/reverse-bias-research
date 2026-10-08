#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 PYTHONUTF8=1 MPLBACKEND=Agg
export MPLCONFIGDIR="$PWD/checkpoints/matplotlib"
mkdir -p data logs figures checkpoints
for mode in PF DD; do
 for nt in 1e14 1e15 1e16; do
  for light in 0 1; do
   python -B code/run_density.py "$mode" "$nt" "$light" 81 .25 30 > "logs/reproduce_${mode}_${nt}_${light}.log" 2>&1
  done
 done
done
for nt in 1e14 1e15 1e16; do
 python -B code/run_density.py PF "$nt" 0 81 .125 30 noforward > "logs/reproduce_halfstep_${nt}.log" 2>&1
 for light in 0 1; do
  python -B code/run_density.py PF "$nt" "$light" 161 .5 30 noforward > "logs/reproduce_mesh_${nt}_${light}.log" 2>&1
  python -B code/refine_windows.py "$nt" "$light" > "logs/reproduce_window_${nt}_${light}.log" 2>&1
 done
done
python -B code/analyze_density.py > logs/analysis.log
python -B code/validate_density.py > logs/validation.log
python -B code/check_invariants.py > logs/assertions.log
python -B code/plot_density.py > logs/figures.log 2>&1
printf 'Complete. Read REPORT.md and data/final_assertions.json\n'
