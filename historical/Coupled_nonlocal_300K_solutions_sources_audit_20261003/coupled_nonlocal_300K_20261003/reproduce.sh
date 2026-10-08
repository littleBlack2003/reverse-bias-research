#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
export OPENBLAS_NUM_THREADS=1
for n in 81 161 321; do
 for d in 0 1; do python3 code/coupled.py "$n" "$d" > "logs/d${d}_N${n}.log" 2>&1; done
done
python3 code/summarize.py
python3 independent_audit/audit.py
python3 -m unittest discover -s code -p 'test_*.py' -v
