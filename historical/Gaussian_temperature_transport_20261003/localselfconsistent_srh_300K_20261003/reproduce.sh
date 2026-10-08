#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/code"
export OPENBLAS_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
python -B integrate.py baseline > ../logs/run_baseline.log 2>&1
python -B -m unittest -v test_local_srh validate.IntegrationTests.test_reference_equivalence validate.IntegrationTests.test_scalar_four_rate_adapter > ../logs/baseline_gate.log 2>&1
python -B integrate.py singlePF > ../logs/run_singlePF.log 2>&1
python -B -m unittest -v test_local_srh validate > ../logs/tests.log 2>&1
python -B forward.py > ../logs/forward.log 2>&1
