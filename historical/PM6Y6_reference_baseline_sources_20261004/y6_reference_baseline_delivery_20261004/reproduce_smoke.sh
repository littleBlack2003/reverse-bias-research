#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
export OPENBLAS_NUM_THREADS=1
export MPLCONFIGDIR="${TMPDIR:-/tmp}/y6_mpl"
python independent_review/test_fd_eos.py
python independent_review/compare_fluxes.py
python independent_review/test_device.py
python code/calibration.py
python code/validate_observables.py
python code/summarize.py
