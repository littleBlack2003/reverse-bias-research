#!/bin/sh
set -eu
export OPENBLAS_NUM_THREADS=1
export MPLCONFIGDIR=/tmp/mpl_gaussian_transport
cd "$(dirname "$0")"
python code/test_gaussian_transport.py
python independent_review/review.py
python code/benchmark_transport.py
python code/audit_saved_states.py
