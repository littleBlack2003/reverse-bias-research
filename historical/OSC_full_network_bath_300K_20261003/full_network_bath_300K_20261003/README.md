# 300 K full-network bath test

One physical point only: explicit 12 nm / 32-state graph, F=1.5 MV/cm, 300 K.
The only physical change is Gaussian -> quantum Ohmic Ec=20 meV slow bath.
Both internal directions, both reservoir integrals, and conditional energy moments are updated.
The numerical refinement changes Fourier tolerance, interpolation spacing, sideband cutoff,
reservoir quadrature order, and stationary arithmetic precision together. Common-mu controls
are validation of equilibrium, not additional field/temperature scans.

Run from any directory:

    OPENBLAS_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 python -B run_comparison.py
    MPLCONFIGDIR=/tmp/bath_mpl python -B build_report.py

Dependencies: Python 3, NumPy, SciPy, mpmath, Pandas, Matplotlib, ReportLab.
Expected run time: several minutes. Outputs replace only this new project's results.
Original copied source SHA256 values are in source_sha256.json; Library archive provenance
is in SOURCE_PROVENANCE.json. No prior deliverable is replaced.

Acceptance criteria: raw detailed-balance log error <1e-6; refinement relative R change
<1e-6; positive probabilities; equal electronic coupling budget; density mass error <1e-7;
node continuity and total heat/chemical-work relative error <1e-15. High precision closure
is algebraic validation and does not imply corresponding accuracy of the physical kernel.
The direct-rate equilibrium current and normalized edge residual are independently retained.
No negative rates or probabilities are clipped; no spectral truncation is renormalized.
The frozen upstream source's results-directory initialization is harmless and its scan entry
point is never invoked. The copied original source files are not edited.
