Independent audit: PASSED within the tested dark-device envelope

Reproduce from the package root:
  PYTHONDONTWRITEBYTECODE=1 python independent_audit/audit.py

24 states checked: N=81/161/321; per-leg distance=0/1 nm;
V=0/-1/-5/-15 V, T=300 K. Production and archived source were not edited here.

Independent checks:
- Archived conservative reaction-adapter equivalence: source error <1.8e-13,
  charge error <1.1e-13, occupancy error <6.9e-14.
- Independent d=0 closed-form FD Marcus source: error <2.6e-14.
- Full colored Jacobian equals individual-column complex steps exactly;
  independent five-point directional check: normalized error <4.5e-10.
- Poisson trap-charge algebra: discrepancy <7.2e-15.
- No sources deposited on contact nodes; trap count 9.6e9 cm^-2 on all grids.

Global diagnostic maxima: current conservation 1.35e-9, Gauss law 1.06e-10,
electrical/free-energy ledger 6.02e-11, all relative. Exact production-source
hash and per-state results are in results.json; equations and interpretation
are in REPORT.txt.

Key limits: this is a conditional one-level FD/lattice-gas/Marcus model, not a
material fit. Mobility is dilute hopping mobility. Orientations are separate
half-Nt trap populations. One spectral state per endpoint and the +q/2 empty-
trap charge reference are explicit assumptions. Distance is per leg, so the
finite H-to-L span is 2 nm. Its proper matched control is d=0 with the same
closure and trap count; the archived MB-SRH baseline is a separate regression.
The energy ledger is dark and isothermal, not a temperature solve. Arbitrary
extreme-state overflow robustness and material predictive accuracy are untested.
