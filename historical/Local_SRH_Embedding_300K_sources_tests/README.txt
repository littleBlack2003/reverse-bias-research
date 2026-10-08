Local reversible SRH/DA interface, 2026-10-03
Scope: conditional local closure + UNIT TESTS only. No device scan, quantum kernel, fit, or claimed material parameters.

Files:
Local_SRH_Embedding_300K_20261003.pdf : Chinese scientific specification
SPEC_ZH.txt : editable specification source
local_srh.py : dependency-free reference equations and bookkeeping
test_local_srh.py : 15 unittest tests
tests.log : execution record
build_artifacts.py : figure and PDF builder (numpy, matplotlib, reportlab)
interface_checks.png : illustrative multiplier sweep, depletion counterexample, DA states
results.json : computed test-fixture values

Reproduce from this directory:
python3 -m unittest -v test_local_srh
MPLCONFIGDIR=/tmp/embedding_mpl python3 build_artifacts.py

All upstream reports are read-only and unmodified. Figures vary dimensionless test multipliers, not material field parameters. LocalSRH assumes positive finite coefficients and densities within ordinary double-precision range; does not claim general overflow protection. steady_sources uses analytic U to avoid cancellation of gross rates. DA helper implements only a reversible equilibrium test construction and stoichiometric accounting, not material DA rates or a driven device solver.

Established: old graph uses spatially separated ideal electron reservoirs at -6 and +6 nm. Keep it as a 12 nm transport benchmark. The local valence/conduction SRH candidate is a new model, not a reinterpretation of old R. STOP: identify the NEW model microscopic orbitals, reversible reset processes and controlled closure before device coupling.
