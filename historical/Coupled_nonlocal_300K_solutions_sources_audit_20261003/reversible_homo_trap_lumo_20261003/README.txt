Reversible HOMO-trap-LUMO candidate, 300 K

Read RATE_SPEC_ZH.txt first. This is a conditional physically explicit sequential transfer candidate, not a Huang2026 reproduction or a material fit.

Run: python3 -m unittest -v test_rates
Dependencies: Python 3, NumPy, SciPy. Exact versions are in SOURCE_PROVENANCE.json.

rates.py uses log-space spectral quadrature and a classical nonadiabatic Marcus kernel. test_rates.py contains 23 diagnostic tests; all input numbers are arithmetic fixtures, not material values.

No changes were made to source PDFs or prior solvers. No device J-V, FF, or knee was computed. The source PDFs are not redistributed in this code package.

The code supports a finite real-valued test envelope; it is not compatible as-is with the prior device solver complex-step Jacobian. Spatial source deposition and path-current verification are mandatory before any nonlocal device coupling.
