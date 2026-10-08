Nonlocal HOMO–trap–LUMO finite-volume reaction adapter, 300 K
2026-10-03

Run: python -m unittest discover -v
Requires Python 3, NumPy, SciPy. Tested with Python 3.12.
41 tests: 18 adapter tests plus unchanged 23 rate tests.

Files
adapter.py: conservative geometry, spectral/fixed kernels, assembly, sparse Jacobian
rates.py: unchanged reversible rate implementation
local_srh_control.py: unchanged prior local SRH implementation for regression only
test_adapter.py / test_rates.py / tests.log: reproducible validation
REPORT_ZH.txt: equations, interface, limitations, measured refinement errors
independent_review.txt: separate numerical review, when included

This is a tested reaction adapter ready to assemble into a new coupled finite-volume
solver. It is NOT an executed self-consistent Poisson/transport solution, a J–V
scan, an onset prediction, or an identified material TAT model. The prior solver
and baseline are untouched. No fitted material rates or parameter scan.

Mesh.faces are control-volume faces in cm; Mesh.nodes hold physical sample
positions strictly inside each volume. Nonuniform and geometric meshes work.
Endpoints are physical cm coordinates; no one-cell-distance replacement is made.
Endpoints outside the first/last resolved sample are rejected, not clipped.
Areal trap counts are Nt*dX in cm^-2. An interpolated point contributes its
shape weight divided by each cell width; do not multiply volume density twice.

State: concatenate phi[N], mu_H[N], mu_L[N], f[M]. Potentials in V,
electron chemical potentials in eV, f dimensionless. mu_H is E_Fp.
Transient outputs: Sn, Sp [cm^-3 s^-1], rho [C cm^-3], df [s^-1],
Jtransfer [A cm^-2] on N+1 faces, drho [C cm^-3 s^-1].
The Jacobian rows concatenate Sn, Sp, rho, df; columns follow the state.
condensed_steady(fields) returns stationary results and a 3N-by-3N reaction
Jacobian after eliminating occupancy. This reaction matrix alone is not a
complete device Jacobian; transport/Poisson/contact terms must still be added.

Model boundaries
The spectral kernel is a conditional nonadiabatic localized sequential Marcus
candidate, not uniquely implied by Huang and not a barrier-shortening model.
The sample energies/couplings/lengths in tests are manufactured numerical data.
Fixed physical paths are used; no ell(F) fit. Material DOS normalization,
thermalized carrier statistics, trap charge reference and counts remain to be
identified before physical device predictions. See the prior rate specification.
