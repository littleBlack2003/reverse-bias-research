FINAL AUDIT: full21_*
These are the final comprehensive numerical audit records: 21 datasets, 168 critical/beyond-open-circuit states and 336 carrier/state EOS checks. full21_summary.json is the scope summary. INDEPENDENT_AUDIT.txt explains the checks, numerical bounds and limitations.

PACKAGE QA: portable_smoke_*
This bounded smoke test checks all three 100 K meshes using the bundled input_with_shunt_saved_states tree after the legacy source directory was moved. It verifies that the audit is portable; it is not a replacement for the full21 audit.

INTERIM CHECKPOINTS: lowT100_* and completed_* (if present)
These names identify intermediate checkpoints, which may cover only a subset. Do not use them as the final audit scope. The final package currently retains lowT100_* as the early 100 K precision checkpoint.

SUPPORTING CHECKS
circuit_identity.csv: direct old/new current and power identities at common biases
same_bias_state_identity.csv: exact reuse of all 1,220 source states
core_reference_identity.csv and immutability_summary.json: code, parameters and original-fit preservation
mesh_convergence.csv: primary and refined-grid changes
root_precision_components.csv: root-region constitutive perturbation indicators

audit_saved_states.py reproduces the full21 saved-state audit without a PDE solve. By default it uses input_with_shunt_saved_states inside the package. See INDEPENDENT_AUDIT.txt for the command and interpretation caveats.
