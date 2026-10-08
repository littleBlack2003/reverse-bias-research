300 K micrograph to device current audit

Run python -B audit_bridge.py then python -B build_report.py.
Dependencies: numpy pandas matplotlib reportlab. No package installation or upstream physical solver execution is needed.
inputs/ contains unchanged cited source excerpts; SOURCE_PROVENANCE.json hashes inputs consumed by the audit. Source files were only read.
results/ includes raw reconstructed directed/gross event fluxes, conditional density/coupling tables, fixed-rate cut capacity upper bounds, existing-field-point extracts, and assertion results.
All density/current conversions are conditional, not experimental fits. 50 mA/cm2 is diagnostic only.
The five-page PDF and image were visually checked.
