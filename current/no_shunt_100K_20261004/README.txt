Current default: no external parallel shunt, 100–300 K.

The only physical output change is Jsh=V/Rsh -> 0; Rs remains0. All material, contact, temperature-extrapolation and field-cap inputs are frozen. No new fit. Main device/FD/field/temperature-law files are byte-identical to the prior model.

Outputs:
figures/PM6Y6_no_shunt_100-300K_LINEAR_JV.png and .pdf
figures/PM6Y6_no_shunt_100-300K_LINEAR_PV.png and .pdf
data/PM6Y6_no_shunt_100-300K_JV.csv: signed current, negative in photovoltaic quadrant
data/PM6Y6_no_shunt_100-300K_PV.csv: main N321 metrics
PM6Y6_去除外部分流_重算说明.txt: full Chinese interpretation and limits
validation/: primary, mesh and independent nodewise precision audit

Code:
python3 code/recompute_no_shunt.py
python3 code/recompute_no_shunt.py --N 641 --metrics-only
python3 code/recompute_no_shunt.py --temperatures 100 125 150 --N 1281 --metrics-only
python3 code/build_artifacts.py

Existing successful results are guarded caches. For a fresh run, first move data/ to a recoverable backup and create a new empty data/ directory. Run N321 before refined meshes. The bundled input_with_shunt_saved_states directory contains historical same-bias PDE states. Their old JSON terminal-current values are provenance only; the new program re-evaluates intrinsic currents with no shunt, requires all state accuracy gates, and solves every uncovered voltage/root/MPP. A copied/interpolated mesh profile is only an initial guess and is never accepted without a self-consistent solve. No previous model directory is needed for recomputation.

N321 curves use 121 evenly spaced bias targets through each new Voc+40mV, plus historical computed points and newly solved critical-point evaluations. Lines in plots connect actual computed states. New Voc is always bracketed by solved points of opposite sign with absolute current larger than its continuity/cancellation bound. Root solver tolerance does not represent mesh or physical-model uncertainty.

Model validity: the low-temperature mobility, beta and gamma laws and idealized contacts are conditional assumptions. Numerical convergence does not establish experimental validity. 100K mesh shifts are non-monotonic at the few-microvolt scale; no Richardson extrapolation or higher-order accuracy claim is made.

Historical room fit is in archive_record/room_fit_with_shunt_original.json. Its targets/residuals are not active claims in this no-shunt version. Current settings are in reference/room_fit_frozen.json. Source inputs and result hashes accompany the archive.
