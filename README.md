# Reverse-bias research backup — 2026-10-08

This is a source-preserving backup of the reverse-bias / PM6:Y6 research. No scientific model, calibration value, numerical result, or native source file was changed for this backup.

## Start here

- `current/no_shunt_100K_20261004/` is the complete adopted no-external-shunt model archive: all 169 original members, including code, calibration inputs, numerical states, figures, numerical audits, dependencies, and historical continuation seeds.
- `expert-review/PM6Y6_core_expert_review_20261005/` is the complete portable expert-review package: six unchanged model source files, four core inputs, 100 K / 300 K reference states, nine tests, and equation-to-code mapping.
- `reports/monograph/OSC_reverse_knee_editable_20261002/main.tex` is the current editable monograph in the recovered source archive. The separate report PDF is retained under `reports/`.
- `historical/` preserves earlier mechanism studies, source code, calibration/input tables, and textual reports. It is a curated source history, not a complete copy of every earlier numerical output directory.
- `revisions/research_sync_20261008/` holds later candidate revisions synced from the expert-review working folder (start with `00_研究进展总览_ZH.md`): revision v2 (density-dependent mobility closure, compared with the 2022 batch) and v3 (2020-batch per-temperature fits, smooth-function hold-out tests, mesh checks), plus an early reverse-bias script and an unconverged μ_eff master-equation study. Each revision folder keeps the original runnable zip and an extracted copy. These are NOT the adopted model; v3 ships its own `extend_device.py` with the temperature guard widened to 70–300 K, while the adopted copies are unchanged. This folder is not covered by `inventory.json`.
- `inventory.json` records source archive paths, Library identities when available, archive SHA-256 values, every retained native member, and individually enumerated omissions with their original hashes.
- `MANIFEST.sha256` covers every file in this backup except itself. `BACKUP_REVIEW.md` states what was checked and what was not.

## Scientific status and interpretation

The adopted model has Rs=0 and zero external parallel-shunt conductance; its terminal current is the intrinsic model current. Removing the previous fitted shunt while freezing other parameters changes the old room-temperature fit. Preserved older fitted residuals are historical provenance, not active fit claims.

100–300 K is the saved simulation range. PM6:Y6 literature-derived calibration and measured-observable constraints are retained with their source-specific temperature ranges; this does not establish experimental validation over the entire 100–300 K range. In particular, the low-temperature mobility, field-coefficient and recombination extrapolations, fixed generation and idealized contacts remain conditional assumptions. The saved 100 K result is an unvalidated model prediction. Neither passing numerical tests nor mesh convergence proves that the model explains this user's sample or its reverse-bias knee.

The `optional_candidate_NOT_ADOPTED/` directory is an isolated density-dependent candidate. It is not the adopted model and is not imported by the normal expert-review entry point. Do not merge its formulas or parameters into the adopted baseline without a separate scientific decision.

## Restore and check the current model

Use a working copy if you want to compute new results. The backup's source hashes describe the frozen files.

The source lock files record NumPy 2.3.5, SciPy 1.17.0 and, for plotting, Matplotlib 3.10.8. Python 3.12.14 was the archived environment. Compatible dependencies must already be installed or installed separately; no environment, executable package cache or credentials are bundled.

For a lightweight offline check:

```sh
cd expert-review/PM6Y6_core_expert_review_20261005
python3 -B review.py hashes
python3 -B -W ignore::ResourceWarning -m unittest discover -s tests -v
python3 -B review.py smoke
```

The smoke/test checks re-evaluate saved states. They are not a fresh temperature scan or an independent experimental validation. For a perturbed-seed nonlinear solve, see the expert package README. Its small input set cannot run the full scan driver.

For full-archive recomputation, use `current/no_shunt_100K_20261004/README.txt`. It contains commands for N321, N641 and selected N1281 runs. Preserve the existing `data/` separately and use an empty output `data/` when explicitly intending a fresh run; successful results are otherwise guarded caches. The bundled `input_with_shunt_saved_states/` is immutable continuation provenance, not the current terminal-current convention.

## Editable report

The original source is unchanged. The seven fonts referenced by current `main.tex` are retained, together with the bundled Noto/Liberation/Latin Modern license notices and provenance. Earlier STIX/fallback renderer fonts were omitted because they are not used by current `main.tex`; their original archive hashes are listed in the omission inventory.

For the current TeX document, run XeLaTeX directly on the preserved `main.tex` in its own directory with shell escape disabled, three times for cross-references, using an existing TeX Live/MiKTeX installation:

```sh
xelatex -no-shell-escape -interaction=nonstopmode -halt-on-error main.tex
```

Do not regenerate `main.tex` from historical fragments merely to restore this snapshot. Some older rendering/packaging scripts retain original absolute workspace paths and older layout assumptions; these are preserved as source evidence, not advertised as portable launchers. The current TeX figure references were checked for availability, but the PDF was not rebuilt during this backup.

## Deliberate exclusions and scope

Excluded: third-party literature PDFs and extracted full-text/OCR/page-image collections, unrelated projects or conversation content, caches, build/execution logs, redundant historical rendered figures/PDFs, superseded font/rendering assets, and historical per-state numerical outputs. Literature URLs, citations, manifests and source hashes remain where present. Original Library archives remain separately addressable by their recorded IDs; this backup does not delete or replace them.

All adopted current-model and expert-review members are retained. Historical scripts may require rerunning their computation to recreate omitted states, restoring required old archive members, and remapping original workspace paths. In particular, archived native source files and their own old manifests were not rewritten to pretend the curated historical selection is a byte-complete historical archive. Use this backup's `inventory.json` and `MANIFEST.sha256` for backup coverage.
