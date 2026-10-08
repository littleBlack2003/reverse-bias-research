# Backup integrity and security review

Date: 2026-10-08 UTC. Scope: local staging preparation only; no external upload was performed by this preparer.

## Verified

- Every retained native file was extracted/copied without byte changes and compared to its source-member SHA-256.
- All 169 members of `PM6Y6_no_shunt_100-300K_reproducible_20261004.zip` are retained. Its original `SHA256.json` validates all 168 other members.
- All 38 members of `PM6Y6_core_expert_review_20261005.zip` are retained. Six core files and four input files match both the adopted archive and its expert-review provenance.
- All nine expert-review unit tests passed on the installed Python / NumPy 2.3.5 / SciPy 1.17.0 environment. The test log is `verification/backup_core_tests.txt`.
- Preserved NPZ arrays were read with `allow_pickle=False`; no object arrays were found.
- Extraction rejected absolute archive members, parent traversal and symlinks. No executable binary, credential configuration, runtime cache, private assistant note or unrelated personal conversation was selected.
- Text scanning found no recognized access-token/private-key signatures, suspicious secret assignments, signed URLs, or private-note references. This is a bounded static review, not a proof that arbitrary historical code is safe to execute.
- Retained PDFs are generated project reports/figures. Bundled third-party literature PDF/OCR trees were omitted with hashes and citations/manifests retained.
- Current TeX `includegraphics` dependencies exist. Exactly seven current TeX font files remain; their original license notices are preserved.

## Source-preserving portability exceptions

The adopted current archive and portable expert-review package are preserved in full. Older research/rendering/packaging source contains original absolute paths, including monograph `source/native_figures.py`, `source/redraw_frozen_serif.py`, historical density-study `code/package_checkpoint.py`, bath provenance and matched-evidence restoration provenance. These identify the original computation workspace and require manual remapping if the historical helpers are reused. They contain no credentials and were not normalized because changing original source/provenance would weaken restoration evidence.

Legacy font-download helpers reference public official font repositories. They were retained as source and were not executed. Existing scientific source scripts may launch local Python processes or XeLaTeX. None was granted new external access during this backup. Rebuild the preserved current TeX directly rather than using an old content-conversion launcher.

## Not claimed

No fresh full 100–300 K scan, no new material calibration, no candidate-model adoption, no new experimental validation, and no monograph PDF rebuild were performed. Historical omitted outputs are not represented as retained. Passing the nine tests establishes their numerical/source checks only.
