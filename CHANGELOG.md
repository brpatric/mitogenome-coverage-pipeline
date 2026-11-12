# Changelog

## 4.2.2 (November 2026)

- `run_pipeline.sh` runs `diagnose_padding_artifact.sh` automatically at the end of every run (Phase 8), writing `diagnostic_report_<specimen>_<date>.txt` to `coverage_run/` and recording the double-count result in `run_provenance.txt`. A problem in this report-only step never stops the run.
- `diagnose_padding_artifact.sh` gains a `--from-driver` mode (report only). Run by hand, it behaves as before and also rebuilds `prelim_filtered_depth.txt`.
- Computed values are unchanged from 4.2.1.

## 4.2.1 (November 2026)

- `run_pipeline.sh` now requires an answer on whether the assembly was rotated: new setting `ROTATED` ("yes"/"no"). If it is empty, the driver asks at the start (and, if rotated, asks for the unrotated file); a run with no terminal stops with a message instead of guessing.
- For a rotated assembly, Phase 0 runs the overlap check on the unrotated assembler output (`overlap_check_original.log`), verifies the rotation (`rotation_check.log`) and checks the rotated file, stopping if any check fails. This replaces the manual pre-rotation overlap check.
- `run_provenance.txt` records `ROTATED` and all three results; the script checksums now include `02_fold_and_plot.py`.
- Computed values are unchanged from 4.2.

## 4.2 (November 2026)

- New end-window check. `run_pipeline.sh` Phase 7 maps the reads, unpadded, to a copy of the assembly starting halfway round the circle and writes `origin_check_depth.txt`; the shift is recorded in `padding_metadata.txt` (`origin_check_shift_S`).
- `02_fold_and_plot.py` compares the folded profile with that profile across the end windows (last and first *K* bp) in 200-bp bins and writes `<prefix>_end_window_check.txt`. A bin fails at >= 10% and >= 3x difference. If any bin fails, the end windows of the reported profile are taken from the origin-check mapping and the folded profile is kept as `<prefix>_folded_depth_before_end_window_correction.txt`.
- The statistics file gains an `End-window check:` line.
- Reason: in two rotated *Hogna carolinensis* assemblies the folded profile overestimated depth near position 1 while the interior comparison passed.
- Profiles that pass the check are unchanged from 4.1.1.

## 4.1.1 (November 2026)

- `run_pipeline.sh`: read files may be `.fastq` or `.fastq.gz` and are used where they are, not copied into `coverage_run/`.
- `ASSEMBLY`, `R1`, `R2` and `ORIGINAL_ASSEMBLY` accept file names, relative paths or absolute paths; the driver stops with a clear message if a file is missing.
- `run_provenance.txt` records full input paths.
- Computed values are unchanged from 4.1.

## 4.1 (November 2026)

- Provenance: `run_pipeline.sh` writes `run_provenance.txt` (date, software versions, settings, MD5 checksums of inputs and scripts, check results, *K*, alignment counts).
- Every script carries a version number; `padding_metadata.txt` records the assembly checksum and pipeline version.
- Statistics and comparison files end with a provenance block (script version, date, library versions, input checksums).
- New optional `ORIGINAL_ASSEMBLY` setting: the driver verifies a rotated assembly before mapping and records the result.
- Computed values are unchanged from 4.0.

## 4.0 (November 2026)

- `run_pipeline.sh` Phase 6 builds the unpadded control depth (`prelim_filtered_depth.txt`) automatically.
- `02_fold_and_plot.py` runs the control comparison when that file is present; `03_compare_to_control.py` kept as an optional standalone script.
- `diagnose_padding_artifact.sh` is now a troubleshooting tool only.
- Folded depth, statistics and figures unchanged from 3.0.

## 3.0 (November 2026)

- `03_compare_to_control.py` added to validate the correction against an unpadded control.
- Overlap check run on the assembler output before rotation, followed by `verify_rotation.py`.
- Reporting guidance split into Methods and Results; Methods template and command-line supplement added.
- Documentation made generic.

## 2.0

- Circularity correction by padding and folding, terminal-overlap check, character validation, MAPQ guidance, `verify_rotation.py` and `diagnose_padding_artifact.sh`.
