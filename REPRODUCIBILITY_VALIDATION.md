# Release validation, 11 November 2026

Environment: MATLAB R2026b, 26.2.0.3386108, 64-bit Windows. Dependency
analysis of the paper drivers identified 31 working-project source files
and base MATLAB only. Four release entry/export/audit/example functions
were added. Working-project solvers and recorded inputs were copied
unchanged; their identities are in SOURCE_DATA_MANIFEST.csv.

## Executed checks

* `reproduce_simax("verify")`: passed on the recorded 24 dense cases,
  16 solvability records, 12 single-state rows, 12 multistate rows,
  4 evaluator rows, two error histories and four floor records.
  Largest dense projection-error/bound ratio: approximately 0.9426.
* `reproduce_simax("report")`: generated all five numerical tables as CSV
  and LaTeX, and both figures as PDF and PNG. Both PNGs were visually
  inspected for labels, curves, legends and export completeness.
* `reproduce_simax("smoke")`: solved independently generated n=80
  single-state and multistate examples. Accepted dimensions were 21 and 27;
  final absolute errors were about 5.73e-8 and 1.51e-9, below their
  estimated bounds (3.22e-7 and 1.77e-7, respectively).
* `reproduce_simax("dense")`: recomputed all 24 n=600 dense bound cases,
  34 stationary-bound records, 16 solvability records (688 candidate
  stationary problems in total), and the two Figure 1 histories.
* `reproduce_simax("evaluators")`: recomputed all four Table 5 evaluators
  using the n=4000 reference extracted from the actual transport ZIP.
* All non-timing fields of the six recomputed CSVs match the archived
  CSVs exactly as decimal strings. Both error histories match at all
  181 and 124 nodes. See validation/numerical_comparison.json, which
  records the comparison for all six files, and validation/dense_rerun.log. Timings were excluded from this comparison.
* All four reference ZIP archives passed ZIP CRC checks. SHA-256 digests
  of each of their eight MAT members were checked against the frozen
  source-data manifest. The code ZIP is also CRC-checked at packaging.

## Scope

Large single-state and multistate benchmarks and the floor experiment
were archived from the completed paper runs. Full `all` and `cached`
numerical reruns were not executed again during release preparation.
Their original completed CSVs and driver logs are included, with all
coarse reference data and code to regenerate the finer-reference checks.

Release preparation does not constitute a journal decision on the badge.

## v1.0.1 update

Documentation and exported table file names were renumbered to the final
manuscript (Tables 1--5), and the empty file validation/evaluators_rerun.log
was removed. Code paths that no experiment uses were removed from the
sources (README, version history); the recorded data are those of v1.0.0.
The cleaned sources were checked as follows (MATLAB R2026b, Windows).

* `checkcode` reports no problem in `nla_core.m` and no error in any file.
* `reproduce_simax("verify")` and `reproduce_simax("report")` passed and
  produced all five tables and both figures.
* `reproduce_simax("smoke")` gave dimensions 21 and 27, errors 5.73e-8 and
  1.51e-9, and bounds 3.22e-7 and 1.77e-7, as in v1.0.0.
* `reproduce_simax("dense")` and `reproduce_simax("evaluators")` (with the
  accepted n=4000 space recomputed) were rerun. All non-timing fields of the
  six recomputed CSVs equal the archived CSVs as decimal strings; see
  validation/numerical_comparison_v1.0.1.json.
* The large drivers were compared with the v1.0.0 sources on the same
  problems: `run_multistate_transport("run",400)` gives identical rows
  for both multistate models, and for the single-state model with n=4000
  the EBA--BDF1 baseline (basis, final state, trajectory) and Algorithm 5.1
  (bases, status, estimate trace) are bitwise identical.
