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
* `reproduce_simax("evaluators")`: recomputed all four Table 6 evaluators
  using the n=4000 reference extracted from the actual transport ZIP.
* All non-timing fields of the six recomputed CSVs match the archived
  CSVs exactly as decimal strings. Both error histories match at all
  181 and 124 nodes. See validation/numerical_comparison.json and the
  sanitized rerun logs. Timings were excluded from this comparison.
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
