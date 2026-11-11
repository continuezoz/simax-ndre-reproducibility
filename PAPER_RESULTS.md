# Numerical result provenance

The numerical tables are generated from the full-precision recorded CSVs.
No experiment data are transcribed from the typeset manuscript.

* Table 2: group `wellposed.csv` by model and scaled/unscaled variant;
  report stationary-failure counts in Leja/increasing order, the minimum
  closed-loop real part and maximum stationary norm across both orders.
  The model nu values are taken from the matching dense-reference rows.
* Table 3: retain dense-test tau values 1e-4, 1e-7 and 1e-9 and report
  projection-error columns. The separate computed-output-error columns
  are retained in the full CSV for inspecting rounding effects.
* Table 4: join the two proposed-evaluator rows and the BDF1 baseline on n.
  Newton backward errors are parsed from the recorded driver log (last
  completed occurrence for each n), rounded in that log to four digits.
* Table 5: join the increment, Radon and baseline rows on n and gamma.
  The paper prints the increment error; both evaluators remain in the data.
* Table 6: direct export of evaluator, maximum error, peak node and time.

The logs retain the history of interrupted and resumed reference runs.
The final CSVs contain the completed paper results. Obsolete
`dense_bounds_v1.csv` is not part of this release. The reference MAT files
contain coarse-space bases and projected trajectories. The refinement
comparisons are made before saving the coarse references; full reruns
recompute the finer references as specified in each driver.

The archived `maxBudgetRatio` in the complete-solve CSV was normalized
using 1e-7 for every n. At n=40000 the actual controller tolerance is
1e-6, which is recorded in the separate tau column and used by the solver.
The manuscript Table 4 uses tau, dimension, relative error and time;
it does not use that legacy `maxBudgetRatio` field.
