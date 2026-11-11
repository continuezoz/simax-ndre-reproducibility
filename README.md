# Reproducibility package for projected nonsymmetric differential Riccati equations

This repository accompanies **Well-Posedness and A Posteriori Error Bounds
for Projected Nonsymmetric Differential Riccati Equations**, by **Ze Zhang
and Jianzhou Liu**. Version **v1.0.0** contains the code, parameter choices,
recorded numerical data, and MATLAB figure scripts for the SIMAX manuscript
dated 11 November 2026. It is provided to support a request for the SIMAX
code-and-data-available reproducibility badge.

Repository: https://github.com/continuezoz/simax-ndre-reproducibility

## Quick start

Open MATLAB in the repository directory:

```matlab
reproduce_simax("verify")       % inspect the recorded numerical results
folder = reproduce_simax("report"); % regenerate all five tables and both figures
reproduce_simax("smoke")        % solve two small instances from their generators
```

The first two commands need only the repository download. `report` creates
a separate directory under `reproduced/` and prints its location. Its five
tables appear as CSV and LaTeX tabular files under `generated/tables/`;
the figures appear as PDF and PNG under `results_scaled/figures/`.
The original recorded results are retained in the repository.

## Full numerical reproduction

```matlab
folder = reproduce_simax("all");
```

This creates a separate working directory, generates both model families,
computes the coarse and finer reference solutions, and runs all experiments
used by the manuscript. No external dataset is needed. All parameters are in
the supplied `.m` files. The original solver routines are copied unchanged
from the working project; `SOURCE_DATA_MANIFEST.csv` records their SHA-256
hashes and the identities of the recorded data.

For a faster rerun using the archived reference solutions, download the four
`simax-reference-*.zip` assets from the **v1.0.0 GitHub Release**, and unzip
each into the repository directory. Each archive has paths beginning with
`results_scaled/`. Then run:

```matlab
folder = reproduce_simax("cached");
```

Both numerical modes recompute the proposed method and the baselines. The
`cached` mode reuses the stored coarse reference trajectories; `all` also
recomputes and compares their finer-space counterparts. The original
fine-space trajectories were used for refinement checks and are regenerated
by `all`; their bases and algorithms are specified in the drivers.

Individual verification runs are also available:

```matlab
reproduce_simax("dense")       % all n=600 bound/solvability tests and Figure 1 data
reproduce_simax("evaluators")  % recompute Table 6; needs the transport reference archive
```

## Paper-to-code map

| Paper item | Numerical driver | Recorded data / figure producer |
|---|---|---|
| Table 2: projected stationary solvability | `run_bound_experiments("wellposed",600)` | `results_scaled/bounds/wellposed.csv` |
| Table 3: dense-reference error bounds | `run_bound_experiments("dense",600)` | `results_scaled/bounds/dense_bounds.csv`; select tolerances 1e-4, 1e-7, 1e-9 |
| Section 6.2: stationary error bounds | `run_bound_experiments("stationary",600)` | `results_scaled/bounds/stationary_bounds.csv` |
| Figure 1: error-bound histories | `run_bound_history(600,1e-7)` | `bounds/history_*.csv`; `plot_bound_history` |
| Table 4: large single-state solves | `run_complete_solves("run")` | `complete_solves/complete_solves.csv`; Newton maxima in `run_run.log` |
| Figure 2: floating-point floors | `run_floor_traces` | `floor/trace_n*.csv`, `floor_summary.csv`; `plot_floor_traces` |
| Table 5: multistate solves | `run_multistate_transport("run")` | `multistate_transport/multistate_run.csv` |
| Table 6: reduced-flow evaluators | `run_evaluator_comparison("transport")` | `evaluator_comparison/transport_evaluators.csv` |

`export_simax_tables` converts these data to five self-contained numeric
tables. Table 1 in the manuscript is a notation table, not a numerical
experiment. The manuscript typesets Tables 2--6 by hand; exported tables
contain the corresponding numerical entries, with their own descriptive
column headers. `PAPER_RESULTS.md` explains aggregation and provenance.

## Environment and cost

* Tested MATLAB: R2026b, 26.2.0.3386108, 64-bit Windows (`PCWIN64`).
* Dependency analysis identifies base MATLAB only. No optional toolbox,
  external executable, Python runtime, or network access is needed to run
  the MATLAB experiments. Python scripts outside this repository were used
  only to prepare publication archives and upload them.
* Paper hardware: Intel Core Ultra 9 185H. Each solve was timed once.
  Historical CSV columns named `medianSeconds`, `minSeconds`, and
  `maxSeconds` consequently contain that same single timing.
* The original four-stage run took about 218 s (dense tests), 1916 s
  (single-state solves), 2596 s (multistate solves including references),
  and 5 s (evaluator comparison). Figure-history and floor runs are extra.
  A full rerun should be budgeted at roughly 1--2 hours; time depends on
  the machine. The report-only mode takes seconds.
* Archived coarse reference data occupy approximately 3.8 GiB. A cached
  rerun copies these files to its isolated working directory; allow space
  for both copies, the downloaded archives, and numerical temporaries.
  Peak memory has not been profiled.
* The deterministic Gauss--Legendre rules used for the paper are included
  under `results/matlab_quadrature_20261007/`; absent rules are generated
  by `gauss_legendre_matlab`. The paper's transport and multistate model
  generators do not use pseudorandom data. Auxiliary model families in
  `nla_core` set their own seed explicitly.

## Parameters and interpretation

The single-state model has c=0.9999, alpha=1e-8, with n=4000, 10000, 20000,
40000; dense tests also use c=alpha=0.5. The multistate model uses rank-10
scattering, gamma=0.9999 or 0.5, asymmetry 0.1, and n=10000 or 100000.
All experiments have X(0)=0 and final time 1.

The proposed controller uses tau_a=1e-12 and normally tau=1e-7; the
single-state n=40000 solve uses tau=1e-6. Dense tests sweep 1e-4 through
1e-9. EBA--BDF1 uses h=0.01, endpoint tolerance 1e-10, and Newton target
1e-11, with stagnated steps accepted below backward error 1e-9. The four
multistate baseline rows labelled BREAKDOWN are reported outcomes of the
published method without deflation, not missing observations. NaN entries
in those rows represent quantities not obtained after breakdown.

The error bound concerns projection error of the exact reduced flow.
Computed-output error also includes floating-point evaluation error.
`dense_bounds.csv` and `history_*.csv` keep those quantities separate.
Error evaluations omit t=0 and compare BDF1 only at its available nodes.
All paper errors use the scaled variables and relative Frobenius norms.

`verify` checks the archived data against the stated outcomes; it does not
rerun the solvers. `smoke`, `dense`, `evaluators`, `all`, and `cached`
perform numerical computations. New runs can differ at rounding level,
especially near the n=40000 accuracy floor, and timings depend on hardware.
The validation record states exactly which modes were executed.

## Files, archive identity, and availability

`FILES.txt` describes the source files; `SOURCE_DATA_MANIFEST.csv` identifies
all copied working sources and data, including large reference files.
`REPRODUCIBILITY_VALIDATION.md` records the release validation. Generated
tables and rerun outputs are intentionally ignored by Git.

GitHub Release v1.0.0 provides a code/data snapshot ZIP and four reference
data ZIPs. The same five archives are prepared for SIMAX supplementary
materials, with a SHA-256 checksum list and supplementary-material index.
Download the code snapshot into an empty directory and extract the data
archives into that same directory to obtain the complete frozen package.
Subsequent code changes should use a new release version.

## Copyright

Copyright (c) 2026 Ze Zhang and Jianzhou Liu. The authors make this code and
data publicly accessible for inspection and reproduction of the reported
experiments. A separate general-purpose open-source reuse license has not
been specified. See `COPYRIGHT.txt`.
