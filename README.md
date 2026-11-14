# Projected nonsymmetric differential Riccati equations

MATLAB code and data for the paper

> Ze Zhang and Jianzhou Liu, *Well-Posedness and A Posteriori Error Bounds
> for Projected Nonsymmetric Differential Riccati Equations*.

The code implements Algorithm 5.1 of the paper, a two-sided rational Krylov
method whose stopping criterion is the a posteriori error bound of the
paper, and the experiments of section 6.

## Requirements

MATLAB R2026b or later; no toolbox is needed.

## Quick start

In MATLAB, in the repository folder:

```matlab
reproduce("check")    % check the recorded results in data/
reproduce("report")   % write Tables 1-5 and Figures 1-3 from data/
reproduce("smoke")    % run Algorithm 5.1 on two small problems
```

`report` writes into a new folder under `reproduced/`.

## Rerunning the experiments

```matlab
reproduce("dense")    % Tables 1-2 and Figure 1, n = 600 (about 4 min)
reproduce("all")      % all experiments, including the reference solutions
```

A complete run took about 1.5 hours on an Intel Core Ultra 9 185H, plus
about 10 minutes for Figure 3.
The reference solutions, unchanged since version 1.1.0, can instead be
downloaded from the
[v1.1.0 release](https://github.com/continuezoz/simax-ndre-reproducibility/releases/tag/v1.1.0)
(four ZIP files, 3.8 GB). Unzip them into the repository folder and run

```matlab
reproduce("evaluators")  % Table 5; needs ndre-reference-transport.zip
reproduce("cached")      % all experiments with the stored references
```

## Tables and figures

| Paper | Script | Data |
|---|---|---|
| Table 1 | `run_bound_experiments("wellposed")` | `data/bounds/wellposed.csv` |
| Table 2 | `run_bound_experiments("dense")` | `data/bounds/dense_bounds.csv` |
| Section 6.2, stationary bound | `run_bound_experiments("stationary")` | `data/bounds/stationary_bounds.csv` |
| Figure 1 | `run_bound_history`, `plot_bound_history` | `data/bounds/history_*.csv` |
| Table 3 | `run_complete_solves` | `data/complete_solves/complete_solves.csv` |
| Figure 2 | `run_floor_traces`, `plot_floor_traces` | `data/floor/` |
| Table 4 | `run_multistate_transport` | `data/multistate_transport/multistate_run.csv` |
| Table 5 | `run_evaluator_comparison` | `data/evaluator_comparison/transport_evaluators.csv` |
| Figure 3 | `run_reflection_design`, `plot_reflection_design` | `data/reflection_design/` |

Errors are relative Frobenius errors in the scaled variables. Times are in
seconds; each solve was run once.

## Contents

- `src/` the method: the scaling of section 2.1 (`ndre_setup`), Algorithm
  5.1 (`rksm_bound_controller`), the error estimate (`error_estimate`), the
  increment evaluation (`incremental_centered_nodes`), scaled Radon
  composition (`scaled_radon_nodes`), the two transport models, and the
  EBA--BDF1 baseline (`ndre_core`, and `eba_deflated_bdf` with deflation for
  section 6.6).
- `experiments/` one script per table or figure.
- `data/` the recorded results and the Gauss--Legendre rules.
- `figures/` Figures 1-3 of the paper.

To call the functions directly, run `addpath src experiments` first.

## License

MIT; see `LICENSE`.
