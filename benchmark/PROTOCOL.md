# DFO Benchmark Suite v2 — evaluation protocol

This protocol applies to every solver and uses only the sequence of evaluated
points and values, so that results obtained by different people on the same
folder are comparable without re-running anything. Everything below is
implemented in `benchmark/`; a solver only has to call the objective.

Indicators: evaluations to targets (Moré–Wild relative levels, and absolute
levels where f* = 0), best value reached, robustness across instances (data
and performance profiles), evaluation cost (`eval_us` in `INDEX.csv`) and,
for folders with categorical variables, path statistics (`analyze_interact`).

## Problems

A *folder* is one sub-folder of `core/` (`C_baseline`, `D`, `CDK_extreme`,
`BK_MK16`, ...; see `docs/CORE_SPEC.md`). A *problem* is one (instance, seed)
pair of a folder. Deterministic solvers run once per instance (seed 1) and
that run stands for every seed; stochastic solvers run with seeds 1..10.

## Interface seen by the solver

`bench_problem` hands the solver a struct `P` with

| field | meaning |
|---|---|
| `P.nvar`, `P.var_types` | dimension and types (`'C'`, `'D'`, `'K'`) |
| `P.cardinality` | number of levels per variable (`Inf` for continuous) |
| `P.ProblemData` | the domain in the cell format of the wrappers |
| `P.lb_C_work`, `P.ub_C_work` | bounds of the continuous variables |
| `P.func_idx(z)` | objective; `z` is a vector of **1-based level indices** (continuous coordinates as values) |
| `P.func_F(xr)` | objective; `xr` is a cell of values / grid values / labels |
| `P.idx2cell`, `P.cell2idx` | conversions (not counted) |
| `P.x0_idx`, `P.x0` | the starting point of the protocol (below) |
| `P.budget`, `P.seed` | budget of this run and seed to use for any randomness |

Both objectives are counted. Nothing else about the problem (optimal value,
optimal solution, construction) is passed to the solver. Categorical labels
carry no order and no distance; a solver must not use the index of a label as
a number.

## Starting point

Every wrapper carries `PB.x0_idx`, a Halton point fixed when the problem was
generated and never chosen after looking at the function (index 1000 + p for
the classical instances, p = index of the base instance; 2000 / 3000 + 10 st +
seed for Family BK, shared by all λ and b of a structure and seed). Nine
further points `PB.x0_alt_idx` serve multi-start protocols; report which ones
were used. `run_benchmark` evaluates `x0` before calling the solver, so every
trace starts with `(1, f(x0))` and `f(x0)` is the same for every solver; a
solver may start elsewhere, but then the convergence test below still uses
`f(x0)` of the protocol. An instance whose `x0` is an optimum would be
*trivial* and excluded from the profiles; no instance of the core is trivial
(`INDEX.csv`, column `x0_is_optimum`).

## Budget

`budget = 1000 n` evaluations per run (option `budget_per_n`; report it when
it differs). Evaluations beyond the budget are not counted: `run_benchmark`
stops the solver through the error `BENCH:budget`. Data profiles are drawn
against `evaluations / (n+1)`, so folders of different dimensions share one
axis.

## Convergence test

Moré & Wild (2009): a run has solved the problem at level τ after `t`
evaluations when its best value satisfies

    f ≤ f_L + τ (f(x0) − f_L),      τ ∈ {1e-1, 1e-3, 1e-5, 1e-7},

where `f_L` is the target of the instance, read from the versioned file
`targets/<folder>_targets.csv` (instance, f_target, source, date). In the core
every target is the known global optimum `PB.f_star` (`source = f_star`), or
the same value certified by enumeration of the finite domain
(`source = certified`). τ = 1e-7 on an instance means the optimum was found.
For Family BK (f* = 0) report also the absolute levels b, b/10, 1e-6, 1e-9
(`analyze_interact`).

## Known optimum

Every instance has a known global optimum: `PB.f_star` and its justification
`PB.f_star_justification` — `E` exact by structure or by construction, `A`
analytic reduction, `L` the established literature value at a refined
minimiser (not proved; for the fully finite instances with at most 2e6
points it is certified by enumeration, `INDEX.csv` column `certified`). The
optimal solution is in `reference_optima.csv` of each folder, used by
`validate_suite` and `analyze_interact` only.

## Invalid evaluations

A point where the function returns NaN or ±Inf counts as one evaluation and
is returned to the solver as `+Inf`. `summary.csv` reports how many such
evaluations each run made (`n_nonfinite`). In the core this happens only for
fletcher_powell (D, K, DK, 3D) on the grid plane x1 = 0, where the original
function is +Inf by definition.

## Held-out part

`holdout_instances.m`: 19 base instances (by base function and dimension, in
every pattern and strategy folder) are held out. Results on them are reported
separately and are never used to choose parameters or variants of a solver.

## Files

```
core/<folder>/INDEX.csv                              per instance: n, cardinalities, f*, certified, x0_is_optimum, ...
core/<folder>/reference_optima.csv                   optimal solution per instance (index and original coordinates)
targets/<folder>_targets.csv                         versioned targets f_L
results/<solver>/<folder>/<instance>__s<seed>.csv    trace: neval, fbest at every improvement
results/<solver>/<folder>/summary.csv                one row per run (nevals = calls, n_distinct = distinct points)
profiles/<folder>_tau<τ>_data.csv / _perf.csv / .png data and performance profiles
```

To publish results, ship the `results/<solver>/<folder>/` trees: anyone can
add a solver and redraw the profiles with `data_profiles`.

## Minimal session

```matlab
addpath('benchmark');
run_benchmark(solver_random_search(), 'core/D');                 % 10 seeds, 1000 n
data_profiles({'random_search'}, 'core/D');                      % uses targets/D_targets.csv
```

A solver is a struct `solver.name`, `solver.fun = @(P) ...` (and
`solver.deterministic = true` when applicable); `solver_random_search.m` is
the smallest example. `validate_suite('core/D')` re-checks a folder
(admissibility, F(x*) = f*, enumeration) and `finalize_core` re-checks all.

## Paired analysis across folders

The same base instances appear in several folders (patterns, strategies).
Compare folders on the base instances common to them: `[bases, per_folder] =
paired_nontrivial(folders)` returns, per folder, the instance names to pass as
the `instances` option of `data_profiles`; otherwise the effect of the folder
is confounded with a change of the function set. In Family BK, instances that
differ only in λ or b share x*, labels and x0, so the effect of the
interaction parameters can be read directly.

## What the protocol does not cover

Wall-clock time (reported in `summary.csv`, not used in the profiles),
constraints other than bounds, and noise (the functions are deterministic).

## References

J. J. Moré, S. M. Wild, Benchmarking derivative-free optimization
algorithms, SIAM J. Optim. 20 (2009) 172–191.
E. D. Dolan, J. J. Moré, Benchmarking optimization software with performance
profiles, Math. Program. 91 (2002) 201–213.
