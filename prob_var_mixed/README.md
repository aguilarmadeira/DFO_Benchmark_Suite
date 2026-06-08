# Mixed-variable benchmark problems

This folder contains the self-contained mixed-variable benchmark problems
used in the GLODS-SI-Mix computational study.

The suite contains **504 MATLAB problem wrappers**, obtained from **63
continuous base problems** and **8 deterministic heterogeneity strategies**.
Each wrapper is fully self-contained: it includes the problem metadata,
variable types, finite domains, categorical labels, categorical embeddings,
decoding function, and the embedded original continuous objective.

**No external generator is required to use these problems.**

Each problem returns a structure `PB`:

```matlab
PB = problem_name();
```

The objective can be evaluated through:

```matlab
f = PB.func_F(x);
```

where `x` is a MATLAB cell array containing continuous values, ordered
discrete values, and categorical labels.

The mixed-variable pattern follows the cyclic structure

```
C, K, D, C, D, K, ...
```

where `C` denotes continuous variables, `D` ordered discrete variables, and
`K` categorical variables.

Categorical variables are exposed as labels, not as ordered numerical
variables. Internally, each wrapper contains the numerical embedding used
only to evaluate the corresponding original continuous benchmark function.

---

## Folder layout

```
prob_var_mixed/
├── README.md          (this file)
├── problemsMix/       504 self-contained wrappers (63 instances × 8 strategies)
│   ├── baseline_1/
│   ├── moderate_1e4/
│   ├── progressive_1e6/
│   ├── extreme_1e8/
│   ├── sobol_oscillatory_1e6/
│   ├── sobol_digit_oscillatory_1e8/
│   ├── spatial_thermal_9e4/
│   ├── halton_oscillatory_1e6/
│   ├── INDEX.txt
│   ├── add_strategy.m
│   └── add_all_strategies.m
└── metadata/
    └── mixed_504_manifest.csv
```

Each strategy folder contains **63** wrapper files named
`<problem>_k<k1>[_k<k2>[_k<k3>]]_<n>D[_A|_B|...].m`, where the `k<·>` fields
list the categorical cardinalities (largest first) and the optional letter
suffix disambiguates several dimensions of the same base problem.

## Quick start

```matlab
addpath(genpath('prob_var_mixed/problemsMix'));   % or run add_all_strategies.m

PB = ackley_k83_k47_10D();        % example wrapper
f  = PB.func_F(PB.x_star);        % objective at a representative optimizer
```

`add_strategy('extreme')` adds a single strategy folder to the path.

## Heterogeneity strategies

The eight strategies are applied **only to the continuous variables** (and
only when an instance has at least two of them; otherwise the wrapper falls
back to `baseline`):

| Folder                          | Strategy                  | Contrast κ |
| ------------------------------- | ------------------------- | ---------- |
| `baseline_1`                    | baseline (no scaling)     | 1          |
| `moderate_1e4`                  | moderate                  | 10⁴        |
| `progressive_1e6`               | progressive               | 10⁶        |
| `extreme_1e8`                   | extreme                   | 10⁸        |
| `sobol_oscillatory_1e6`         | Sobol oscillatory         | 10⁶        |
| `sobol_digit_oscillatory_1e8`   | Sobol digit-oscillatory   | 10⁸        |
| `halton_oscillatory_1e6`        | Halton oscillatory        | 10⁶        |
| `spatial_thermal_9e4`           | spatial thermal           | ≈ 9 × 10⁴  |

## Notes

- **Counts.** The **63** base problems are (problem, dimension) *instances*,
  not function files: some functions appear at several dimensions (e.g.
  `fifteenn_local_minima` at n = 2, 4, 6, 8, 10; `tenn_local_minima` at
  n = 2, 4, 6, 8; `cauchy`, `griewank`, `hartman_4`, `salomon`,
  `shekel_foxholes`, `cosine_mixture`, `epistatic_michalewicz`,
  `exponencial` at two dimensions each). `63 × 8 = 504`.
- **Domains.** Discrete and categorical cardinalities lie in `[20, 130]`;
  discrete grids are strictly increasing and non-uniform; categorical
  embeddings are a non-uniform base scrambled by a fixed permutation, so no
  ordinal information leaks through the labels. Where a continuous optimizer
  coordinate is documented, it is injected exactly into the corresponding
  grid/embedding so the original optimum remains reachable.
- **Bounds.** Wrappers use the **original (literature) bounds**; the 0.35
  upper-bound perturbation of the continuous category is not applied here.
- **metadata/** holds `mixed_504_manifest.csv`, a per-instance inventory
  (problem, `n`, `n_C`/`n_D`/`n_K`, cardinalities, strategy, κ, documented
  `f*`, wrapper file).

## License

Distributed under the **GNU Lesser General Public License v3**
(LGPL-3.0-or-later); see the top-level `LICENSE`. The Cat-Suite comparison
in the paper uses the publicly available CatMADS data
(`github.com/bbopt/CatMADS_prototype`) and is not redistributed here.
