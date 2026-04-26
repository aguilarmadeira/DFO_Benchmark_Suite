# Continuous-variable problems

Bound-constrained test problems with continuous decision variables, in
two forms:

- `originals/` — 48 distinct problem definitions (function files) collected
  from the standard derivative-free benchmarking literature, in their
  original (unscaled) form.
- `scaled/` — 504 self-contained wrappers obtained by applying eight
  deterministic scaling strategies to the 63 (problem, dimension)
  instances generated from the 48 originals. Each wrapper embeds the
  scaled bounds and the relative-position mapping used to recover the
  original objective. These are the instances that produce the
  experimental results reported in the accompanying paper
  (JCDE-2026-065).
- `scaling_toolkit/` — three MATLAB functions that expose the same
  scaling methodology so it can be applied to user-supplied
  bound-constrained problems.

---

## Note on counts: problem definitions vs instances

The accompanying paper refers to **63 problems** in the test set. This
counts **(problem, dimension) instances** rather than function files,
because some problem definitions are evaluated at multiple dimensions:

| Function file               | Instances | Dimensions covered     |
|-----------------------------|----------:|------------------------|
| `cauchy.m`                  |         2 | n = 4, 10              |
| `cosine_mixture.m`          |         2 | n = 2, 4               |
| `epistatic_michalewicz.m`   |         2 | n = 5, 10              |
| `exponencial.m`             |         2 | n = 2, 4               |
| `fifteenn_local_minima.m`   |         5 | n = 2, 4, 6, 8, 10     |
| `griewank.m`                |         2 | n = 5, 10              |
| `hartman_4.m`               |         2 | n = 3, 6               |
| `salomon.m`                 |         2 | n = 5, 10              |
| `shekel_foxholes.m`         |         2 | n = 5, 10              |
| `tenn_local_minima.m`       |         4 | n = 2, 4, 6, 8         |
| (38 other functions)        |        38 | one dimension each     |
| **Total**                   |    **63** |                        |

This explains the apparent count mismatch between folders:

- `originals/` contains **48 files** (one per problem definition).
- `scaled/<strategy>/` contains **63 files** (one per instance, with the
  dimension fixed and embedded in the wrapper bounds).
- `scaled/` total: **63 × 8 = 504 wrappers**.

---

## Test problems

The 48 distinct problem definitions yield 63 (problem, dimension)
instances spanning dimensions `n ∈ {2, 3, 4, 5, 6, 8, 9, 10}`. Coverage
includes:

| Category         | Examples                                                    |
|------------------|-------------------------------------------------------------|
| Unimodal         | Sphere, Powell, Wood, Rosenbrock                            |
| Multimodal       | Ackley, Rastrigin, Griewank, Schwefel, Langerman            |
| Low-dimensional  | Branin–Hoo, Goldstein–Price, Six-Hump Camel, Hosaki         |
| Shekel family    | Shekel 4-5, 4-7, 4-10, Foxholes (5D, 10D)                   |
| Parametric       | Fifteen / Ten Local Minima (n = 2, 4, 6, 8, 10)             |
| Heterogeneous    | gulf (asymmetric bounds, multiscale by construction)        |

Full per-problem specifications (dimension, bounds, number of local and
global minimizers) are listed in **Appendix A of the accompanying paper**
(JCDE-2026-065). The original references for each problem are
embedded in the docstring of every wrapper file.

---

## Scaling strategies

Eight strategies, organized in three categories:

| Strategy                    | Category               | Contrast (κ)         | Pattern              |
|-----------------------------|------------------------|----------------------|----------------------|
| `baseline`                  | (no scaling)           | 1                    | identity             |
| `moderate`                  | I — Systematic         | 10⁴                  | Cyclic geometric     |
| `progressive`               | I — Systematic         | 10⁶                  | Cyclic geometric     |
| `extreme`                   | I — Systematic         | 10⁸                  | Binary partition     |
| `sobol_oscillatory`         | II — Quasi-random      | 10⁶                  | Continuous QMC       |
| `sobol_digit_oscillatory`   | II — Quasi-random      | 10⁸                  | Binary QMC (digit)   |
| `halton_oscillatory`        | II — Quasi-random      | 10⁶                  | Alternative QMC      |
| `spatial_thermal`           | III — Application      | ≈ 9 × 10⁴            | Multiphysics         |

The mathematical definition of each strategy is given in **Appendix B
of the accompanying paper**. The implementation that produces the
scale factors is in `scaling_toolkit/generate_scale_factors.m`.

---

## Layout of `scaled/`

```
scaled/
├── INDEX.txt                          (inventory of all wrappers)
├── add_strategy.m                     (helper: add one strategy folder to path)
├── add_all_strategies.m               (helper: add all eight folders to path)
├── baseline_1/
├── moderate_1e4/
├── progressive_1e6/
├── extreme_1e8/
├── sobol_oscillatory_1e6/
├── sobol_digit_oscillatory_1e8/
├── halton_oscillatory_1e6/
└── spatial_thermal_9e4/
```

Each strategy folder contains 63 wrapper files of the form
`<problem>_<n>D[_A|_B|...].m`, where the optional letter suffix
disambiguates multiple instances of the same problem at different
dimensions. Wrappers are written in OLD-STYLE MATLAB function syntax
(no closing `end`) to remain compatible with the original problem
implementations they embed.

---

## Self-contained wrapper API

Every wrapper in `scaled/` exposes the same three-mode interface:

```matlab
info       = ackley_10D();        % full metadata struct
[lb, ub]   = ackley_10D(10);      % work-space bounds for the wrapper
f          = ackley_10D(x);       % evaluate at point x ∈ work space
```

The `info` struct includes:

| Field                 | Meaning                                                     |
|-----------------------|-------------------------------------------------------------|
| `name`                | Wrapper file name                                           |
| `problem`             | Source problem name                                         |
| `dimension`           | n                                                           |
| `strategy`            | Scaling strategy used                                       |
| `kappa`               | Target contrast factor                                      |
| `lb_orig`, `ub_orig`  | Original bounds (after the paper's 0.35 perturbation)       |
| `lb_work`, `ub_work`  | Scaled (work-space) bounds passed to the optimizer          |
| `scale_factors`       | Per-coordinate scale factors                                |
| `contrast_ratio`      | Effective max(s)/min(s)                                     |
| `f_global_min`        | Known global minimum value, when documented in literature   |
| `x_global_min_*`      | Known global minimizer(s) in original/work coordinates      |
| `mapping`             | Symbolic description of the work → orig affine mapping      |

Internally, every evaluation maps the input from work space back to
original coordinates via the relative-position formula (paper Eq. B2):

```
t      = clip01( (x - lb_work) ./ (ub_work - lb_work) )
x_orig = lb_orig + t .* (ub_orig - lb_orig)
f      = <embedded original objective>(x_orig)
```

This preserves the landscape structure (multimodality, basins of
attraction, optima locations) up to the affine scaling, while
introducing controlled scale contrast in the work-space coordinates
seen by the optimizer.

---

## Quick start

```matlab
% Add all strategies to the MATLAB path
addpath(genpath('prob_var_continuous/scaled'));

% Inspect a problem
info       = ackley_10D();
[lb, ub]   = ackley_10D(10);

% Evaluate
f = ackley_10D(rand(10,1));
```

Or add a single strategy folder:

```matlab
run('prob_var_continuous/scaled/add_strategy.m');
add_strategy('extreme');
```

---

## Applying the methodology to your own problems

The `scaling_toolkit/` folder provides the three MATLAB functions that
implement the scaling pipeline described in Appendix B of the paper:

```matlab
addpath('prob_var_continuous/scaling_toolkit');

% Your bound-constrained objective, in original variables
fhandle = @my_objective;
lb      = [0; 0; 0; 0];
ub      = [1; 10; 100; 1000];

% Apply the digit-based oscillatory scaling
s             = generate_scale_factors(numel(lb), 'sobol_digit_oscillatory');
[lb_s, ub_s]  = generate_scaled_bounds(lb, ub, s);
wrapper       = create_scaled_wrapper(fhandle, lb, ub, lb_s, ub_s);

% wrapper is a function handle: wrapper(x_work) returns my_objective(x_orig)
f_value = wrapper((lb_s + ub_s)/2);
```

The wrapper preserves the relative landscape of the original
objective; only the geometric representation of the search space is
modified.

---

## Notes on bounds

The `scaled/` wrappers apply two transformations to the original
bounds:

1. **Upper-bound perturbation** (paper Section 4.1): the upper bound
   `u` is replaced by `u − 0.35 (u − ℓ)` to avoid trivial
   initialization at the center of symmetric domains. This convention
   follows Custódio & Madeira (2015).
2. **Scaling**: the perturbed bounds are then transformed according to
   the chosen scaling strategy.

Therefore, `info.lb_orig` and `info.ub_orig` in each wrapper contain
the **perturbed** original bounds, and `info.lb_work` / `info.ub_work`
contain the corresponding scaled bounds. The unperturbed bounds from
the source literature are listed in Appendix A of the paper.

---

## References

Original problem references are embedded in the docstring of each
wrapper file. The most frequent ones are:

- Ali, M. M., Khompatraporn, C., Zabinsky, Z. B. (2005).
  *A numerical evaluation of several stochastic algorithms on selected
  continuous global optimization test problems.*
  Journal of Global Optimization, 31:635–672.
- Brachetti, P., et al. (1997).
  *A new version of the Price's algorithm for global optimization.*
  Journal of Global Optimization, 10:165–184.
- Storn, R., Price, K. (1997).
  *Differential evolution — a simple and efficient heuristic for global
  optimization over continuous spaces.*
  Journal of Global Optimization, 11:341–359.
- Huyer, W., Neumaier, A. (1999).
  *Global optimization by multilevel coordinate search.*
  Journal of Global Optimization, 14:331–355.
