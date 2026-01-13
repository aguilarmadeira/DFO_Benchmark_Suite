# problems_scaled

Self-contained MATLAB benchmark functions with pre-applied scaling transformations for derivative-free optimization research.

## Overview

This directory contains **441 ready-to-use test functions** (63 base problems × 7 scaling strategies). Each `.m` file is completely self-contained: it includes the objective function, bounds, scaling parameters, and known optimum information.

**Key features:**
- **Zero configuration**: Just call the function — no setup required
- **Self-documented**: Each file contains full metadata accessible via `info = problem()`
- **Consistent interface**: All functions follow the same calling convention
- **Reproducible**: Deterministic scaling transformations

## Directory Structure

```
problems_scaled/
├── baseline_1/              # No scaling (κ = 1)
├── progressive_1e6/         # Cyclic geometric scaling (κ = 10⁶)
├── extreme_1e8/             # Binary partition scaling (κ = 10⁸)
├── sobol_oscillatory_1e6/   # Continuous QMC pattern (κ = 10⁶)
├── sobol_digit_oscillatory_1e8/  # Binary QMC pattern (κ = 10⁸)
├── halton_oscillatory_1e6/  # Alternative QMC pattern (κ = 10⁶)
├── spatial_thermal_9e4/     # Multiphysics-inspired (κ ≈ 9×10⁴)
└── INDEX.txt                # Complete file listing
```

Each subfolder contains 63 problem files named `<problem>_<n>D.m`.

## Scaling Strategies

| Strategy | Contrast (κ) | Description |
|----------|--------------|-------------|
| `baseline_1` | 1 | Original bounds, no transformation |
| `progressive_1e6` | 10⁶ | Cyclic geometric progression |
| `extreme_1e8` | 10⁸ | Binary partition (half small, half large) |
| `sobol_oscillatory_1e6` | 10⁶ | Sobol-sequence driven continuous pattern |
| `sobol_digit_oscillatory_1e8` | 10⁸ | Sobol digit-based binary assignment |
| `halton_oscillatory_1e6` | 10⁶ | Halton-sequence driven pattern |
| `spatial_thermal_9e4` | 9×10⁴ | Multiphysics-inspired (spatial + thermal) |

## Usage

### Evaluate the function

```matlab
% Add strategy folder to path
addpath('problems_scaled/extreme_1e8');

% Evaluate at a point
x = rand(10, 1) .* 1e8;  % Point in scaled space
f = ackley_10D(x);
```

### Get bounds

```matlab
% Returns bounds in scaled (working) space
[lb, ub] = ackley_10D(10);  % Pass dimension as argument
```

### Get full metadata

```matlab
% Returns structure with all problem information
info = ackley_10D();

% Available fields:
%   info.name           - Function name
%   info.problem        - Base problem name
%   info.dimension      - Problem dimension
%   info.strategy       - Scaling strategy used
%   info.kappa          - Contrast ratio
%   info.lb_orig        - Original lower bounds
%   info.ub_orig        - Original upper bounds
%   info.lb_work        - Scaled (working) lower bounds
%   info.ub_work        - Scaled (working) upper bounds
%   info.scale_factors  - Per-variable scale factors
%   info.f_global_min   - Known global minimum value
%   info.x_global_min_orig - Global minimizer (original space)
%   info.x_global_min_work - Global minimizer (scaled space)
```

### Run all problems in a strategy

```matlab
% Example: benchmark loop
strategy = 'extreme_1e8';
files = dir(fullfile('problems_scaled', strategy, '*.m'));

for i = 1:length(files)
    fname = files(i).name(1:end-2);  % Remove .m
    func = str2func(fname);
    info = func();
    
    [lb, ub] = func(info.dimension);
    % ... run your optimizer ...
end
```

## Coordinate Mapping

Each wrapper preserves the landscape structure while applying scaling to the decision space. The mapping from scaled coordinates `x` to original coordinates `x_orig` is:

```
t      = clip01((x - lb_work) ./ (ub_work - lb_work))
x_orig = lb_orig + t .* (ub_orig - lb_orig)
f      = f_original(x_orig)
```

This ensures that:
- Global and local minima are preserved in relative coordinates
- The landscape structure (multimodality, basins) remains unchanged
- Only the geometric scaling of the search space is modified

## Test Problems

The 63 base problems include classic benchmarks from the derivative-free optimization literature:

| Dimension | Problems |
|-----------|----------|
| 2D | Ackley, Aluffi-Pentini, Becker-Lago, Bohachevsky, Branin-Hoo, ... |
| 3–4D | Cauchy, Cosine Mixture, Fletcher-Powell, Goldstein-Price, ... |
| 5–6D | Epistatic Michalewicz, Griewank, Shekel, ... |
| 9–10D | Ackley, Cauchy, Langerman, Rastrigin, Schwefel, ... |

See `INDEX.txt` for the complete listing of all 441 instances.

## References

Wrapper/scaling formulation:
> J. F. A. Madeira, "Global and Local Optimization using Direct Search — A Scale-Invariant Approach (GLODS-SI)", *Journal of Global Optimization*, 2026.

Original test functions from:
- Ali, M.M., Khompatraporn, C., Zabinsky, Z.B. (2005). J. Global Optim. 31, 635–672.
- Brachetti, P., et al. (1997). J. Global Optim. 10, 165–184.

## License

Code: [MIT License](../LICENSE.md)
