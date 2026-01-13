# Scaling Transformation Utilities

This directory contains MATLAB functions for applying deterministic variable scaling transformations to benchmark problems.

## Overview

Scale heterogeneity across variables is a major challenge in derivative-free optimization. These utilities allow systematic testing of algorithm robustness by applying controlled scaling transformations while preserving problem structure.

## Functions

### `generate_scale_factors.m`

Generates per-variable scale factors for a given strategy.

```matlab
s = generate_scale_factors(n, strategy)
s = generate_scale_factors(n, strategy, 'MaxContrast', kappa)
```

**Inputs:**
- `n` — Problem dimension
- `strategy` — Scaling strategy name (see below)
- `kappa` — Target contrast ratio (optional)

**Output:**
- `s` — Vector of scale factors (n × 1)

### `generate_scaled_bounds.m`

Applies scale factors to transform problem bounds.

```matlab
[lb_scaled, ub_scaled] = generate_scaled_bounds(lb, ub, s)
```

**Inputs:**
- `lb`, `ub` — Original lower and upper bounds
- `s` — Scale factors from `generate_scale_factors`

**Output:**
- `lb_scaled`, `ub_scaled` — Transformed bounds

## Scaling Strategies

### Category I: Systematic Geometric Scaling

| Strategy | Contrast (κ) | Pattern |
|----------|--------------|---------|
| `baseline` | 1 | No scaling |
| `progressive` | 10⁶ | Cyclic geometric: exponents cycle through {0,1,2,3,-1,-2,-3} |
| `extreme` | 10⁸ | Binary partition: half variables at κ^(-1/2), half at κ^(1/2) |

### Category II: Quasi-Random Oscillatory Scaling

| Strategy | Contrast (κ) | Sequence |
|----------|--------------|----------|
| `sobol_oscillatory` | 10⁶ | Sobol low-discrepancy, continuous mapping |
| `sobol_digit_oscillatory` | 10⁸ | Sobol with digit-based binary assignment |
| `halton_oscillatory` | 10⁶ | Halton low-discrepancy sequence |

### Category III: Application-Inspired Scaling

| Strategy | Contrast (κ) | Motivation |
|----------|--------------|------------|
| `spatial_thermal` | ≈9×10⁴ | Multiphysics: spatial coordinates vs. thermal properties |

## Mathematical Details

### Scale Factor Generation

For quasi-random strategies, scale factors are computed as:

```
s_i = 10^(6*q_i - 3)
```

where `q_i ∈ [0,1]` comes from the low-discrepancy sequence, yielding `s_i ∈ [10^(-3), 10^3]`.

### Bounds Transformation

Given original bounds `[ℓ, u]` and scale factors `s`, the scaled bounds are:

```
c_i = (ℓ_i + u_i) / 2           % center
w_i = u_i - ℓ_i                  % width
ℓ̃_i = c_i * s_i - 0.5 * s_i * w_i
ũ_i = c_i * s_i + 0.5 * s_i * w_i
```

### Objective Wrapper

To preserve problem structure, evaluate via relative coordinates:

```matlab
function f = wrapper(x_scaled, f_orig, lb, ub, lb_scaled, ub_scaled)
    t = (x_scaled - lb_scaled) ./ (ub_scaled - lb_scaled);  % relative position
    x_orig = lb + t .* (ub - lb);                            % map to original
    f = f_orig(x_orig);
end
```

## Example Usage

```matlab
% Original problem setup
lb = -10 * ones(5, 1);
ub =  10 * ones(5, 1);
f_orig = @rosenbrock;

% Apply progressive scaling (κ = 10^6)
s = generate_scale_factors(5, 'progressive', 'MaxContrast', 1e6);
[lb_s, ub_s] = generate_scaled_bounds(lb, ub, s);

% Evaluate at scaled point
x_scaled = lb_s + rand(5,1) .* (ub_s - lb_s);
t = (x_scaled - lb_s) ./ (ub_s - lb_s);
x_orig = lb + t .* (ub - lb);
f = f_orig(x_orig);
```

## References

1. Hansen, N., Auger, A., Finck, S., Ros, R. (2009). Real-parameter black-box optimization benchmarking 2009: Experimental setup. *GECCO Workshop on BBOB*.

2. Madeira, J.F.A. (2025). GLODS-SI: Global and Local Optimization using Direct Search – A Scale-Invariant Approach. *J. Global Optim.* (submitted)
