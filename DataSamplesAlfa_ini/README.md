# DataSamplesAlfa_ini

Deterministic initial samples for derivative-free optimization benchmark problems under various scaling strategies.

## Overview

This dataset provides pre-computed initial samples for 63 bound-constrained test problems commonly used in derivative-free global optimization research. The samples are generated using quasi-random Sobol sequences in normalized coordinates [0,1]^n and include objective function evaluations under different variable scaling strategies.

**Key features:**
- **Deterministic and reproducible**: All samples use fixed Sobol sequences with consistent skip/leap parameters
- **Budget-neutral**: 60×n sample points per problem (n = problem dimension)
- **Multiple scaling regimes**: 7 strategies ranging from baseline (κ=1) to extreme heterogeneity (κ=10^8)
- **Ready-to-use**: HDF5 format readable in MATLAB, Python, Julia, R, and other languages

## Directory Structure

```
DataSamplesAlfa_ini/
├── README.md                          # This file
├── index.xlsx                         # Master index with all instances
├── baseline_1/                        # Baseline scaling (κ=1)
│   ├── ackley_10D.h5
│   ├── rosenbrock_2D.h5
│   └── ...
├── progressive_1e6/                   # Progressive scaling (κ=10^6)
├── extreme_1e8/                       # Extreme scaling (κ=10^8)
├── sobol_oscillatory_1e6/             # Sobol oscillatory (κ=10^6)
├── sobol_digit_oscillatory_1e8/       # Sobol digit-driven (κ=10^8)
├── halton_oscillatory_1e6/            # Halton oscillatory (κ=10^6)
└── spatial_thermal_9e4/               # Spatial-thermal (κ≈9×10^4)
```

## Scaling Strategies

| Strategy                  | Contrast (κ) | Description                          |
|---------------------------|--------------|--------------------------------------|
| `baseline`                | 1            | No scaling (original bounds)         |
| `progressive`             | 10^6         | Cyclic geometric pattern             |
| `extreme`                 | 10^8         | Binary partition (worst-case)        |
| `sobol_oscillatory`       | 10^6         | Continuous quasi-random pattern      |
| `sobol_digit_oscillatory` | 10^8         | Binary quasi-random pattern          |
| `halton_oscillatory`      | 10^6         | Alternative quasi-random sequence    |
| `spatial_thermal`         | ≈9×10^4      | Application-inspired (multiphysics)  |

## HDF5 File Structure

Each `.h5` file contains:

```
/meta/
    problem          (attr)  : Problem name (e.g., "ackley")
    instance         (attr)  : Instance identifier (e.g., "ackley_10D")
    strategy         (attr)  : Scaling strategy name
    dimension        (attr)  : Problem dimension n
    kappa            (attr)  : Target contrast factor
    contrast_ratio   (attr)  : Actual max(s)/min(s) ratio
    lb_orig          (n,)    : Original lower bounds
    ub_orig          (n,)    : Original upper bounds (after 0.35 perturbation)
    lb_work          (n,)    : Scaled lower bounds
    ub_work          (n,)    : Scaled upper bounds
    scale_factors    (n,)    : Per-variable scale factors

/sample/
    Y                (m, n)  : Sample points in normalized [0,1]^n
    X_work           (m, n)  : Sample points in scaled domain
    F                (m,)    : Objective values f(x_orig)

/results/
    alpha_merge      (attr)  : Computed initial step size (merge-based)
    radius_merge     (attr)  : Corresponding comparison radius
    alpha_status     (attr)  : "ok", "fallback_small_sample", or "fallback_error"
```

## Usage Examples

### Python (h5py)
```python
import h5py
import numpy as np

with h5py.File('baseline_1/ackley_10D.h5', 'r') as f:
    Y = f['/sample/Y'][:]           # (m, n) normalized sample
    F = f['/sample/F'][:]           # (m,) objective values
    n = f['/meta'].attrs['dimension']
    problem = f['/meta'].attrs['problem'].decode()
    
print(f"Problem: {problem}, n={n}, m={len(F)}")
print(f"f_min={F.min():.4f}, f_max={F.max():.4f}")
```

### MATLAB
```matlab
info = h5info('baseline_1/ackley_10D.h5');
Y = h5read('baseline_1/ackley_10D.h5', '/sample/Y');
F = h5read('baseline_1/ackley_10D.h5', '/sample/F');
n = h5readatt('baseline_1/ackley_10D.h5', '/meta', 'dimension');

fprintf('Dimension: %d, Sample size: %d\n', n, numel(F));
```

### Julia
```julia
using HDF5

h5open("baseline_1/ackley_10D.h5", "r") do f
    Y = read(f, "/sample/Y")
    F = read(f, "/sample/F")
    n = read_attribute(f["/meta"], "dimension")
    println("n=$n, m=$(length(F))")
end
```

## Index File

The `index.xlsx` file provides a master table with all 441 instances (63 problems × 7 strategies):

| Column           | Description                                      |
|------------------|--------------------------------------------------|
| `strategy_folder`| Subdirectory name                                |
| `instance`       | Instance identifier                              |
| `problem`        | Problem function name                            |
| `D`              | Dimension                                        |
| `strategy`       | Scaling strategy                                 |
| `kappa`          | Target contrast                                  |
| `contrast_ratio` | Actual contrast ratio                            |
| `m`              | Sample size (60×D)                               |
| `f_min`, `f_med`, `f_max`, `f_std` | Sample statistics    |
| `alpha_merge`    | Computed step size                               |
| `radius_merge`   | Computed comparison radius                       |
| `alpha_status`   | Computation status                               |
| `h5_path`        | Relative path to HDF5 file                       |

## Sample Generation Details

- **Sequence**: Sobol quasi-random sequence
- **Sample size**: m = 60 × n points per problem
- **Skip**: n + 1 (dimension-dependent)
- **Leap**: 2^n
- **Bound perturbation**: Upper bounds reduced by 35% of range to avoid trivial central minima

## Test Problems

The benchmark includes 63 problems with dimensions n ∈ {2, 3, 4, 5, 6, 8, 9, 10}:

Ackley, Aluffi-Pentini, Becker-Lago, Bohachevsky, Branin-Hoo, Cauchy (n=4,10), 
Cosine Mixture (n=2,4), Dekkers-Aarts, Epistatic Michalewicz (n=5,10), 
Exponential (n=2,4), Fifteen Local Minima (n=2,4,6,8,10), Fletcher-Powell, 
Goldstein-Price, Griewank (n=5,10), Gulf, Hartman (n=3,6), Hosaki, Kowalik, 
Langerman, McCormick, Meyer-Roth, Miele-Cantrell, Multi-Gaussian, Neumaier 2, 
Neumaier 3, Paviani, Periodic, Poissonian, Powell, Rastrigin, Rosenbrock, 
Salomon (n=5,10), Schaffer 1, Schaffer 2, Schwefel, Shekel (4-5, 4-7, 4-10), 
Shekel Foxholes (n=5,10), Shubert, Sinusoidal, Six-Hump Camel, Sphere, 
Storn-Tchebychev, Ten Local Minima (n=2,4,6,8), Three-Hump Camel, Transistor, Wood.

## References

If you use this dataset, please cite:

```bibtex
@article{madeira2025glods_si,
  title={{GLODS-SI}: Global and Local Optimization using Direct Search 
         -- A Scale-Invariant Approach},
  author={Madeira, J. F. A.},
  year={2025},
  note={Submitted}
}

@article{madeira2025autoinit,
  title={Automatic Step-Size Initialization for Derivative-Free Global Optimization
         via Fitness Landscape Analysis},
  author={Madeira, J. F. A. and Thomson, Sarah L.},
  year={2025},
  note={In preparation}
}
```

## License

This repository uses a dual license:
- **Code** (`.m` scripts): [MIT License](https://opensource.org/licenses/MIT)
- **Data** (`.h5`, `.xlsx` files): [CC-BY 4.0](https://creativecommons.org/licenses/by/4.0/)

**Attribution is required when using the data.** See LICENSE.md for full details.

## Contact

J. F. A. Madeira  
IDMEC, Instituto Superior Técnico, Universidade de Lisboa  
Email: aguilarmadeira@tecnico.ulisboa.pt

## Version History

- **v1.0** (January 2026): Initial release with 63 problems × 7 scaling strategies
