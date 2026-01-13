# DFO Benchmark Suite

A comprehensive benchmark suite for **Derivative-Free Global Optimization**, containing test problems, scaling strategies, and pre-computed initial samples.

[![License: MIT](https://img.shields.io/badge/Code-MIT-blue.svg)](https://opensource.org/licenses/MIT)
[![License: CC BY 4.0](https://img.shields.io/badge/Data-CC%20BY%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)
[![DOI](https://img.shields.io/badge/DOI-pending-orange.svg)](https://doi.org/)

## Overview

This repository provides resources for benchmarking derivative-free optimization algorithms:

- **63 test problems** with dimensions n ∈ {2, 3, 4, 5, 6, 8, 9, 10}
- **7 scaling strategies** with contrast ratios from κ=1 to κ=10⁸
- **Pre-computed initial samples** (60×n Sobol points per problem)
- **HDF5 format** readable in MATLAB, Python, Julia, R

## Repository Structure

```
DFO_Benchmark_Suite/
├── README.md                 # This file
├── LICENSE.md                # Dual license (MIT + CC-BY 4.0)
├── CITATION.cff              # Citation metadata for GitHub
├── problems/                 # Test problem definitions
│   ├── README.md
│   └── *.m                   # MATLAB function files
├── scaling/                  # Scaling transformation utilities
│   ├── README.md
│   ├── generate_scale_factors.m
│   └── generate_scaled_bounds.m
└── DataSamplesAlfa_ini/      # Pre-computed samples
    ├── README.md
    ├── index.xlsx
    ├── baseline_1/
    ├── progressive_1e6/
    ├── extreme_1e8/
    ├── sobol_oscillatory_1e6/
    ├── sobol_digit_oscillatory_1e8/
    ├── halton_oscillatory_1e6/
    └── spatial_thermal_9e4/
```

## Quick Start

### Python
```python
import h5py

with h5py.File('DataSamplesAlfa_ini/baseline_1/ackley_10D.h5', 'r') as f:
    Y = f['/sample/Y'][:]    # Normalized sample points [0,1]^n
    F = f['/sample/F'][:]    # Objective values
    n = f['/meta'].attrs['dimension']
```

### MATLAB
```matlab
Y = h5read('DataSamplesAlfa_ini/baseline_1/ackley_10D.h5', '/sample/Y');
F = h5read('DataSamplesAlfa_ini/baseline_1/ackley_10D.h5', '/sample/F');
n = h5readatt('DataSamplesAlfa_ini/baseline_1/ackley_10D.h5', '/meta', 'dimension');
```

## Test Problems

The benchmark includes 63 problems commonly used in derivative-free global optimization:

| Category | Problems |
|----------|----------|
| Unimodal | Sphere, Powell, Wood, Rosenbrock |
| Multimodal | Ackley, Rastrigin, Griewank, Schwefel, Langerman |
| Low-dimensional | Branin-Hoo, Goldstein-Price, Six-Hump Camel, Hosaki |
| Shekel family | Shekel 4-5, 4-7, 4-10, Foxholes |
| Parametric | Fifteen/Ten Local Minima (n=2,4,6,8,10) |

See [`problems/README.md`](problems/README.md) for complete specifications.

## Scaling Strategies

| Strategy | Contrast (κ) | Description |
|----------|--------------|-------------|
| `baseline` | 1 | Original bounds (no scaling) |
| `progressive` | 10⁶ | Cyclic geometric pattern |
| `extreme` | 10⁸ | Binary partition (worst-case) |
| `sobol_oscillatory` | 10⁶ | Continuous quasi-random |
| `sobol_digit_oscillatory` | 10⁸ | Binary quasi-random |
| `halton_oscillatory` | 10⁶ | Alternative QMC sequence |
| `spatial_thermal` | ≈9×10⁴ | Multiphysics-inspired |

See [`scaling/README.md`](scaling/README.md) for implementation details.

## Citation

If you use this benchmark suite, please cite:

```bibtex
@misc{madeira2026dfobenchmark,
  author       = {Madeira, Jos{\'e} F. Aguilar},
  title        = {{DFO} Benchmark Suite: Test Problems and Initial Samples 
                  for Derivative-Free Optimization},
  year         = {2026},
  publisher    = {GitHub},
  url          = {https://github.com/aguilarmadeira/DFO_Benchmark_Suite},
  note         = {Version 1.0.0}
}
```

### Related Publications

- **GLODS**: Custódio, A.L., Madeira, J.F.A. (2015). GLODS: Global and Local Optimization using Direct Search. *J. Global Optim.* 62, 1–28. [DOI:10.1007/s10898-014-0224-9](https://doi.org/10.1007/s10898-014-0224-9)

- **GLODS-SI**: Madeira, J.F.A. (2025). GLODS-SI: Global and Local Optimization using Direct Search – A Scale-Invariant Approach. *J. Global Optim.* (submitted)

- **Auto-Init**: Madeira, J.F.A., Thomson, S.L. (2025). Automatic Step-Size Initialization for Derivative-Free Global Optimization via Fitness Landscape Analysis. *J. Global Optim.* (in preparation)

## License

This repository uses a **dual license**:

| Component | License | Files |
|-----------|---------|-------|
| Source code | [MIT](https://opensource.org/licenses/MIT) | `.m`, `.py` |
| Data & documentation | [CC-BY 4.0](https://creativecommons.org/licenses/by/4.0/) | `.h5`, `.xlsx`, `.md` |

**Attribution is required when using the data.** See [LICENSE.md](LICENSE.md) for details.

## Acknowledgments

This work was supported by Fundação para a Ciência e a Tecnologia (FCT) through LAETA (project [UID/50022/2025](https://doi.org/10.54499/UID/50022/2025)).

## Contact

**José F. Aguilar Madeira**  
IDMEC, Instituto Superior Técnico, Universidade de Lisboa  
ISEL, Instituto Politécnico de Lisboa  
Email: aguilarmadeira@tecnico.ulisboa.pt  
ORCID: [0000-0001-9523-3808](https://orcid.org/0000-0001-9523-3808)
