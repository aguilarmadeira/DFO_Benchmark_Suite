# DFO Benchmark Suite

A benchmark suite for **Derivative-Free Optimization (DFO)**, organized by the
type of decision variables. The current release provides bound-constrained
problems with **continuous variables** under controlled scale heterogeneity;
additional categories (discrete, categorical, mixed) are planned for future
releases.

[![License: LGPL v3](https://img.shields.io/badge/License-LGPL_v3-blue.svg)](https://www.gnu.org/licenses/lgpl-3.0)
[![DOI](https://img.shields.io/badge/DOI-pending-orange.svg)](https://doi.org/)

---

## Repository structure

```
DFO_Benchmark_Suite/
├── README.md                    (this file)
├── LICENSE                      (GNU LGPL v3)
├── CITATION.cff                 (citation metadata)
│
└── prob_var_continuous/         Category: continuous variables
    ├── README.md                (category-specific documentation)
    ├── originals/               48 original problem definitions (Ali, Brachetti, ...)
    ├── scaled/                  504 self-contained wrappers (63 instances × 8 strategies)
    └── scaling_toolkit/         3 MATLAB functions to scale new problems
```

Future releases may add:

```
├── prob_var_discrete/           (planned)
├── prob_var_categorical/        (planned)
└── prob_var_mixed/              (planned)
```

---

## Current release: continuous-variable problems

The `prob_var_continuous/` category contains:

- **48 original problem definitions** corresponding to **63 (problem,
  dimension) instances** (some functions are evaluated at multiple
  dimensions, e.g., `cauchy` at `n=4` and `n=10`). Dimensions span
  `n ∈ {2, 3, 4, 5, 6, 8, 9, 10}`. Problems are collected from the
  standard derivative-free benchmarking literature
  (Ali et al. 2005; Brachetti et al. 1997; Storn & Price 1997; and others).
- **504 self-contained MATLAB wrappers** (63 instances × 8 scaling strategies),
  in which heterogeneous variable scales have been embedded by construction.
  Each wrapper is fully self-contained at runtime: no external scaling
  utilities or path setup are required to evaluate the objective.
- **A standalone scaling toolkit** of three MATLAB functions
  (`generate_scale_factors`, `generate_scaled_bounds`,
  `create_scaled_wrapper`) that allows users to apply the same scaling
  methodology to their own bound-constrained problems.

For full details, scaling strategies, file layout, and usage examples,
see [`prob_var_continuous/README.md`](prob_var_continuous/README.md).

---

## Quick start

Add the relevant strategy folder and (optionally) the scaling toolkit to
the MATLAB path, and call any wrapper as a regular MATLAB function:

```matlab
addpath(genpath('prob_var_continuous/scaled'));

% Get info, bounds, evaluate
info       = ackley_10D();
[lb, ub]   = ackley_10D(10);
f          = ackley_10D(rand(10,1));
```

---

## License

This repository is distributed under the **GNU Lesser General Public
License version 3** (LGPL-3.0-or-later). This license is inherited from
the original GLODS framework (Custódio & Madeira, 2015) on which the
embedded problem definitions and the algorithmic context of this suite
are based.

See [`LICENSE`](LICENSE) for the full license text.

---

## Citation

If you use this benchmark suite in academic work, please cite the
accompanying paper:

```bibtex
@article{Madeira2026GLODSSI,
  author  = {Madeira, J. F. A.},
  title   = {{GLODS-SI}: Scale-Invariant {Global--Local} Direct Search
             for Engineering Design Optimization},
  journal = {Journal of Computational Design and Engineering},
  year    = {2026},
  note    = {Manuscript ID JCDE-2026-065}
}
```

A `CITATION.cff` file is also provided for tools that consume that
metadata format (GitHub, Zenodo, etc.).

### Underlying framework

The original GLODS framework on which this suite is methodologically based:

> Custódio, A. L., Madeira, J. F. A. (2015).
> *GLODS: Global and Local Optimization using Direct Search.*
> Journal of Global Optimization, 62, 1–28.
> [doi:10.1007/s10898-014-0224-9](https://doi.org/10.1007/s10898-014-0224-9)

---

## Acknowledgments

This work was supported by *Fundação para a Ciência e a Tecnologia*
(FCT) through LAETA (project
[UID/50022/2025](https://doi.org/10.54499/UID/50022/2025)).

---

## Contact

**José F. Aguilar Madeira**
IDMEC, Instituto Superior Técnico, Universidade de Lisboa
ISEL, Instituto Politécnico de Lisboa
Email: aguilarmadeira@tecnico.ulisboa.pt
ORCID: [0000-0001-9523-3808](https://orcid.org/0000-0001-9523-3808)
