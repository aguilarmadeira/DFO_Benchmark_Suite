# DFO Benchmark Suite

A benchmark suite for **Derivative-Free Optimization (DFO)**, organized by
the type of decision variables. The current release provides two
categories: bound-constrained **continuous-variable** problems under
controlled scale heterogeneity, and **mixed-variable** problems combining
continuous, ordered discrete, and categorical variables.

[![License: LGPL v3](https://img.shields.io/badge/License-LGPL_v3-blue.svg)](https://www.gnu.org/licenses/lgpl-3.0)
[![Paper DOI](https://img.shields.io/badge/Paper-10.1093%2Fjcde%2Fqwag049-blue.svg)](https://doi.org/10.1093/jcde/qwag049)

---

## Repository structure

```text
DFO_Benchmark_Suite/
├── README.md                   (this file)
├── LICENSE                     (GNU LGPL v3)
├── CITATION.cff                (citation metadata)
│
├── prob_var_continuous/        Category: continuous variables
│   ├── README.md
│   ├── originals/              48 original problem definitions
│   ├── scaled/                 504 self-contained scaled wrappers (63 × 8)
│   └── scaling_toolkit/        3 MATLAB functions to scale new problems
│
└── prob_var_mixed/             Category: mixed variables (C / D / K)
    ├── README.md
    ├── problemsMix/            504 self-contained mixed wrappers (63 × 8)
    └── metadata/               machine-readable inventory of instances
```

---

## Current release

The repository contains two benchmark categories.

### Continuous-variable benchmark suite

The `prob_var_continuous/` category contains bound-constrained continuous
problems under controlled scale heterogeneity. It includes the original
continuous problem definitions, self-contained scaled wrappers, and the
scaling toolkit used to generate heterogeneous continuous test instances.
See [`prob_var_continuous/README.md`](prob_var_continuous/README.md).

### Mixed-variable benchmark suite

The `prob_var_mixed/` category contains the mixed-variable benchmark
instances used in the GLODS-SI-Mix computational study. These problems
combine continuous, ordered discrete, and categorical variables. The main
synthetic suite contains **504 self-contained MATLAB wrappers** generated
from **63 classical continuous test-function instances** under **8
deterministic heterogeneity strategies**.

Each mixed wrapper is self-contained at runtime and exposes a problem
structure with variable types, domains, decoding rules, categorical
embeddings, and an objective handle. The mixed-variable construction
follows a fixed cyclic variable pattern of continuous, categorical, and
ordered discrete variables, while categorical variables are represented
through permutation-invariant labels and evaluated through embedded
numerical grids.

The category also reports the **16-problem Cat-Suite** benchmark used for
comparison with CatMADS and published baselines; those baseline values are
taken from the publicly available CatMADS data
(`github.com/bbopt/CatMADS_prototype`) and are not redistributed here. See
[`prob_var_mixed/README.md`](prob_var_mixed/README.md).

---

## Quick start

Continuous wrappers expose a three-mode interface (info / bounds / value):

```matlab
addpath(genpath('prob_var_continuous/scaled'));

info     = ackley_10D();        % metadata struct
[lb, ub] = ackley_10D(10);      % work-space bounds
f        = ackley_10D(rand(10,1));
```

Mixed wrappers return a problem struct and evaluate a cell array of values:

```matlab
addpath(genpath('prob_var_mixed/problemsMix'));

PB = ackley_k83_k47_10D();      % problem struct
f  = PB.func_F(PB.x_star);      % objective at a representative optimizer
```

---

## License

This repository is distributed under the **GNU Lesser General Public
License version 3** (LGPL-3.0-or-later), inherited from the original GLODS
framework (Custódio & Madeira, 2015) on which the embedded problem
definitions and the algorithmic context of this suite are based. See
[`LICENSE`](LICENSE) for the full text.

---

## Citation

If you use this benchmark suite in academic work, please cite the
accompanying papers.

Archived suite (software/dataset):

```bibtex
@dataset{madeira_dfo_benchmark_suite_2026,
  author    = {Madeira, J. F. A.},
  title     = {{DFO Benchmark Suite: Continuous and Mixed-Variable Test Problems for Derivative-Free Optimization}},
  year      = {2026},
  version   = {1.0.0},
  publisher = {Zenodo},
}
```

Continuous-variable suite:

```bibtex
@article{Madeira2026GLODSSI,
  author  = {Madeira, J. F. A.},
  title   = {{GLODS-SI}: Scale-Invariant {Global--Local} Direct Search
             for Engineering Design Optimization},
  journal = {Journal of Computational Design and Engineering},
  year    = {2026},
  note    = {qwag049},
  doi     = {10.1093/jcde/qwag049}
}
```

Mixed-variable suite (GLODS-SI-Mix):

```bibtex
@unpublished{MadeiraGLODSSIMix,
  author = {Madeira, J. F. A.},
  title  = {{GLODS-SI-Mix}: Scale-Invariant Direct Search for
            Mixed-Variable Derivative-Free Optimization},
  note   = {Manuscript under review}
}
```

A `CITATION.cff` file is also provided for tools that consume that
metadata format (GitHub, Zenodo, etc.).

### Underlying framework

> Custódio, A. L., Madeira, J. F. A. (2015). *GLODS: Global and Local
> Optimization using Direct Search.* Journal of Global Optimization, 62,
> 1–28. [doi:10.1007/s10898-014-0224-9](https://doi.org/10.1007/s10898-014-0224-9)

The Cat-Suite problems used in the mixed-variable comparison are due to the
CatMADS authors (`github.com/bbopt/CatMADS_prototype`; report G-2025-39).

---

## Acknowledgments

This work was supported by *Fundação para a Ciência e a Tecnologia* (FCT)
through LAETA (project [UID/50022/2025](https://doi.org/10.54499/UID/50022/2025)).

---

## Contact

**J. F. A. Madeira** — IDMEC, Instituto Superior Técnico,
Universidade de Lisboa; ISEL, Instituto Politécnico de Lisboa.
Email: <aguilarmadeira@tecnico.ulisboa.pt> ·
ORCID: [0000-0001-9523-3808](https://orcid.org/0000-0001-9523-3808)
