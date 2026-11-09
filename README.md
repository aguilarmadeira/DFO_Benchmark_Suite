# DFO Benchmark Suite — version 2.1.2

[![License: LGPL v3](https://img.shields.io/badge/License-LGPL_v3-blue.svg)](https://www.gnu.org/licenses/lgpl-3.0)
[![DOI (all versions)](https://img.shields.io/badge/DOI_(all_versions)-10.5281%2Fzenodo.20782892-blue.svg)](https://doi.org/10.5281/zenodo.20782892)

A reproducible collection of deterministic black-box optimisation problems
with continuous (C), ordered discrete (D) and categorical (K) variables,
independent of any solver. Every problem has a **known global optimum** with a
stated justification, and the evaluation protocol uses only the points and
values a solver evaluates. The domains are chosen so that the centre of the
box is never the optimum (no advantage for methods that evaluate the centre
first).

## Three languages, same 841 problems

Every problem is a single self-contained file, with no dependencies, in each
of three languages. The three versions give the same values (checked at
26 071 points; largest relative difference 2·10⁻¹⁴).

| language | folder | requires | details |
|---|---|---|---|
| **MATLAB / Octave** | `core/` | nothing | `docs/CORE_SPEC.md` |
| **Python** | `python/core/` | standard library only | [`python/README.md`](python/README.md) |
| **Julia** | `julia/core/` | Base only | [`julia/README.md`](julia/README.md) |

The folder names inside the three (`C_baseline`, `CDK_baseline`, `BK_binary`,
...) and the instance names are the same. The optima (`core/<folder>/reference_optima.csv`),
the instance index (`core/<folder>/INDEX.csv`) and the targets (`targets/`) are
plain CSV files, shared by the three versions.

**MATLAB / Octave**

```matlab
addpath('core/CDK_baseline');
PB = ackley_cdk_10D();                 % one problem
f  = PB.func_idx(PB.x0_idx);           % evaluate at the protocol's starting point
PB.var_types, PB.cardinality, PB.f_star
```

**Python**

```python
import sys; sys.path.insert(0, 'python/core/CDK_baseline')
from ackley_cdk_10D import ackley_cdk_10D
PB = ackley_cdk_10D()
f  = PB.func_idx(PB.x0_idx)
```

**Julia**

```julia
include("julia/core/CDK_baseline/ackley_cdk_10D.jl")
PB = ackley_cdk_10D.problem()
f  = PB.func_idx(PB.x0_idx)
```

The benchmarking tools (`benchmark/`: run a solver, data and performance
profiles) are in MATLAB:

```matlab
addpath('benchmark');
run_benchmark(solver_random_search(), 'core/D');   % 1000 n evaluations, 10 seeds
data_profiles({'random_search'}, 'core/D');
```

To benchmark your own solver in MATLAB write a struct `solver.name`, `solver.fun =
@(P) ...` (see `benchmark/solver_random_search.m` and `benchmark/PROTOCOL.md`).
The solver receives the domain, the objective and the starting point only.
In Python or Julia, any solver can call `PB.func_idx` directly; the protocol
(`benchmark/PROTOCOL.md`) says what to record.

## Contents

| path | what |
|---|---|
| `core/` | MATLAB/Octave version: 17 folders, 841 problems (below) |
| `core/<folder>/reference_optima.csv` | optimal solution of every instance (not used by solvers) |
| `core/<folder>/INDEX.csv` | per instance: dimension, cardinalities, f*, enumeration certificate, ... |
| `targets/` | versioned targets of the protocol (one file per folder) |
| `benchmark/` | protocol (`PROTOCOL.md`) and tools: run, validate, targets, profiles |
| `python/` | Python version: same 17 folders, example and MATLAB–Python check |
| `julia/` | Julia version: same 17 folders, example and MATLAB–Julia check |
| `docs/CORE_SPEC.md` | specification: construction, interface, starting points |
| `docs/OPTIMA.md` | the 55 base instances, their optimum and its justification |
| `LICENSE`, `CITATION.cff`, `CHANGELOG.md` | licence (LGPL-3.0-or-later), citation metadata, changes |

| folders | variable types | problems |
|---|---|---|
| `C_baseline`, `C_moderate`, `C_extreme` | continuous | 55 + 55 + 55 |
| `D`, `K`, `DK` | ordered discrete, categorical, alternating | 55 each |
| `CD_*`, `CK_*`, `CDK_*` (baseline / moderate / extreme) | mixed | 55 / 35 / 35, 55 / 35 / 35, 55 / 31 / 31 |
| `BK_MK16`, `BK_binary` | continuous + nominal categories (synthetic, controlled interaction) | 54 + 90 |

The 697 problems of the first rows derive from 55 classical base instances
(dimensions 2 to 10); finite coordinates take values on uniform grids and the
objective is F = f ∘ T, with T a deformation of at most half a grid cell that
keeps the global optimum value f* exactly (`docs/CORE_SPEC.md`). The moderate
and extreme folders rescale the continuous coordinates and contain only the
instances with at least two continuous variables.

## Verification

841 instances, 0 problems, 0 trivial starting points (`benchmark/finalize_core`);
F(x*) = f* in every file; generation reproduced identically in MATLAB and
Octave; optimum of the fully finite versions certified by enumeration where
the grid has at most 2·10⁶ points (39 of 55 instances in each of D, K, DK).
Python version against MATLAB at 26 071 points: 93.1% of the values identical,
largest relative difference 2.0·10⁻¹⁴ (`python/README.md`); Julia version:
93.5% identical, largest relative difference 2.0·10⁻¹⁴ (`julia/README.md`).

## Held-out part

`benchmark/holdout_instances.m` lists 19 base instances whose results are to
be reported separately and never used to tune a solver.

## Versions

**Version 2.1.2** adds the Julia version (`julia/`). **Version 2.1.1** added
the Python version (`python/`) and made the MATLAB
wrappers reject a non-finite value for an ordered discrete variable; the
problems, starting points, optima and targets are those of version 2.1.0
(doi:10.5281/zenodo.24041183), and results on 2.1.0, 2.1.1 and 2.1.2 are comparable.

**Version 2.0.0** (tag [`v2.0.0`](https://github.com/aguilarmadeira/DFO_Benchmark_Suite/tree/v2.0.0),
doi:10.5281/zenodo.24037856) is superseded: it had the same construction on
the classical boxes, where the optimum is the centre of the box for 105 of
the 841 problems. Version 2.1.0 crops each classical box asymmetrically
(`docs/CORE_SPEC.md`); optimum values and minimisers are unchanged.

**Version 1.0.0** — the continuous-variable suite (504 scaled wrappers and
the scaling toolkit) and the mixed-variable suite used in the GLODS-SI
(JCDE, doi:10.1093/jcde/qwag049) and GLODS-SI-Mix studies — remains
available, unchanged, under the tag
[`v1.0.0`](https://github.com/aguilarmadeira/DFO_Benchmark_Suite/tree/v1.0.0)
and on Zenodo (doi:10.5281/zenodo.20782893). Results obtained on different
versions are not comparable.

## License

Distributed under the **GNU Lesser General Public License version 3**
(LGPL-3.0-or-later), as version 1, inherited from the GLODS framework
(Custódio & Madeira, 2015) from which the embedded original test-function
definitions come. See [`LICENSE`](LICENSE).

## Citation

```bibtex
@software{madeira_dfo_benchmark_suite_v2,
  author    = {Madeira, J. F. A.},
  title     = {{DFO Benchmark Suite}},
  year      = {2026},
  version   = {2.1.2},
  publisher = {Zenodo},
  doi       = {10.5281/zenodo.20782892},
  url       = {https://doi.org/10.5281/zenodo.20782892}
}
```

The DOI 10.5281/zenodo.20782892 cites all versions and resolves to the latest;
each version has its own DOI on Zenodo (v2.1.1: 10.5281/zenodo.24055130;
v2.1.0: 10.5281/zenodo.24041183;
v2.0.0: 10.5281/zenodo.24037856;
v1.0.0: 10.5281/zenodo.20782893).
`CITATION.cff` carries the same metadata for GitHub and Zenodo.

> Custódio, A. L., Madeira, J. F. A. (2015). *GLODS: Global and Local
> Optimization using Direct Search.* Journal of Global Optimization, 62,
> 1–28. [doi:10.1007/s10898-014-0224-9](https://doi.org/10.1007/s10898-014-0224-9)

## Acknowledgments

This work was supported by *Fundação para a Ciência e a Tecnologia* (FCT)
through LAETA (project [UID/50022/2025](https://doi.org/10.54499/UID/50022/2025)).

## Contact

**J. F. A. Madeira** — IDMEC, Instituto Superior Técnico, Universidade de
Lisboa; ISEL, Instituto Politécnico de Lisboa.
Email: <aguilarmadeira@tecnico.ulisboa.pt> ·
ORCID: [0000-0001-9523-3808](https://orcid.org/0000-0001-9523-3808)
