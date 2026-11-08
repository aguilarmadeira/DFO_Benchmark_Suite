# DFO Benchmark Suite v2.1.1 — Python version

The 841 problems of `core/` (MATLAB) as self-contained Python files: one file per instance, standard library only
(`math`), no data files. The mathematical instances are those of version 2.1.0: same domains, starting points and
optima. Values agree with MATLAB within the tolerance documented below (not necessarily to the last bit at every
point), so results obtained in Python and in MATLAB are comparable.

```python
import sys; sys.path.insert(0, 'python/core/CDK_baseline')
from ackley_cdk_10D import ackley_cdk_10D
PB = ackley_cdk_10D()
f = PB.func_idx(PB.x0_idx)        # protocol starting point
```

`example.py` runs the quasi-random (Halton) search baseline of `benchmark/solver_random_search.m` with the
protocol of `benchmark/PROTOCOL.md`.

## Interface

| field | meaning |
|---|---|
| `PB.nvar`, `PB.var_types` | dimension and types (`'C'`, `'D'`, `'K'`) |
| `PB.cardinality` | number of levels per variable (`math.inf` for continuous) |
| `PB.lb_C_work`, `PB.ub_C_work` | bounds of the continuous variables |
| `PB.func_idx(z)` | objective; `z` = list of **1-based level indices** (continuous coordinates as values), as in MATLAB |
| `PB.func_F(xr)` | objective; `xr` = list of values / grid values / labels |
| `PB.idx2cell`, `PB.cell2idx` | conversions |
| `PB.x0_idx`, `PB.x0_alt_idx` | protocol starting point and the 9 alternatives |
| `PB.f_star`, `PB.f_star_justification` | known global optimum and its justification (E / A / L) |
| `PB.ProblemData` | domain: `[lb, ub]` (C), `[grid]` (D), list of labels (K) |

Family BK: `PB.lambda_` (`lambda` is a Python keyword). Optimal solutions are in
`core/<folder>/reference_optima.csv` of the MATLAB part (not used by solvers). Invalid input (wrong length,
non-finite or non-numeric values, indices outside 1..m, unknown labels, values off the grid) raises `ValueError`.

## How the files were made and checked

The 42 original test functions were translated by hand, keeping the order of operations of the MATLAB code. The
construction data (domains, grids, permutations of the categorical levels, breakpoints of T, anchors, starting
points) were read from the MATLAB wrappers of v2.1.0 and written as literals; there is no second generator.

Comparison with MATLAB at 26 071 points (x0, the 9 alternative starting points, x* and 20 random points of every
instance; `verify/`): 93.1% of the values identical to the last bit, largest relative difference 2.0e-14, none above
1e-12. F(x*) is bit-for-bit equal to MATLAB in 828 of the 841 instances; the other 13 are gulf, where
f* ≈ 1e-30 (rounding noise around zero). To repeat the check: run `verify/eval_points.m` in MATLAB and then
`python verify/compare_values.py` (exit code 1 if the points differ or any relative difference exceeds 1e-12).
