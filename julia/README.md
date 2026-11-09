# DFO Benchmark Suite v2.1.2 — Julia version

The 841 problems of `core/` (MATLAB) as self-contained Julia files: one file per instance, one module per file,
Base only, no data files. The mathematical instances are those of version 2.1.0 (same domains, starting points and
optima); values agree with MATLAB within the tolerance documented below, so results obtained in Julia, Python and
MATLAB are comparable.

```julia
include("julia/core/CDK_baseline/ackley_cdk_10D.jl")
PB = ackley_cdk_10D.problem()
f  = PB.func_idx(PB.x0_idx)        # protocol starting point
```

`example.jl` runs the quasi-random (Halton) search baseline of `benchmark/solver_random_search.m` with the protocol
of `benchmark/PROTOCOL.md` (it gives the same result as `python/example.py`).

## Interface

`problem()` returns a named tuple with the fields of the MATLAB structure: `nvar`, `var_types` (`"C"`, `"D"`,
`"K"`), `cardinality` (`Inf` for continuous), `lb_C_work`, `ub_C_work`, `ProblemData`, `x0_idx`, `x0_alt_idx`,
`f_star`, `f_star_justification`, and the functions `func_idx(z)` (z = vector of **1-based level indices**,
continuous coordinates as values, as in MATLAB), `func_F(xr)` (values, grid values, labels), `idx2cell`,
`cell2idx`. Invalid input (wrong length, non-finite or non-numeric values, indices outside 1..m, unknown labels,
values off the grid) throws an `ArgumentError`. Optimal solutions: `core/<folder>/reference_optima.csv` of the
MATLAB part.

When the files are loaded inside a function (not at top level), call the returned functions through
`Base.invokelatest`, as `verify/eval_points.jl` does (world-age rule of Julia ≥ 1.12).

## How the files were made and checked

Made like the Python version: the 42 original functions translated by hand keeping the order of operations, the
construction data read from the MATLAB wrappers and written as literals (no second generator). Checked with Julia
1.13 at the 26 071 points of `python/verify` (x0, the 9 alternative starting points, x* and 20 random points of
every problem): 93.5% of the values identical to MATLAB to the last bit, largest relative difference 2.0e-14, none
above 1e-12 (97.9% identical to Python); F(x*) bit-for-bit equal to MATLAB in 828 of 841 problems (the other 13:
gulf, f* ≈ 1e-30). Invalid-input test: 841 problems, 33 684 checks, 0 failures.

To repeat: `julia verify/compare_values.jl --julia` (exit code 1 on any difference above 1e-12) and
`julia verify/test_inputs.jl core`.
