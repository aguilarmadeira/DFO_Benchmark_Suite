# DFO Benchmark Suite — changes

## 2.0.0

New collection (the v1 folders are unchanged and remain available under tag `v1.0.0`).

- Known global optimum for every problem: value f*, an admissible optimal
  solution (in `reference_optima.csv`) and a justification class E / A / L
  (`docs/OPTIMA.md`). Base instances whose optimum cannot be justified are
  not included (55 of the 63 v1 base instances are).
- Construction F = f ∘ T on uniform grids of odd cardinality fixed by the
  number of finite variables; F attains f* exactly; the visible values carry
  no information about the location of the optimum.
- Seven type patterns (C, D, K, DK, CD, CK, CDK) and three scale strategies
  (baseline, moderate, extreme) where the strategy changes the problem.
- Family BK: synthetic problems with nominal categorical variables and a
  controlled interaction between the continuous and categorical parts.
- Protocol starting point in every file (`PB.x0_idx`, Halton points) and nine
  alternatives; no instance has its starting point at an optimum.
- Index interface `PB.func_idx(z)` in addition to `PB.func_F(xr)`.
- Fully deterministic: no pseudo-random generator; identical files from
  MATLAB and Octave.
- Evaluation protocol and tools in `benchmark/`; versioned targets in
  `targets/`; per-folder `INDEX.csv`; held-out part.

## 1.0.0

Published version used in the GLODS-SI and GLODS-SI-Mix studies
(DOI 10.5281/zenodo.20782893).
