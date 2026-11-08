# DFO Benchmark Suite — changes

## 2.1.1

- Python version of the 841 problems (`python/core/`): one self-contained file
  per instance, standard library only, same interface as the MATLAB wrappers
  (level indices 1..m; `lambda_` in Family BK). The construction data were
  read from the MATLAB wrappers; the 42 original functions were translated.
- Checked against MATLAB at 26 071 points (x0, alternative starting points,
  x*, random points): 93.1% identical to the last bit, largest relative
  difference 2.0e-14; F(x*) bit-for-bit equal in 828 of 841 instances (the
  other 13: gulf, f* ≈ 1e-30). Check repeatable with `python/verify/`.
- MATLAB wrappers: `func_F` / `cell2idx` now reject a non-finite value for
  an ordered discrete variable (in 2.1.0 a NaN was mapped to the first grid
  value). One line per wrapper; the objective, domains, starting points,
  optima and targets are unchanged, so results on 2.1.0 and 2.1.1 are
  comparable. The Python files apply the same rule.

## 2.1.0

- Decentred domains: the classical box of each of the 55 base instances is
  cropped asymmetrically (deterministic Halton rule, 13–20% of each width
  removed on the side that keeps the minimiser inside; each interval keeps
  80–87% of its classical width), so that the centre
  of the domain is at least 7% of the width away from the optimum in every
  coordinate. In 2.0.0 the optimum was the centre of the box for 105 of the
  841 problems; in 2.1.0 for none. Optimum values
  and reference minimisers are unchanged (x* ∈ new box ⊆ classical box);
  known minimisers outside the new box are dropped.
- Everything else unchanged: construction F = f ∘ T, uniform grids,
  cardinalities, starting-point rule, interfaces, Family BK (whose
  continuous optimum is never at the centre by construction).
- Continuous bounds are no longer round numbers: wrappers accept a value
  within 1e-12 of the width outside a bound and clamp it.
- Results on 2.0.0 and 2.1.0 are not comparable; 2.0.0 stays available
  (tag v2.0.0, doi:10.5281/zenodo.24037856).

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
