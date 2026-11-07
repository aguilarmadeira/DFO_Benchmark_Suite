# DFO Benchmark Suite — version 2.1.0 core: public specification

Version 2.1.0 differs from 2.0.0 only in the domains of the 697 classical
problems (section 2, *Domain*); the construction, grids, cardinalities,
interfaces, starting-point rule and Family BK are unchanged. Results obtained
on 2.0.0 and 2.1.0 are not comparable.

A reproducible collection of black-box optimisation problems with continuous
(C), ordered discrete (D) and categorical (K) variables, independent of any
solver, with a known and justified global optimum for every instance.

## 1. Contents

| folder | pattern of variable types | instances | scale strategy |
|---|---|---|---|
| `C_baseline`, `C_moderate`, `C_extreme` | all continuous | 55 each | baseline / moderate / extreme |
| `D` | all ordered discrete | 55 | — |
| `K` | all categorical | 55 | — |
| `DK` | D, K, D, K, … | 55 | — |
| `CD_baseline`, `CD_moderate`, `CD_extreme` | C, D, C, D, … | 55 / 35 / 35 | baseline / moderate / extreme |
| `CK_baseline`, `CK_moderate`, `CK_extreme` | C, K, C, K, … | 55 / 35 / 35 | idem |
| `CDK_baseline`, `CDK_moderate`, `CDK_extreme` | C, K, D, C, D, K, … | 55 / 31 / 31 | idem |
| `BK_MK16`, `BK_binary` | C + nominal categories (synthetic, Family BK) | 54 + 90 | — |

The 55 base instances are classical test functions (54 functions, dimensions 2
to 10; `reference_optima.csv` of every folder lists them). The moderate and
extreme strategies rescale the continuous coordinates (contrast 10^4 and
10^8 between coordinates) and are generated only for instances with at least
two continuous variables, where they change the problem.

## 2. Construction

* **Domain (v2.1.0).** Many classical functions have their optimum at the
  centre of their classical box Ω₀ (14 of the 55 base instances: ackley,
  griewank, rastrigin, sphere, …), which favours any method that evaluates
  the centre. Each classical box is therefore cropped asymmetrically to a
  sub-box Ω ⊆ Ω₀ that contains x*: since x* ∈ Ω ⊆ Ω₀, min over Ω of f = f* at
  the same x*, so the optimum value and the reference minimiser are
  unchanged. The crop uses the location of a known minimiser: the new boxes
  are a construction choice and are not independent of x*. Rule, for base
  instance p (index in the original list of 63) and coordinate i, with
  H = Halton point 7000 + p in dimension 2n: ρ_i = 0.13 + 0.07 H_i ∈
  [0.13, 0.20] is the fraction of the width removed (0.13 is the smallest
  fraction that leaves an optimum at the classical centre 7% of the new
  width away: ρ/(2(1−ρ)) = 0.0747); the candidates are the upper crop
  [a_i, b_i − ρ_i w_i] and the lower crop [a_i + ρ_i w_i, b_i]; a candidate is
  admissible if it contains x*_i and its centre is at least 7% of its width
  away from x*_i; if both are admissible the side is H_{n+i} < 0.5 → upper;
  if neither, ρ_i is adjusted in steps of 0.01 within [0.13, 0.20] (never
  needed for the 55 instances). Result: each interval keeps 80–87% of its
  classical width (volume kept, median: 69% for n = 2, 43% for n = 5, 16% for
  n = 10); in every coordinate of every classical instance the centre of the
  domain is at least 7% of the width from x* (minimum 7.2%), so x* is never
  placed at the central grid index (half a cell is at most 6.25% of the width,
  m ≥ 9); every other known global minimiser that remains in Ω is at least
  9.0% (∞-norm) from the centre (for shubert one coordinate of one of them is
  close to the central value, the other is not). Known minimisers outside Ω
  are dropped (branin 3 → 2, cosine_mixture 2ⁿ → 1, dekkers_aarts 2 → 1,
  shubert 18 → 12). The domain is the same for all patterns and strategies of
  a base instance. Check on the 841 instances: the centre never reaches the
  convergence level τ = 1e−7 (v2.0.0: 105 instances). Axial points (centre
  with one coordinate at a bound) are only a diagnostic: none is optimal; one
  reaches τ = 1e−7 relative to f(x0) — miele_cantrell K 4D, F = 24.3 against
  f* = 0, because f(x0) ≈ 1.1·10¹⁰ for that function (a property of the
  relative test, see benchmark/PROTOCOL.md); v2.0.0: 13 instances.

For a base function f on the box Ω = Π[a_i, b_i] (the domain above) with a
known global minimiser x* and value f* = f(x*):

* **Continuous coordinates** keep their interval (baseline) or an affinely
  rescaled one (moderate, extreme); the wrapper maps the work interval back
  to [a_i, b_i] affinely.
* **Finite coordinates** (D and K) take values on the uniform grid of m levels
  `g_k = a_i + (k−1)(b_i − a_i)/(m−1)`. m is odd and depends only on the
  number of finite variables n_f of the instance: 1 → 65, 2 → 151, 3 → 41,
  4 → 25, 5–6 → 17, 7–8 → 13, 9–10 → 9. Discrete variables expose the grid
  values; categorical variables expose nominal labels `k<i>_1 … k<i>_m`,
  each carrying one grid value through a fixed permutation (the labels have
  no order and no distance; the value behind a label is not exposed).
* **Objective** `F(z) = f(T(z))`. On each finite coordinate, `z*_i` is the grid
  point nearest to `x*_i`, index `k = floor(t + 1/2) + 1` with
  `t = (x*_i − a_i)(m − 1)/(b_i − a_i)` (an exact tie, e.g. x* = −1 on [−2, 2]
  with m = 151, goes to the upper index; the rule uses only arithmetic
  operations, so it gives the same index on every platform), and T_i is the piecewise-affine bijection of
  [a_i, b_i] with `T_i(z*_i) = x*_i` that fixes the end points (if `z*_i` is
  an end point, T_i maps [a_i, b_i] onto [x*_i, b_i] or [a_i, x*_i]). On the
  continuous coordinates T is the identity. Hence
  `‖T_i − Id‖∞ = |x*_i − z*_i| ≤ (b_i − a_i)/(2(m−1))`, `T(domain) ⊆ Ω`, and
  `min F = f* = F(z*)` exactly. F is f up to a deformation of at most half a
  grid cell; it is not f restricted to a grid (which would in general not
  attain f*).
* **Optimum.** `PB.f_star` is f* and `PB.f_star_justification` says how it is
  justified: `E` exact by structure (a sum of non-negative terms attaining 0,
  or an identity giving the bound), `A` analytic reduction (separable or
  symmetric function, 1-D root computed to machine precision), `L` literature
  (classical multimodal functions whose optimum is the value established in
  the literature, at a refined minimiser). Of the 55 base instances, 36 are
  E, 6 are A and 13 are L (goldstein_price, hartman, hosaki, kowalik,
  mccormick, multi_gaussian, shekel_45/47/410, shekel_foxholes ×2, shubert);
  the L optima are not proved and are reported as such (11 of the 13 are
  certified by enumeration on the fully finite versions; `docs/OPTIMA.md`).
  `f_star` is the value computed by the function at x* (e.g. 4.44e−16 for
  ackley at the origin), so that F(z*) = f_star holds exactly in floating
  point. The optimal solution is in `reference_optima.csv` (index
  representation and original coordinates), not in the wrapper's struct `PB`:
  the construction data (grid breakpoints of T, values carried by the labels)
  are kept in a private struct captured by the function handles. The wrapper
  is self-contained source code, so a reader of the file can of course
  reconstruct the optimum; what matters for the protocol is that nothing in
  `PB` reveals it and that nothing is printed.
* **Starting point.** `PB.x0_idx` is the Halton point of index 1000 + p (p =
  index of the base instance; Family BK: 2000 or 3000 + 10 st + seed, see
  section 4), scaled to the box and rounded to the grid; it is chosen without
  reference to the function. Nine further points
  (`PB.x0_alt_idx`, Halton indices 1000 + p + 63k) serve multi-start
  protocols. An instance whose x0 happens to be optimal is reported as
  trivial by the validation, never altered.
* **Determinism.** Grids, permutations, starting points and Family BK are
  produced by deterministic rules and the Halton sequence (radical inverse in
  the first n primes); no pseudo-random generator is used anywhere. The same
  files are obtained in MATLAB and Octave up to last-digit differences of the
  elementary functions (file-by-file comparison, tolerance 1e−12;
  the only visible effect is `f_star` of gulf, ≈ 1e−30 with the true value 0,
  whose last digits depend on the exp of the platform). Verified on
  2026-11-07: 0 differing files between a MATLAB and an Octave generation.

## 3. Interface of a wrapper

```matlab
PB = ackley_cdk_10D();          % prints one line; PB is a struct
f  = PB.func_idx(z);            % z: 1 x n, continuous coordinates as values, finite ones as 1-based level indices
f  = PB.func_F(xr);             % xr: 1 x n cell of values (C), grid values (D) and labels (K)
PB.var_types, PB.cardinality    % {'C','K','D',...}, levels per variable (Inf for continuous)
PB.ProblemData                  % [lb;ub] for C, {grid values} for D, {labels} for K
PB.lb_C_work, PB.ub_C_work      % bounds of the continuous coordinates
PB.x0_idx, PB.x0_alt_idx        % starting point(s) of the protocol
PB.f_star, PB.f_star_justification
PB.idx2cell(z), PB.cell2idx(xr) % conversions
PB.base_function, PB.pattern, PB.strategy, PB.m, PB.n_finite, PB.M_K   % descriptive
```

Every wrapper is self-contained (the original function is embedded) and
needs nothing but MATLAB or Octave.

## 4. Family BK (synthetic, nominal categories)

`F(x_C, κ) = ‖x_C − x_C* − λ s(κ)‖² + b (1 − h(κ)/n_K)`, with n_K nominal
variables of m categories, h(κ) the number of coordinates of κ equal to the
optimal combination, s(κ) = s only for the combinations with h = 0. The
optimum is (x_C*, κ*), f* = 0, by construction (`PB.f_star_justification =
'E'`; the solution is in `reference_optima.csv`); the parameter q = λ²σ²/b
(`PB.obstacle_ratio`) and `PB.j_min` (minimal number of simultaneous
categorical changes that improves from an antipodal conditional minimiser
without continuous re-optimisation) describe the interaction. `BK_MK16` fixes
the number of combinations at 16 with 1 × 16, 2 × 4 and 4 × 2 categories;
`BK_binary` uses 2, 4, 6 and 8 binary variables; nC = 2 continuous
variables in [0, 1], σ = 0.15, λ ∈ {0, 1}, seeds 1–3, b chosen so that j_min takes the
values 1, ⌈nK/2⌉, nK and ∞ at λ = 1. Starting point: Halton index
2000 + 10 st + seed (`BK_MK16`) or 3000 + 10 st + seed (`BK_binary`), st =
position of the structure in the design list; instances that differ only in
λ or b share x*, the labels and x0, so that the effect of the interaction
parameters is isolated. Because x0 is chosen without reference to the function, its
categorical part coincides with the optimal combination in 10 of the 144
instances (`BK_MK16`: the 4 variants of 1 × 16, seed 2; `BK_binary`: the 6
variants of 2 × 2, seed 1); they are not trivial (the continuous part of x0 is
not optimal) and are kept, but in them the categorical search plays no role
from x0. `analyze_interact` reports h(x0) for every run (column `h_x0`).

## 5. Protocol

`benchmark/PROTOCOL.md`: budget 1000 n, start at `PB.x0_idx`, Moré–Wild
convergence test with the targets of `targets/`, invalid evaluations count
and return +Inf, data and performance profiles, path statistics, held-out
part (`benchmark/holdout_instances.m`). Results on different versions of the
suite are not comparable; version 1 (published) and version 2 are both
frozen once released.

## 6. What is not in the core

The complementary research collections (controlled cardinalities and
proportions, label permutations, further interaction families, diagnostic
designs) are not part of version 2. The generators are not distributed;
the published files are the reference.

## 7. Verification (2026-11-07)

841 instances, 0 problems, 0 trivial (`benchmark/validate_suite`, MATLAB and
Octave); F(x*) = f* in every file; every evaluated starting point and sample
point admissible; MATLAB and Octave generations identical up to 1e−12
relative (0 differing files); enumeration of the finite domain (≤ 2·10⁶
points) in 39 of 55 instances of each of D, K and DK, all agreeing with f*.
