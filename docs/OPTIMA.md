# Global optima of the 55 base instances (DFO Benchmark Suite v2.1.0 core)

Every instance of the core derives from one of these base instances and has the same global optimum value f*
(for F = f o T, F(z*) = f* exactly). Justification classes: **E** exact by structure (a sum of non-negative
terms attaining its bound, or an identity giving the bound); **A** analytic reduction (separable or symmetric
function reduced to a 1-D problem solved to machine precision); **L** the value established in the literature,
at a minimiser refined numerically (not proved). The domain of every instance is a sub-box of the classical box
that contains the minimiser (docs/CORE_SPEC.md), so f* is unchanged; *global minimisers* counts the known ones
that lie in that domain. Column *enumerated*: the optimum of the fully finite versions (folders D, K, DK) was
certified by exhaustive enumeration of the grid (domains of at most 2e6 points); this certifies min F on the
grid and is evidence, not proof, for the continuous optimum of f. The optimal solutions are in
`reference_optima.csv` of each folder. The justification notes describe the function on its classical box.

| # | base instance | n | class | f* | global minimisers | enumerated (D/K/DK) | justification |
|---|---|---|---|---|---|---|---|
| 1 | ackley | 10 | E | 4.440892099e-16 | 1 | no | f >= 0, = 0 only at 0 (both terms vanish only at 0) |
| 2 | aluffi_pentini | 2 | A | -0.3523860738 | 1 | yes | separable: 0.25x^4-0.5x^2+0.1x + 0.5y^2; 1-D quartic, root refined (-1.04668, 0); globality by 1-D analysis |
| 3 | becker_lago | 2 | E | 0 | 4 | yes | (/x/-5)^2+(/y/-5)^2 = 0 at (±5,±5); all four exact |
| 4 | bohachevsky | 2 | E | 0 | 1 | yes | sum of nonneg terms = 0 only at 0 |
| 5 | branin_hoo | 2 | A | 0.3978873577 | 2 | yes | classical; three minimizers at (-pi,12.275),(pi,2.275),(3pi,2.475), f* = 5/(4pi); globality by the standard analysis of the quadratic-in-x2 form |
| 6 | cosine_mixture | 2 | E | -2.2 | 1 | yes | f = 0.1 sum cos(5 pi x) - sum x^2 on [-1,1]: f >= -1.1 n with equality iff x_i = ±1 (2^n corner minimizers), f* = -2.2 (n=2), -4.4 (n=4) |
| 7 | cosine_mixture | 4 | E | -4.4 | 1 | yes | f = 0.1 sum cos(5 pi x) - sum x^2 on [-1,1]: f >= -1.1 n with equality iff x_i = ±1 (2^n corner minimizers), f* = -2.2 (n=2), -4.4 (n=4) |
| 8 | dekkers_aarts | 2 | A | -24776.51834 | 1 | yes | x1 = 0 by the 1e5 x1^2 term; 1-D in x2: root refined to ±14.94511, f* = -24776.518 (literature -24777 at ±15 is rounded); globality by 1-D analysis after showing x1 = 0 |
| 9 | exponencial | 2 | E | -1 | 1 | yes | -exp(-0.5/x/^2) >= -1, = -1 only at 0 |
| 10 | exponencial | 4 | E | -1 | 1 | yes | -exp(-0.5/x/^2) >= -1, = -1 only at 0 |
| 11 | fifteenn_local_minima | 2 | E | 1.349783804e-32 | 1 | yes | sum of nonneg terms, = 0 only at ones (sin(3 pi) = 0, x_i = 1) |
| 12 | fifteenn_local_minima | 4 | E | 1.349783804e-32 | 1 | yes | sum of nonneg terms, = 0 only at ones (sin(3 pi) = 0, x_i = 1) |
| 13 | fifteenn_local_minima | 6 | E | 1.349783804e-32 | 1 | no | sum of nonneg terms, = 0 only at ones (sin(3 pi) = 0, x_i = 1) |
| 14 | fifteenn_local_minima | 8 | E | 1.349783804e-32 | 1 | no | sum of nonneg terms, = 0 only at ones (sin(3 pi) = 0, x_i = 1) |
| 15 | fifteenn_local_minima | 10 | E | 1.349783804e-32 | 1 | no | sum of nonneg terms, = 0 only at ones (sin(3 pi) = 0, x_i = 1) |
| 16 | fletcher_powell | 3 | E | 0 | 1 | yes | sum of squares of residuals, = 0 at (1,0,0); +Inf on the plane x1 = 0 by its original definition (excluded by the v2.1.0 domain) |
| 17 | goldstein_price | 2 | L | 3 | 1 | yes | f* = 3 at (0,-1): classical, verified in the literature; polynomial, globality by the standard analysis |
| 18 | griewank | 10 | E | 0 | 1 | no | sum x^2/4000 - prod cos + 1 >= 0, = 0 only at 0 |
| 19 | griewank | 5 | E | 0 | 1 | yes | sum x^2/4000 - prod cos + 1 >= 0, = 0 only at 0 |
| 20 | gulf | 3 | E | 9.709124425e-31 | 1 | yes | least squares with zero residual at (50,25,1.5) (exact fit), f* = 0 |
| 21 | hartman_4 | 3 | L | -3.862782148 | 1 | yes | Hartman 3 (-3.86278) and 6 (-3.32237): literature only; refined point; 4 local minima |
| 22 | hartman_4 | 6 | L | -3.322368011 | 1 | no | Hartman 3 (-3.86278) and 6 (-3.32237): literature only; refined point; 4 local minima |
| 23 | hosaki | 2 | L | -2.345811576 | 1 | yes | f* = -2.3458 at (4,2): literature |
| 24 | kowalik | 4 | L | 0.0003074859878 | 1 | yes | least squares with nonzero residual 3.07486e-4; globality from the literature; refined point |
| 25 | mccormick | 2 | L | -1.913222955 | 1 | yes | -1.91322 at (-0.54719,-1.54719): classical; refined |
| 26 | miele_cantrell | 4 | E | 0 | 1 | yes | sum of nonneg terms = 0 at (0,1,1,1) |
| 27 | multi_gaussian | 2 | L | -1.296954046 | 1 | yes | -1.29695 at (-0.01356,-0.01356): literature |
| 28 | neumaier2 | 4 | E | 0 | 1 | yes | sum of squares = 0 at (1,2,2,3) |
| 29 | neumaier3 | 10 | E | -210 | 1 | no | Trid: f* = -n(n+4)(n-1)/6 = -210, x*_i = i(n+1-i), exact |
| 30 | paviani | 10 | A | -45.77845205 | 1 | no | symmetric; reduces to 1-D in the common coordinate 9.351; globality by symmetry + 1-D analysis (written proof pending) |
| 31 | periodic | 2 | E | 0.9 | 1 | yes | 1 + sin^2 x + sin^2 y - 0.1 exp(-x^2-y^2) >= 0.9, = 0.9 only at 0 |
| 32 | powell | 4 | E | 0 | 1 | yes | sum of nonneg terms = 0 at 0 |
| 33 | rastrigin | 10 | E | 0 | 1 | no | sum (x^2 - 10 cos(2 pi x) + 10) >= 0, = 0 only at 0 |
| 34 | rosenbrock | 2 | E | 0 | 1 | yes | = 0 only at ones |
| 35 | salomon | 5 | E | 0 | 1 | yes | 1 - cos(2 pi r) + 0.1 r >= 0, = 0 only at r = 0 |
| 36 | salomon | 10 | E | 0 | 1 | no | 1 - cos(2 pi r) + 0.1 r >= 0, = 0 only at r = 0 |
| 37 | schaffer1 | 2 | E | 0 | 1 | yes | 0.5 + (sin^2 sqrt - 0.5)/(1+0.001 r^2)^2 >= 0, = 0 at 0 |
| 38 | schaffer2 | 2 | E | 0 | 1 | yes | as schaffer1 (variant) |
| 39 | schwefel | 10 | A | -4189.828873 | 1 | no | separable: n copies of the 1-D problem -x sin(sqrt/x/) on [-500,500], 1-D global at 420.9687 (fine 1-D search certifies); f* = -418.9829 n |
| 40 | shekel_45 | 4 | L | -10.15319968 | 1 | yes | -10.1532 at ~(4,4,4,4): literature (5 local minima); refined |
| 41 | shekel_47 | 4 | L | -10.40294057 | 1 | yes | -10.4029: literature (7 local minima); refined |
| 42 | shekel_410 | 4 | L | -10.53640982 | 1 | yes | -10.5364: literature (10 local minima); refined |
| 43 | shekel_foxholes | 5 | L | -10.40561724 | 1 | yes | Bersini et al. 5-D / 10-D (30 foxholes): literature; refined |
| 44 | shekel_foxholes | 10 | L | -10.20879279 | 1 | no | Bersini et al. 5-D / 10-D (30 foxholes): literature; refined |
| 45 | shubert | 2 | L | -186.7309088 | 12 | yes | -186.7309 at 18 points (e.g. (-7.0835,4.8580)): classical, extensively verified |
| 46 | sinusoidal | 10 | E | -3.5 | 1 | no | -(A prod sin + prod sin) >= -(A+1) = -3.5, attained iff all sines = 1: x_i = 120 (unique in [0,180]) |
| 47 | sixhumpcamel | 2 | A | -1.031628453 | 2 | yes | -1.031628 at ±(0.0898,-0.7127): classical; polynomial, globality by the standard analysis |
| 48 | sphere | 3 | E | 0 | 1 | yes | = 0 only at 0 |
| 49 | tenn_local_minima | 2 | E | 2.355817158e-31 | 1 | yes | sum of nonneg terms = 0 only at ones |
| 50 | tenn_local_minima | 4 | E | 1.177908579e-31 | 1 | yes | sum of nonneg terms = 0 only at ones |
| 51 | tenn_local_minima | 6 | E | 7.852723859e-32 | 1 | no | sum of nonneg terms = 0 only at ones |
| 52 | tenn_local_minima | 8 | E | 5.889542894e-32 | 1 | no | sum of nonneg terms = 0 only at ones |
| 53 | threehumpcamel | 2 | E | 0 | 1 | yes | 2x^2-1.05x^4+x^6/6+xy+y^2 >= 0, = 0 at 0 (standard) |
| 54 | transistor | 9 | E | 1.600493267e-24 | 1 | no | sum of squares of a system with a root; refined root gives 1.6e-24 (documented point gives 1.8e-7): f* = 0 up to machine precision; the root is the global minimizer |
| 55 | wood | 4 | E | 0 | 1 | yes | sum of squares = 0 at ones |

Totals: E 36 (enumerated 24), A 6 (enumerated 4), L 13 (enumerated 11). Instances not enumerated have 6 or more finite variables.
