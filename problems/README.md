# Test Problems

This directory contains MATLAB implementations of 63 bound-constrained test problems for derivative-free global optimization.

## Problem List

| # | Problem | n | Bounds | Local | Global |
|---|---------|---|--------|-------|--------|
| 1 | Ackley | 10 | [-30, 30] | many | 1 |
| 2 | Aluffi-Pentini | 2 | [-10, 10] | 2 | 1 |
| 3 | Becker-Lago | 2 | [-10, 10] | 4 | 4 |
| 4 | Bohachevsky | 2 | [-50, 50] | many | 1 |
| 5 | Branin-Hoo | 2 | [-5,10]×[0,15] | 3 | 3 |
| 6 | Cauchy | 4 | [3, 17] | many | — |
| 7 | Cauchy | 10 | [2, 26] | many | — |
| 8 | Cosine Mixture | 2 | [-1, 1] | many | — |
| 9 | Cosine Mixture | 4 | [-1, 1] | many | — |
| 10 | Dekkers-Aarts | 2 | [-20, 20] | 3 | 2 |
| 11 | Epistatic Michalewicz | 5 | [0, π] | many | 1 |
| 12 | Epistatic Michalewicz | 10 | [0, π] | many | 1 |
| 13 | Exponential | 2 | [-1, 1] | — | 1 |
| 14 | Exponential | 4 | [-1, 1] | — | 1 |
| 15 | Fifteen Local Minima | 2 | [-10, 10] | 15² | 1 |
| 16 | Fifteen Local Minima | 4 | [-10, 10] | 15⁴ | 1 |
| 17 | Fifteen Local Minima | 6 | [-10, 10] | 15⁶ | 1 |
| 18 | Fifteen Local Minima | 8 | [-10, 10] | 15⁸ | 1 |
| 19 | Fifteen Local Minima | 10 | [-10, 10] | 15¹⁰ | 1 |
| 20 | Fletcher-Powell | 3 | [-10, 10] | many | 1 |
| 21 | Goldstein-Price | 2 | [-2, 2] | 4 | 1 |
| 22 | Griewank | 10 | [-400, 400] | many | 1 |
| 23 | Griewank | 5 | [-600, 600] | many | 1 |
| 24 | Gulf | 3 | [0.1,100]×[0,25.6]×[0,5] | — | 1 |
| 25 | Hartman 4 | 3 | [0, 1] | 4 | 1 |
| 26 | Hartman 4 | 6 | [0, 1] | 4 | 1 |
| 27 | Hosaki | 2 | [0,5]×[0,6] | 2 | 1 |
| 28 | Kowalik | 4 | [0, 0.42] | — | 1 |
| 29 | Langerman | 10 | [0, 10] | many | 1 |
| 30 | McCormick | 2 | [-1.5,4]×[-3,3] | 2 | 1 |
| 31 | Meyer-Roth | 3 | [-10, 10] | — | 1 |
| 32 | Miele-Cantrell | 4 | [-10, 10] | — | 1 |
| 33 | Multi-Gaussian | 2 | [-2, 2] | 5 | 1 |
| 34 | Neumaier 2 | 4 | [0, 4] | — | 1 |
| 35 | Neumaier 3 | 10 | [-100, 100] | — | 1 |
| 36 | Paviani | 10 | [2.001, 9.999] | — | 1 |
| 37 | Periodic | 2 | [-10, 10] | 50 | 1 |
| 38 | Poissonian | 2 | [1,21]×[1,8] | — | — |
| 39 | Powell | 4 | [-10, 10] | 1 | 1 |
| 40 | Rastrigin | 10 | [-5.12, 5.12] | many | 1 |
| 41 | Rosenbrock | 2 | [-5.12, 5.12] | 1 | 1 |
| 42 | Salomon | 5 | [-100, 100] | — | 1 |
| 43 | Salomon | 10 | [-100, 100] | — | 1 |
| 44 | Schaffer 1 | 2 | [-100, 100] | many | 1 |
| 45 | Schaffer 2 | 2 | [-100, 100] | many | 1 |
| 46 | Schwefel | 10 | [-500, 500] | many | 1 |
| 47 | Shekel 4-5 | 4 | [0, 10] | 5 | 1 |
| 48 | Shekel 4-7 | 4 | [0, 10] | 7 | 1 |
| 49 | Shekel 4-10 | 4 | [0, 10] | 10 | 1 |
| 50 | Shekel Foxholes | 5 | [0, 10] | many | 1 |
| 51 | Shekel Foxholes | 10 | [0, 10] | many | 1 |
| 52 | Shubert | 2 | [-10, 10] | 760 | 18 |
| 53 | Sinusoidal | 10 | [0, 180] | — | 1 |
| 54 | Six-Hump Camel | 2 | [-3,3]×[-2,2] | 6 | 2 |
| 55 | Sphere | 3 | [-5.12, 5.12] | 1 | 1 |
| 56 | Storn-Tchebychev | 9 | [-128, 128] | — | 1 |
| 57 | Ten Local Minima | 2 | [-10, 10] | 10² | 1 |
| 58 | Ten Local Minima | 4 | [-10, 10] | 10⁴ | 1 |
| 59 | Ten Local Minima | 6 | [-10, 10] | 10⁶ | 1 |
| 60 | Ten Local Minima | 8 | [-10, 10] | 10⁸ | 1 |
| 61 | Three-Hump Camel | 2 | [-5, 5] | 3 | 1 |
| 62 | Transistor | 9 | [-10, 10] | — | 1 |
| 63 | Wood | 4 | [-10, 10] | — | 1 |

**Note:** "—" indicates unknown or very large number of local minima.

## Usage

Each problem is implemented as a MATLAB function:

```matlab
% Example: Ackley function
f = ackley(x);  % x is a column vector of dimension n
```

## Function Signature

All problem functions follow the same interface:

```matlab
function f = problem_name(x)
% PROBLEM_NAME  Evaluates the problem_name test function
%   f = problem_name(x) returns the objective value at x
%
% Input:
%   x - column vector of dimension n
%
% Output:
%   f - scalar objective value
```

## References

1. Ali, M.M., Khompatraporn, C., Zabinsky, Z.B. (2005). A numerical evaluation of several stochastic algorithms on selected continuous global optimization test problems. *J. Global Optim.* 31, 635–672.

2. Brachetti, P., Ciccoli, M.F., Di Pillo, G., Lucidi, S. (1997). A new version of the Price's algorithm for global optimization. *J. Global Optim.* 10, 165–184.

3. Custódio, A.L., Madeira, J.F.A. (2015). GLODS: Global and local optimization using direct search. *J. Global Optim.* 62, 1–28.
