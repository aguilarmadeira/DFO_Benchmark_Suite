"""Minimal use of the Python version of the DFO Benchmark Suite: evaluate one instance and run the quasi-random
search baseline of benchmark/solver_random_search.m with the protocol of benchmark/PROTOCOL.md (x0 first, budget
1000 n). The points are those of the Halton sequence from index 2^20 * seed + 1 (no random generator), as in MATLAB."""
import os, sys, math
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), 'core', 'CDK_baseline'))
from ackley_cdk_10D import ackley_cdk_10D

def primes(k):
    p, c = [], 2
    while len(p) < k:
        if all(c % q for q in p if q * q <= c): p.append(c)
        c += 1
    return p

def halton_point(idx, d):
    """point idx (1-based) of the d-dimensional Halton sequence, bases = first d primes (as benchmark/halton_point.m)"""
    x = []
    for b in primes(d):
        i, f, v = idx, 1 / b, 0.0
        while i > 0:
            v = v + f * (i % b); i = i // b; f = f / b
        x.append(v)
    return x

PB = ackley_cdk_10D(); n = PB.nvar; seed = 1; budget = 1000 * n
print(PB.name, n, PB.var_types, PB.cardinality, 'f* =', PB.f_star)
best = PB.func_idx(PB.x0_idx)                    # protocol: x0 is the first evaluation
for t in range(1, budget):
    u = halton_point(1048576 * seed + t, n)
    z = [PB.lb_C_work[i] + u[i] * (PB.ub_C_work[i] - PB.lb_C_work[i]) if c == 'C'
         else min(1 + math.floor(u[i] * PB.cardinality[i]), PB.cardinality[i]) for i, c in enumerate(PB.var_types)]
    best = min(best, PB.func_idx(z))
print('best value after', budget, 'evaluations:', best)
