"""interBK_nC2_K1x16_s1_b0p018_l1  DFO Benchmark Suite v2.1.0 core, Python version (self-contained; standard library only).
INTERBK_NC2_K1X16_S1_B0P018_L1  DFO Benchmark Suite v2.1.0 core, Family BK: f = ||x_C - x_C* - lambda s(kappa)||^2 + b (1 - h(kappa)/nK).

 Copyright (C) 2026 J. F. A. Madeira

 nC = 2 continuous variables in [0,1]; nK = 1 categorical variables with m = 16 categories each (M_K = 16 combinations)
 lambda = 1, sigma = 0.15, b = 0.018, q = lambda^2 sigma^2 / b = 1.25, j_min = Inf (minimal number of simultaneous categorical
 changes from an antipodal conditional minimiser that improves f without continuous re-optimisation; Inf: none; 0: lambda = 0)
 h(kappa) = number of categorical coordinates equal to the optimal combination; s(kappa) = s if h = 0, else 0.
 Global optimum f* = 0 (E: by construction); the optimal solution is listed in reference_optima.csv of this
 folder, not here. Starting point of the protocol: PB.x0_idx (Halton).


Python version of the MATLAB wrapper interBK_nC2_K1x16_s1_b0p018_l1.m: same instance, same values (level indices 1..m as in MATLAB).
Usage:  from interBK_nC2_K1x16_s1_b0p018_l1 import interBK_nC2_K1x16_s1_b0p018_l1;  PB = interBK_nC2_K1x16_s1_b0p018_l1();  f = PB.func_idx(PB.x0_idx)
"""
# Copyright (C) 2026 J. F. A. Madeira. LGPL-3.0-or-later.
import math
from types import SimpleNamespace

_NAME = 'interBK_nC2_K1x16_s1_b0p018_l1'
_N = 3
_TYPES = ['C', 'C', 'K']
_CARD = [math.inf, math.inf, 16.0]
_NC = 2
_NK = 1
_MK = 16
_XS = [0.6693599700927735, 0.31318170784715516]
_S = [-0.07721356559429048, 0.12860040936255296]
_L2S = [[15, 1, 3, 4, 13, 9, 12, 14, 6, 0, 11, 8, 10, 7, 5, 2]]
_KS = [14]
_LAMBDA = 1.0
_B = 0.018

def _linspace(a, b, m):
    # MATLAB linspace: a + (k*(b-a))/(m-1), first and last values exactly a and b
    n1 = m - 1
    y = [a + (k * (b - a)) / n1 for k in range(m)]
    y[0] = a; y[-1] = b
    return y

def _err(msg):
    raise ValueError(_NAME + ': ' + msg)

def _vector(z):
    try:
        z = [float(v) for v in z]
    except Exception:
        _err('z must be a finite numeric vector with %d entries.' % _N)
    if len(z) != _N or not all(math.isfinite(v) for v in z):
        _err('z must be a finite numeric vector with %d entries.' % _N)
    return z

def _check_idx(v, m, i):
    if not math.isfinite(v) or v != math.floor(v) or v < 1 or v > m:
        _err('z%d = %g is not an index in 1..%d.' % (i + 1, v, m))
    return int(v)

def _objective(z):
    z = _vector(z)
    u = z[:_NC]
    if any(v < 0 or v > 1 for v in u):
        _err('continuous coordinates must lie in [0,1].')
    kk = z[_NC:]
    if any(k != math.floor(k) or k < 1 or k > _MK for k in kk):
        _err('label indices must be integers in 1..%d.' % _MK)
    h = 0
    for i in range(_NK):
        h += 1 if _L2S[i][int(kk[i]) - 1] == _KS[i] else 0
    shift = [_LAMBDA * s for s in _S] if h == 0 else [0.0] * _NC
    f = 0.0
    for j in range(_NC):
        d = u[j] - _XS[j] - shift[j]
        f += d * d
    return f + _B * (1 - h / _NK)

def _problem_data():
    return [[0.0, 1.0] for _ in range(_NC)] + [['k%d_%d' % (i + 1, j) for j in range(1, _MK + 1)] for i in range(_NC, _N)]

def _idx2cell(z):
    z = _vector(z); pd = _problem_data(); xr = []
    for i, t in enumerate(_TYPES):
        if t == 'C': xr.append(z[i])
        elif t == 'D': xr.append(pd[i][0][_check_idx(z[i], _CARD[i], i) - 1])
        else: xr.append(pd[i][_check_idx(z[i], _CARD[i], i) - 1])
    return xr

def _cell2idx(xr):
    if not isinstance(xr, (list, tuple)) or len(xr) != _N:
        _err('xr must be a list with %d entries.' % _N)
    pd = _problem_data(); z = [0.0] * _N
    for i, t in enumerate(_TYPES):
        v = xr[i]
        if t == 'C':
            try:
                v = float(v)
            except (TypeError, ValueError):
                _err('x%d must be a finite scalar.' % (i + 1))
            if isinstance(xr[i], (str, bytes, bool)) or not math.isfinite(v):
                _err('x%d must be a finite scalar.' % (i + 1))
            z[i] = v
        elif t == 'D':
            try:
                v = float(v)
            except (TypeError, ValueError):
                _err('x%d must be a grid value.' % (i + 1))
            if isinstance(xr[i], (str, bytes, bool)) or not math.isfinite(v):
                _err('x%d must be a grid value.' % (i + 1))
            grid = pd[i][0]
            d, k = min((abs(gv - v), j) for j, gv in enumerate(grid))
            if d > 1e-10 * max(1, max(abs(gv) for gv in grid)):
                _err('x%d = %g is not a grid value.' % (i + 1, v))
            z[i] = k + 1
        else:
            if not isinstance(v, str) or v not in pd[i]:
                _err('x%d unknown label %s.' % (i + 1, v))
            z[i] = pd[i].index(v) + 1
    return z

def interBK_nC2_K1x16_s1_b0p018_l1():
    PB = SimpleNamespace()
    PB.name = _NAME
    PB.nvar = 3
    PB.var_types = ['C', 'C', 'K']
    PB.cardinality = [math.inf, math.inf, 16.0]
    PB.f_star = 0.0
    PB.f_star_justification = 'E'
    PB.n_minimizers_documented = 1
    PB.description = 'DFO Benchmark Suite v2.1.0 core, Family BK (controlled categorical dimension): 1 categorical variables x 16 categories.'
    PB.x0_idx = [0.85888671875, 0.5098308184727939, 5.0]
    PB.x0_alt_idx = [[0.343994140625, 0.5386374028349337, 16.0], [0.601806640625, 0.4686785550983082, 8.0], [0.097900390625, 0.4816338972717573, 1.0], [0.918212890625, 0.5433622923334858, 11.0], [0.408447265625, 0.4734034445968602, 4.0], [0.666259765625, 0.5022100289590001, 14.0], [0.162353515625, 0.5200426764212772, 7.0], [0.795166015625, 0.4500838286846517, 3.0], [0.282470703125, 0.5118122237463801, 10.0]]
    PB.suite_release = '2.1.0'
    PB.configuration = 'interBK_nC2_K1x16'
    PB.lambda_ = 1.0
    PB.sigma = 0.15
    PB.b = 0.018
    PB.q = 1.25
    PB.j_min = math.inf
    PB.obstacle_ratio = 1.25
    PB.seed = 1.0
    PB.relabel = 1.0
    PB.base_function = 'familyBK'
    PB.family = 'BK'
    PB.pattern = 'BK'
    PB.nK = 1
    PB.m = 16
    PB.lb_orig = [0.0, 0.0, 0.0]
    PB.ub_orig = [1.0, 1.0, 1.0]
    PB.lb_C_work = [0.0, 0.0, 0.0]
    PB.ub_C_work = [1.0, 1.0, 1.0]
    PB.ProblemData = _problem_data()
    PB.func_idx = _objective
    PB.func_F = lambda xr: _objective(_cell2idx(xr))
    PB.idx2cell = _idx2cell
    PB.cell2idx = _cell2idx
    return PB
