"""interBK_nC2_K6x2_s3_b0p054_l0  DFO Benchmark Suite v2.1.0 core, Python version (self-contained; standard library only).
INTERBK_NC2_K6X2_S3_B0P054_L0  DFO Benchmark Suite v2.1.0 core, Family BK: f = ||x_C - x_C* - lambda s(kappa)||^2 + b (1 - h(kappa)/nK).

 Copyright (C) 2026 J. F. A. Madeira

 nC = 2 continuous variables in [0,1]; nK = 6 categorical variables with m = 2 categories each (M_K = 64 combinations)
 lambda = 0, sigma = 0.15, b = 0.05399999999999999, q = lambda^2 sigma^2 / b = 0, j_min = 0 (minimal number of simultaneous categorical
 changes from an antipodal conditional minimiser that improves f without continuous re-optimisation; Inf: none; 0: lambda = 0)
 h(kappa) = number of categorical coordinates equal to the optimal combination; s(kappa) = s if h = 0, else 0.
 Global optimum f* = 0 (E: by construction); the optimal solution is listed in reference_optima.csv of this
 folder, not here. Starting point of the protocol: PB.x0_idx (Halton).


Python version of the MATLAB wrapper interBK_nC2_K6x2_s3_b0p054_l0.m: same instance, same values (level indices 1..m as in MATLAB).
Usage:  from interBK_nC2_K6x2_s3_b0p054_l0 import interBK_nC2_K6x2_s3_b0p054_l0;  PB = interBK_nC2_K6x2_s3_b0p054_l0();  f = PB.func_idx(PB.x0_idx)
"""
# Copyright (C) 2026 J. F. A. Madeira. LGPL-3.0-or-later.
import math
from types import SimpleNamespace

_NAME = 'interBK_nC2_K6x2_s3_b0p054_l0'
_N = 8
_TYPES = ['C', 'C', 'K', 'K', 'K', 'K', 'K', 'K']
_CARD = [math.inf, math.inf, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0]
_NC = 2
_NK = 6
_MK = 2
_XS = [0.3561700820922852, 0.434638275932794]
_S = [0.12024173457180934, 0.08967678220789681]
_L2S = [[1, 0], [1, 0], [0, 1], [0, 1], [0, 1], [0, 1]]
_KS = [1, 1, 1, 0, 0, 0]
_LAMBDA = 0.0
_B = 0.05399999999999999

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

def interBK_nC2_K6x2_s3_b0p054_l0():
    PB = SimpleNamespace()
    PB.name = _NAME
    PB.nvar = 8
    PB.var_types = ['C', 'C', 'K', 'K', 'K', 'K', 'K', 'K']
    PB.cardinality = [math.inf, math.inf, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0]
    PB.f_star = 0.0
    PB.f_star_justification = 'E'
    PB.n_minimizers_documented = 1
    PB.description = 'DFO Benchmark Suite v2.1.0 core, Family BK (controlled categorical dimension): 6 categorical variables x 2 categories.'
    PB.x0_idx = [0.608642578125, 0.05410760554793477, 2.0, 1.0, 2.0, 1.0, 1.0, 2.0]
    PB.x0_alt_idx = [[0.094482421875, 0.08291418991007468, 1.0, 1.0, 2.0, 1.0, 1.0, 2.0], [0.914794921875, 0.0019814052735863436, 2.0, 1.0, 1.0, 1.0, 2.0, 1.0], [0.410888671875, 0.06370980033531473, 2.0, 1.0, 2.0, 2.0, 2.0, 2.0], [0.668701171875, 0.09251638469745464, 1.0, 1.0, 2.0, 2.0, 1.0, 2.0], [0.158935546875, 0.022557536960829135, 2.0, 1.0, 1.0, 2.0, 2.0, 1.0], [0.791748046875, 0.04039018442310623, 1.0, 1.0, 1.0, 1.0, 2.0, 2.0], [0.287841796875, 0.10211857948483462, 2.0, 1.0, 2.0, 1.0, 1.0, 2.0], [0.545654296875, 0.03215973174820912, 1.0, 1.0, 2.0, 1.0, 1.0, 1.0], [0.032958984375, 0.06096631611034903, 1.0, 1.0, 1.0, 2.0, 2.0, 2.0]]
    PB.suite_release = '2.1.0'
    PB.configuration = 'interBK_nC2_K6x2'
    PB.lambda_ = 0.0
    PB.sigma = 0.15
    PB.b = 0.05399999999999999
    PB.q = 0.0
    PB.j_min = 0.0
    PB.obstacle_ratio = 0.0
    PB.seed = 3.0
    PB.relabel = 1.0
    PB.base_function = 'familyBK'
    PB.family = 'BK'
    PB.pattern = 'BK'
    PB.nK = 6
    PB.m = 2
    PB.lb_orig = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    PB.ub_orig = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
    PB.lb_C_work = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    PB.ub_C_work = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
    PB.ProblemData = _problem_data()
    PB.func_idx = _objective
    PB.func_F = lambda xr: _objective(_cell2idx(xr))
    PB.idx2cell = _idx2cell
    PB.cell2idx = _cell2idx
    return PB
