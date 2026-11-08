"""interBK_nC2_K1x16_s3_b0p018_l0  DFO Benchmark Suite v2.1.0 core, Python version (self-contained; standard library only).
INTERBK_NC2_K1X16_S3_B0P018_L0  DFO Benchmark Suite v2.1.0 core, Family BK: f = ||x_C - x_C* - lambda s(kappa)||^2 + b (1 - h(kappa)/nK).

 Copyright (C) 2026 J. F. A. Madeira

 nC = 2 continuous variables in [0,1]; nK = 1 categorical variables with m = 16 categories each (M_K = 16 combinations)
 lambda = 0, sigma = 0.15, b = 0.018, q = lambda^2 sigma^2 / b = 0, j_min = 0 (minimal number of simultaneous categorical
 changes from an antipodal conditional minimiser that improves f without continuous re-optimisation; Inf: none; 0: lambda = 0)
 h(kappa) = number of categorical coordinates equal to the optimal combination; s(kappa) = s if h = 0, else 0.
 Global optimum f* = 0 (E: by construction); the optimal solution is listed in reference_optima.csv of this
 folder, not here. Starting point of the protocol: PB.x0_idx (Halton).


Python version of the MATLAB wrapper interBK_nC2_K1x16_s3_b0p018_l0.m: same instance, same values (level indices 1..m as in MATLAB).
Usage:  from interBK_nC2_K1x16_s3_b0p018_l0 import interBK_nC2_K1x16_s3_b0p018_l0;  PB = interBK_nC2_K1x16_s3_b0p018_l0();  f = PB.func_idx(PB.x0_idx)
"""
# Copyright (C) 2026 J. F. A. Madeira. LGPL-3.0-or-later.
import math
from types import SimpleNamespace

_NAME = 'interBK_nC2_K1x16_s3_b0p018_l0'
_N = 3
_TYPES = ['C', 'C', 'K']
_CARD = [math.inf, math.inf, 16.0]
_NC = 2
_NK = 1
_MK = 16
_XS = [0.6114986419677735, 0.6926970256340776]
_S = [0.09832831896164666, -0.11327639511467813]
_L2S = [[4, 15, 10, 5, 9, 11, 3, 14, 8, 7, 0, 12, 13, 1, 2, 6]]
_KS = [5]
_LAMBDA = 0.0
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

def interBK_nC2_K1x16_s3_b0p018_l0():
    PB = SimpleNamespace()
    PB.name = _NAME
    PB.nvar = 3
    PB.var_types = ['C', 'C', 'K']
    PB.cardinality = [math.inf, math.inf, 16.0]
    PB.f_star = 0.0
    PB.f_star_justification = 'E'
    PB.n_minimizers_documented = 1
    PB.description = 'DFO Benchmark Suite v2.1.0 core, Family BK (controlled categorical dimension): 1 categorical variables x 16 categories.'
    PB.x0_idx = [0.73388671875, 0.28760859625057156, 11.0]
    PB.x0_alt_idx = [[0.218994140625, 0.31641518061271146, 4.0], [0.851806640625, 0.24645633287608593, 14.0], [0.347900390625, 0.2594116750495351, 7.0], [0.605712890625, 0.3211400701112635, 3.0], [0.095947265625, 0.25118122237463797, 10.0], [0.916259765625, 0.2799878067367779, 6.0], [0.412353515625, 0.297820454199055, 13.0], [0.670166015625, 0.2278616064624295, 9.0], [0.157470703125, 0.28959000152415787, 2.0]]
    PB.suite_release = '2.1.0'
    PB.configuration = 'interBK_nC2_K1x16'
    PB.lambda_ = 0.0
    PB.sigma = 0.15
    PB.b = 0.018
    PB.q = 0.0
    PB.j_min = 0.0
    PB.obstacle_ratio = 0.0
    PB.seed = 3.0
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
