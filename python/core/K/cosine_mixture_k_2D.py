"""cosine_mixture_k_2D  DFO Benchmark Suite v2.1.0 core, Python version (self-contained; standard library only).
COSINE_MIXTURE_K_2D  DFO Benchmark Suite v2.1.0 core: cosine_mixture, dimension 2, pattern K, strategy baseline.

 Copyright (C) 2026 J. F. A. Madeira

 Objective F(z) = f(T(z)), where f is the original function (embedded below, translated from the MATLAB original) and T acts on the finite
 coordinates only: uniform grid of 151 levels on the interval of the domain, piecewise-affine map with |T - Id| <=
 half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
 Global optimum of F on this domain: f* = -2.2 (E: exact by structure); the optimal solution is listed in
 reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
 Variables: x1:K x2:K


Python version of the MATLAB wrapper cosine_mixture_k_2D.m: same instance, same values (level indices 1..m as in MATLAB).
Usage:  from cosine_mixture_k_2D import cosine_mixture_k_2D;  PB = cosine_mixture_k_2D();  f = PB.func_idx(PB.x0_idx)
"""
# Copyright (C) 2026 J. F. A. Madeira. LGPL-3.0-or-later.
import math
from types import SimpleNamespace

_NAME = 'cosine_mixture_k_2D'
_N = 2
_TYPES = ['K', 'K']
_CARD = [151.0, 151.0]
_LB = [-0.73625732421875, -0.7006238886348626]
_UB = [1.0, 1.0]
_LBW = [-0.73625732421875, -0.7006238886348626]
_UBW = [1.0, 1.0]
_TB = [1.0, 1.0]
_TI = [1.0, 1.0]
_CW = [math.nan, math.nan]
_CX = [math.nan, math.nan]
_M = 151
_PERM = {0: [78, 17, 3, 19, 67, 88, 73, 105, 80, 143, 92, 62, 59, 90, 8, 82, 101, 123, 102, 44, 97, 50, 23, 36, 145, 52, 72, 116, 16, 40, 79, 136, 64, 148, 141, 65, 47, 151, 144, 18, 86, 122, 56, 94, 126, 104, 7, 89, 14, 137, 121, 111, 128, 98, 77, 119, 63, 46, 70, 134, 125, 31, 130, 113, 71, 5, 135, 147, 30, 53, 114, 11, 138, 87, 66, 106, 21, 91, 61, 150, 32, 48, 57, 15, 42, 115, 43, 54, 76, 132, 110, 55, 142, 34, 49, 1, 103, 6, 2, 120, 127, 28, 112, 117, 45, 9, 41, 107, 33, 108, 131, 133, 149, 140, 100, 38, 99, 75, 29, 139, 69, 27, 124, 12, 24, 39, 85, 4, 95, 35, 37, 68, 51, 118, 13, 146, 20, 81, 22, 58, 96, 83, 93, 84, 60, 26, 25, 109, 74, 10, 129], 1: [49, 128, 70, 86, 72, 133, 63, 109, 48, 89, 113, 93, 110, 114, 69, 9, 75, 12, 42, 29, 84, 50, 41, 132, 95, 129, 76, 126, 17, 74, 124, 61, 78, 62, 104, 2, 83, 67, 16, 143, 77, 28, 103, 81, 40, 3, 134, 98, 106, 21, 31, 94, 96, 65, 22, 30, 36, 117, 54, 7, 140, 102, 37, 118, 101, 127, 14, 27, 120, 97, 139, 121, 34, 91, 68, 11, 23, 43, 64, 13, 35, 105, 92, 18, 111, 123, 46, 57, 79, 71, 45, 59, 87, 125, 130, 149, 136, 145, 6, 147, 90, 33, 131, 1, 51, 15, 138, 141, 119, 146, 44, 58, 82, 66, 108, 47, 107, 85, 122, 32, 80, 135, 25, 144, 38, 142, 99, 20, 112, 53, 56, 88, 73, 8, 39, 55, 148, 116, 151, 100, 5, 137, 19, 24, 10, 4, 26, 52, 150, 115, 60]}

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

def _orig(x):
    s1 = 0.0
    for v in x: s1 += math.cos(5 * math.pi * v)
    s2 = 0.0
    for v in x: s2 += v ** 2
    return 0.1 * s1 - s2

def _T_map(g, a, b, zs, xs):
    # piecewise-affine T with T(zs) = xs; fixes the bounds when zs is interior
    if abs(g - zs) <= 1e-12 * max(1, abs(b - a)):
        return xs
    if abs(zs - a) <= 1e-14 * max(1, abs(a)):
        x = xs + (g - a) * (b - xs) / (b - a)
    elif abs(zs - b) <= 1e-14 * max(1, abs(b)):
        x = a + (g - a) * (xs - a) / (b - a)
    elif g <= zs:
        x = a + (g - a) * (xs - a) / (zs - a)
    else:
        x = xs + (g - zs) * (b - xs) / (b - zs)
    return min(max(x, min(a, xs)), max(b, xs))

def _map_work_to_orig(xw, lbw, ubw, lbo, ubo):
    if abs(ubw - lbw) <= 1e-14:
        t = 0
    else:
        t = (xw - lbw) / (ubw - lbw)
    t = max(0, min(1, t))
    return lbo + t * (ubo - lbo)

_E = {i: [g[p - 1] for p in perm] for i, perm in _PERM.items() for g in [_linspace(_LB[i], _UB[i], _M)]}

def _decode(z):
    z = _vector(z)
    xh = [0.0] * _N
    for i, t in enumerate(_TYPES):
        if t == 'C':
            v = z[i]; tl = 1e-12 * (_UBW[i] - _LBW[i])
            if v < _LBW[i] - tl or v > _UBW[i] + tl:
                _err('x%d = %g outside [%g, %g].' % (i + 1, v, _LBW[i], _UBW[i]))
            v = min(max(v, _LBW[i]), _UBW[i])
            if abs(v - _CW[i]) <= 1e-12 * (_UBW[i] - _LBW[i]):
                xh[i] = _CX[i]
            else:
                xh[i] = _map_work_to_orig(v, _LBW[i], _UBW[i], _LB[i], _UB[i])
        else:
            k = _check_idx(z[i], _M, i)
            if t == 'K':
                g = _E[i][k - 1]
            else:
                g = _LB[i] + (k - 1) * (_UB[i] - _LB[i]) / (_M - 1)
            xh[i] = _T_map(g, _LB[i], _UB[i], _TB[i], _TI[i])
    return xh

def _objective(z):
    return _orig(_decode(z))

def _problem_data():
    pd = []
    for i, t in enumerate(_TYPES):
        if t == 'C': pd.append([_LBW[i], _UBW[i]])
        elif t == 'D': pd.append([_linspace(_LB[i], _UB[i], _M)])
        else: pd.append(['k%d_%d' % (i + 1, j) for j in range(1, _M + 1)])
    return pd

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

def cosine_mixture_k_2D():
    PB = SimpleNamespace()
    PB.name = _NAME
    PB.nvar = 2
    PB.var_types = ['K', 'K']
    PB.cardinality = [151.0, 151.0]
    PB.f_star = -2.2
    PB.f_star_justification = 'E'
    PB.n_minimizers_documented = 1
    PB.description = 'DFO Benchmark Suite v2.1.0 core (Python): K adaptation of cosine_mixture (F = f o T).'
    PB.x0_idx = [10.0, 9.0]
    PB.x0_alt_idx = [[144.0, 13.0], [70.0, 3.0], [107.0, 12.0], [33.0, 14.0], [126.0, 4.0], [52.0, 8.0], [89.0, 18.0], [15.0, 5.0], [135.0, 9.0]]
    PB.suite_release = '2.1.0'
    PB.base_function = 'cosine_mixture'
    PB.m = 151
    PB.pattern = 'K'
    PB.strategy = 'baseline'
    PB.scale_factors_C = []
    PB.lb_orig = [-0.73625732421875, -0.7006238886348626]
    PB.ub_orig = [1.0, 1.0]
    PB.lb_C_work = [-0.73625732421875, -0.7006238886348626]
    PB.ub_C_work = [1.0, 1.0]
    PB.ProblemData = _problem_data()
    PB.func_idx = _objective
    PB.func_F = lambda xr: _objective(_cell2idx(xr))
    PB.idx2cell = _idx2cell
    PB.cell2idx = _cell2idx
    PB.decode = lambda xr: _decode(_cell2idx(xr))
    return PB
