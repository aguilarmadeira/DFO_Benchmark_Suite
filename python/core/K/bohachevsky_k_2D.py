"""bohachevsky_k_2D  DFO Benchmark Suite v2.1.0 core, Python version (self-contained; standard library only).
BOHACHEVSKY_K_2D  DFO Benchmark Suite v2.1.0 core: bohachevsky, dimension 2, pattern K, strategy baseline.

 Copyright (C) 2026 J. F. A. Madeira

 Objective F(z) = f(T(z)), where f is the original function (embedded below, translated from the MATLAB original) and T acts on the finite
 coordinates only: uniform grid of 151 levels on the interval of the domain, piecewise-affine map with |T - Id| <=
 half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
 Global optimum of F on this domain: f* = 0 (E: exact by structure); the optimal solution is listed in
 reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
 Variables: x1:K x2:K


Python version of the MATLAB wrapper bohachevsky_k_2D.m: same instance, same values (level indices 1..m as in MATLAB).
Usage:  from bohachevsky_k_2D import bohachevsky_k_2D;  PB = bohachevsky_k_2D();  f = PB.func_idx(PB.x0_idx)
"""
# Copyright (C) 2026 J. F. A. Madeira. LGPL-3.0-or-later.
import math
from types import SimpleNamespace

_NAME = 'bohachevsky_k_2D'
_N = 2
_TYPES = ['K', 'K']
_CARD = [151.0, 151.0]
_LB = [-35.3909912109375, -31.920083320632017]
_UB = [50.0, 50.0]
_LBW = [-35.3909912109375, -31.920083320632017]
_UBW = [50.0, 50.0]
_TB = [-0.09604817708333258, -0.24431776998763866]
_TI = [0.0, 0.0]
_CW = [math.nan, math.nan]
_CX = [math.nan, math.nan]
_M = 151
_PERM = {0: [105, 55, 25, 50, 93, 64, 140, 143, 23, 54, 76, 127, 11, 148, 72, 82, 65, 150, 151, 149, 68, 99, 10, 75, 36, 77, 79, 120, 119, 43, 102, 95, 52, 12, 108, 16, 21, 117, 109, 33, 69, 146, 97, 5, 9, 118, 94, 103, 29, 61, 71, 122, 101, 67, 73, 42, 83, 134, 132, 17, 35, 41, 74, 89, 84, 13, 110, 116, 37, 124, 136, 114, 34, 115, 123, 145, 111, 113, 57, 4, 112, 90, 47, 2, 142, 3, 104, 32, 125, 39, 129, 20, 130, 87, 7, 98, 1, 51, 40, 60, 18, 131, 81, 22, 62, 14, 44, 106, 30, 92, 19, 24, 38, 15, 45, 88, 28, 100, 49, 141, 63, 139, 66, 58, 53, 59, 85, 86, 8, 46, 126, 6, 78, 135, 70, 137, 121, 80, 96, 133, 26, 27, 107, 31, 91, 56, 128, 48, 144, 147, 138], 1: [139, 12, 57, 149, 126, 51, 104, 84, 143, 106, 8, 117, 142, 92, 135, 82, 75, 3, 124, 112, 32, 9, 134, 78, 83, 115, 63, 26, 145, 89, 44, 123, 127, 148, 66, 137, 125, 29, 49, 105, 38, 52, 58, 80, 141, 70, 101, 87, 85, 61, 46, 1, 19, 21, 88, 62, 48, 2, 7, 18, 121, 43, 74, 54, 122, 131, 77, 140, 39, 150, 65, 17, 118, 47, 69, 10, 144, 40, 119, 103, 132, 90, 93, 138, 30, 6, 111, 35, 56, 42, 95, 24, 107, 91, 11, 110, 5, 55, 113, 128, 146, 130, 86, 28, 68, 20, 50, 4, 36, 98, 147, 31, 45, 22, 53, 97, 120, 109, 59, 37, 73, 108, 79, 71, 67, 72, 99, 102, 25, 64, 15, 27, 100, 94, 96, 13, 16, 116, 133, 34, 76, 81, 33, 114, 151, 23, 129, 136, 14, 60, 41]}

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
    return x[0] ** 2 + 2 * x[1] ** 2 - 0.3 * math.cos(3 * math.pi * x[0]) * math.cos(4 * math.pi * x[1]) + 0.3

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

def bohachevsky_k_2D():
    PB = SimpleNamespace()
    PB.name = _NAME
    PB.nvar = 2
    PB.var_types = ['K', 'K']
    PB.cardinality = [151.0, 151.0]
    PB.f_star = 0.0
    PB.f_star_justification = 'E'
    PB.n_minimizers_documented = 1
    PB.description = 'DFO Benchmark Suite v2.1.0 core (Python): K adaptation of bohachevsky (F = f o T).'
    PB.x0_idx = [34.0, 120.0]
    PB.x0_alt_idx = [[125.0, 124.0], [51.0, 133.0], [88.0, 123.0], [15.0, 126.0], [135.0, 130.0], [61.0, 119.0], [98.0, 129.0], [24.0, 131.0], [116.0, 120.0]]
    PB.suite_release = '2.1.0'
    PB.base_function = 'bohachevsky'
    PB.m = 151
    PB.pattern = 'K'
    PB.strategy = 'baseline'
    PB.scale_factors_C = []
    PB.lb_orig = [-35.3909912109375, -31.920083320632017]
    PB.ub_orig = [50.0, 50.0]
    PB.lb_C_work = [-35.3909912109375, -31.920083320632017]
    PB.ub_C_work = [50.0, 50.0]
    PB.ProblemData = _problem_data()
    PB.func_idx = _objective
    PB.func_F = lambda xr: _objective(_cell2idx(xr))
    PB.idx2cell = _idx2cell
    PB.cell2idx = _cell2idx
    PB.decode = lambda xr: _decode(_cell2idx(xr))
    return PB
