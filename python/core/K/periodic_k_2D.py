"""periodic_k_2D  DFO Benchmark Suite v2.1.0 core, Python version (self-contained; standard library only).
PERIODIC_K_2D  DFO Benchmark Suite v2.1.0 core: periodic, dimension 2, pattern K, strategy baseline.

 Copyright (C) 2026 J. F. A. Madeira

 Objective F(z) = f(T(z)), where f is the original function (embedded below, translated from the MATLAB original) and T acts on the finite
 coordinates only: uniform grid of 151 levels on the interval of the domain, piecewise-affine map with |T - Id| <=
 half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
 Global optimum of F on this domain: f* = 0.9 (E: exact by structure); the optimal solution is listed in
 reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
 Variables: x1:K x2:K


Python version of the MATLAB wrapper periodic_k_2D.m: same instance, same values (level indices 1..m as in MATLAB).
Usage:  from periodic_k_2D import periodic_k_2D;  PB = periodic_k_2D();  f = PB.func_idx(PB.x0_idx)
"""
# Copyright (C) 2026 J. F. A. Madeira. LGPL-3.0-or-later.
import math
from types import SimpleNamespace

_NAME = 'periodic_k_2D'
_N = 2
_TYPES = ['K', 'K']
_CARD = [151.0, 151.0]
_LB = [-10.0, -10.0]
_UB = [6.3563232421875, 6.055621602398009]
_LBW = [-10.0, -10.0]
_UBW = [6.3563232421875, 6.055621602398009]
_TB = [0.03187825520833343, -0.04551460651323325]
_TI = [0.0, 0.0]
_CW = [math.nan, math.nan]
_CX = [math.nan, math.nan]
_M = 151
_PERM = {0: [29, 55, 117, 132, 31, 140, 43, 25, 131, 78, 18, 8, 11, 27, 115, 52, 73, 108, 75, 86, 76, 145, 126, 97, 4, 41, 46, 17, 21, 96, 1, 67, 150, 65, 49, 15, 141, 34, 51, 72, 7, 13, 151, 23, 60, 125, 109, 101, 102, 22, 136, 47, 139, 69, 74, 116, 122, 148, 39, 147, 144, 91, 114, 82, 10, 28, 40, 118, 83, 128, 38, 146, 33, 134, 9, 120, 54, 95, 35, 77, 30, 53, 66, 44, 84, 12, 42, 58, 2, 32, 36, 56, 5, 70, 119, 100, 45, 26, 68, 99, 103, 110, 92, 59, 123, 24, 71, 61, 138, 130, 149, 129, 57, 87, 90, 111, 127, 104, 93, 79, 88, 14, 112, 48, 63, 20, 85, 81, 142, 113, 107, 19, 50, 143, 105, 37, 137, 3, 106, 62, 64, 80, 133, 98, 135, 89, 121, 94, 6, 124, 16], 1: [69, 22, 44, 19, 112, 71, 49, 90, 109, 120, 28, 122, 126, 140, 143, 53, 6, 33, 151, 136, 91, 46, 20, 30, 145, 92, 13, 18, 21, 76, 147, 131, 74, 75, 108, 16, 24, 94, 55, 85, 8, 14, 17, 101, 115, 63, 150, 107, 57, 142, 102, 149, 25, 66, 51, 78, 105, 70, 41, 59, 27, 138, 43, 54, 111, 106, 42, 88, 83, 81, 40, 116, 35, 52, 11, 80, 82, 7, 144, 79, 137, 56, 146, 47, 141, 15, 45, 61, 5, 36, 132, 124, 9, 99, 2, 133, 50, 31, 73, 3, 119, 62, 128, 110, 86, 29, 125, 67, 10, 118, 12, 87, 64, 95, 121, 139, 93, 113, 103, 89, 98, 26, 135, 60, 77, 34, 100, 129, 97, 1, 127, 39, 72, 23, 130, 65, 32, 38, 4, 104, 114, 134, 58, 148, 84, 68, 123, 37, 48, 96, 117]}

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
    return 1 + math.sin(x[0]) ** 2 + math.sin(x[1]) ** 2 - 0.1 * math.exp(-x[0] ** 2 - x[1] ** 2)

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

def periodic_k_2D():
    PB = SimpleNamespace()
    PB.name = _NAME
    PB.nvar = 2
    PB.var_types = ['K', 'K']
    PB.cardinality = [151.0, 151.0]
    PB.f_star = 0.9
    PB.f_star_justification = 'E'
    PB.n_minimizers_documented = 1
    PB.description = 'DFO Benchmark Suite v2.1.0 core (Python): K adaptation of periodic (F = f o T).'
    PB.x0_idx = [104.0, 111.0]
    PB.x0_alt_idx = [[30.0, 115.0], [124.0, 104.0], [50.0, 107.0], [86.0, 116.0], [12.0, 106.0], [133.0, 110.0], [59.0, 112.0], [95.0, 102.0], [21.0, 111.0]]
    PB.suite_release = '2.1.0'
    PB.base_function = 'periodic'
    PB.m = 151
    PB.pattern = 'K'
    PB.strategy = 'baseline'
    PB.scale_factors_C = []
    PB.lb_orig = [-10.0, -10.0]
    PB.ub_orig = [6.3563232421875, 6.055621602398009]
    PB.lb_C_work = [-10.0, -10.0]
    PB.ub_C_work = [6.3563232421875, 6.055621602398009]
    PB.ProblemData = _problem_data()
    PB.func_idx = _objective
    PB.func_F = lambda xr: _objective(_cell2idx(xr))
    PB.idx2cell = _idx2cell
    PB.cell2idx = _cell2idx
    PB.decode = lambda xr: _decode(_cell2idx(xr))
    return PB
