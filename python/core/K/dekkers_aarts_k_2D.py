"""dekkers_aarts_k_2D  DFO Benchmark Suite v2.1.0 core, Python version (self-contained; standard library only).
DEKKERS_AARTS_K_2D  DFO Benchmark Suite v2.1.0 core: dekkers_aarts, dimension 2, pattern K, strategy baseline.

 Copyright (C) 2026 J. F. A. Madeira

 Objective F(z) = f(T(z)), where f is the original function (embedded below, translated from the MATLAB original) and T acts on the finite
 coordinates only: uniform grid of 151 levels on the interval of the domain, piecewise-affine map with |T - Id| <=
 half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
 Global optimum of F on this domain: f* = -24776.51834231769 (A: analytic reduction, 1-D root to machine precision); the optimal solution is listed in
 reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
 Variables: x1:K x2:K


Python version of the MATLAB wrapper dekkers_aarts_k_2D.m: same instance, same values (level indices 1..m as in MATLAB).
Usage:  from dekkers_aarts_k_2D import dekkers_aarts_k_2D;  PB = dekkers_aarts_k_2D();  f = PB.func_idx(PB.x0_idx)
"""
# Copyright (C) 2026 J. F. A. Madeira. LGPL-3.0-or-later.
import math
from types import SimpleNamespace

_NAME = 'dekkers_aarts_k_2D'
_N = 2
_TYPES = ['K', 'K']
_CARD = [151.0, 151.0]
_LB = [-20.0, -12.145811106030584]
_UB = [14.025146484375, 20.0]
_LBW = [-20.0, -12.145811106030584]
_UBW = [14.025146484375, 20.0]
_TB = [-0.03858072916666799, 14.856670223035108]
_TI = [0.0, 14.94511215764898]
_CW = [math.nan, math.nan]
_CX = [math.nan, math.nan]
_M = 151
_PERM = {0: [104, 16, 90, 131, 76, 136, 122, 139, 87, 132, 60, 74, 130, 66, 103, 91, 70, 46, 19, 1, 55, 44, 33, 134, 41, 138, 101, 22, 51, 92, 2, 135, 64, 96, 58, 98, 18, 54, 28, 49, 109, 40, 89, 118, 143, 59, 95, 25, 31, 105, 48, 78, 144, 26, 73, 69, 84, 11, 117, 30, 50, 57, 62, 150, 88, 34, 149, 147, 99, 15, 83, 45, 67, 110, 145, 151, 61, 100, 24, 53, 77, 115, 20, 63, 106, 75, 125, 21, 43, 137, 102, 97, 56, 7, 23, 72, 128, 82, 81, 3, 68, 13, 36, 86, 141, 6, 39, 107, 142, 5, 47, 148, 140, 120, 121, 71, 146, 4, 79, 80, 12, 114, 108, 112, 129, 29, 85, 27, 10, 94, 123, 32, 37, 8, 52, 116, 38, 14, 124, 93, 133, 35, 9, 119, 126, 17, 113, 65, 42, 127, 111], 1: [94, 104, 82, 92, 65, 55, 64, 115, 107, 21, 46, 67, 74, 130, 29, 105, 96, 137, 95, 2, 40, 116, 146, 34, 84, 15, 37, 71, 111, 39, 61, 3, 98, 85, 68, 147, 19, 48, 54, 60, 36, 120, 57, 5, 10, 62, 88, 81, 32, 38, 138, 122, 25, 28, 128, 52, 124, 106, 18, 145, 51, 93, 101, 113, 149, 35, 44, 144, 143, 102, 22, 119, 69, 7, 97, 50, 63, 108, 117, 56, 79, 151, 123, 66, 136, 135, 118, 24, 47, 103, 16, 30, 133, 131, 129, 76, 31, 150, 86, 148, 73, 17, 109, 91, 126, 110, 45, 6, 41, 12, 132, 72, 127, 87, 20, 80, 26, 13, 142, 90, 125, 4, 121, 1, 139, 141, 100, 43, 27, 112, 11, 53, 59, 33, 78, 9, 70, 49, 23, 134, 42, 89, 77, 58, 75, 140, 114, 8, 83, 99, 14]}

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
    return 10 ** 5 * x[0] ** 2 + x[1] ** 2 - (x[0] ** 2 + x[1] ** 2) ** 2 + 1 / 10 ** 5 * (x[0] ** 2 + x[1] ** 2) ** 4

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

def dekkers_aarts_k_2D():
    PB = SimpleNamespace()
    PB.name = _NAME
    PB.nvar = 2
    PB.var_types = ['K', 'K']
    PB.cardinality = [151.0, 151.0]
    PB.f_star = -24776.51834231769
    PB.f_star_justification = 'A'
    PB.n_minimizers_documented = 1
    PB.description = 'DFO Benchmark Suite v2.1.0 core (Python): K adaptation of dekkers_aarts (F = f o T).'
    PB.x0_idx = [48.0, 109.0]
    PB.x0_alt_idx = [[83.0, 113.0], [9.0, 103.0], [145.0, 112.0], [71.0, 114.0], [107.0, 104.0], [33.0, 108.0], [126.0, 118.0], [52.0, 105.0], [88.0, 109.0]]
    PB.suite_release = '2.1.0'
    PB.base_function = 'dekkers_aarts'
    PB.m = 151
    PB.pattern = 'K'
    PB.strategy = 'baseline'
    PB.scale_factors_C = []
    PB.lb_orig = [-20.0, -12.145811106030584]
    PB.ub_orig = [14.025146484375, 20.0]
    PB.lb_C_work = [-20.0, -12.145811106030584]
    PB.ub_C_work = [14.025146484375, 20.0]
    PB.ProblemData = _problem_data()
    PB.func_idx = _objective
    PB.func_F = lambda xr: _objective(_cell2idx(xr))
    PB.idx2cell = _idx2cell
    PB.cell2idx = _cell2idx
    PB.decode = lambda xr: _decode(_cell2idx(xr))
    return PB
