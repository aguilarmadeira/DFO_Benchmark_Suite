"""branin_hoo_k_2D  DFO Benchmark Suite v2.1.0 core, Python version (self-contained; standard library only).
BRANIN_HOO_K_2D  DFO Benchmark Suite v2.1.0 core: branin_hoo, dimension 2, pattern K, strategy baseline.

 Copyright (C) 2026 J. F. A. Madeira

 Objective F(z) = f(T(z)), where f is the original function (embedded below, translated from the MATLAB original) and T acts on the finite
 coordinates only: uniform grid of 151 levels on the interval of the domain, piecewise-affine map with |T - Id| <=
 half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
 Global optimum of F on this domain: f* = 0.3978873577297382 (A: analytic reduction, 1-D root to machine precision); the optimal solution is listed in
 reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
 Variables: x1:K x2:K


Python version of the MATLAB wrapper branin_hoo_k_2D.m: same instance, same values (level indices 1..m as in MATLAB).
Usage:  from branin_hoo_k_2D import branin_hoo_k_2D;  PB = branin_hoo_k_2D();  f = PB.func_idx(PB.x0_idx)
"""
# Copyright (C) 2026 J. F. A. Madeira. LGPL-3.0-or-later.
import math
from types import SimpleNamespace

_NAME = 'branin_hoo_k_2D'
_N = 2
_TYPES = ['K', 'K']
_CARD = [151.0, 151.0]
_LB = [-5.0, 2.128654168571864]
_UB = [7.283648681640625, 15.0]
_LBW = [-5.0, 2.128654168571864]
_UBW = [7.283648681640625, 15.0]
_TB = [-3.116507202148438, 12.254112889295332]
_TI = [-3.141592653589793, 12.275]
_CW = [math.nan, math.nan]
_CX = [math.nan, math.nan]
_M = 151
_PERM = {0: [37, 57, 67, 119, 24, 96, 18, 105, 23, 124, 56, 144, 50, 33, 28, 95, 128, 61, 42, 145, 31, 117, 13, 65, 151, 113, 84, 146, 83, 133, 75, 11, 14, 141, 20, 21, 123, 121, 36, 142, 51, 7, 139, 93, 74, 147, 114, 72, 107, 132, 3, 71, 106, 86, 150, 90, 102, 12, 17, 29, 89, 126, 59, 2, 16, 137, 35, 48, 73, 4, 9, 32, 109, 63, 10, 94, 91, 58, 78, 148, 138, 62, 70, 26, 92, 129, 49, 101, 81, 69, 39, 108, 111, 27, 110, 112, 131, 88, 80, 100, 149, 120, 135, 68, 6, 66, 97, 52, 127, 41, 122, 87, 103, 82, 1, 47, 104, 64, 134, 115, 34, 79, 53, 125, 46, 54, 85, 136, 30, 143, 38, 55, 5, 15, 22, 98, 130, 116, 8, 76, 45, 77, 140, 40, 43, 118, 99, 25, 60, 44, 19], 1: [110, 79, 44, 135, 119, 53, 17, 88, 144, 8, 145, 21, 36, 127, 80, 77, 128, 147, 96, 43, 87, 32, 14, 16, 38, 23, 61, 2, 51, 70, 71, 132, 15, 50, 20, 22, 107, 35, 136, 40, 57, 7, 151, 108, 140, 26, 62, 117, 134, 141, 98, 10, 66, 148, 46, 69, 3, 138, 19, 30, 68, 92, 78, 4, 18, 76, 37, 81, 83, 6, 105, 34, 54, 101, 12, 89, 41, 60, 1, 28, 49, 64, 104, 29, 116, 33, 52, 63, 84, 106, 146, 56, 82, 31, 121, 109, 42, 93, 85, 102, 123, 120, 133, 73, 11, 72, 103, 131, 122, 129, 86, 142, 111, 90, 9, 55, 113, 74, 24, 5, 139, 91, 143, 125, 59, 149, 99, 114, 47, 94, 58, 75, 27, 39, 48, 126, 25, 13, 45, 118, 95, 130, 65, 124, 137, 97, 100, 112, 67, 150, 115]}

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
    return (x[1] - 1.275 * (x[0] / math.pi) ** 2 + 5 * x[0] / math.pi - 6) ** 2 + 10 * (1 - 1 / (8 * math.pi)) * math.cos(x[0]) + 10

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

def branin_hoo_k_2D():
    PB = SimpleNamespace()
    PB.name = _NAME
    PB.nvar = 2
    PB.var_types = ['K', 'K']
    PB.cardinality = [151.0, 151.0]
    PB.f_star = 0.39788735772973816
    PB.f_star_justification = 'A'
    PB.n_minimizers_documented = 2
    PB.description = 'DFO Benchmark Suite v2.1.0 core (Python): K adaptation of branin_hoo (F = f o T).'
    PB.x0_idx = [109.0, 36.0]
    PB.x0_alt_idx = [[32.0, 41.0], [126.0, 50.0], [51.0, 40.0], [90.0, 42.0], [13.0, 47.0], [136.0, 36.0], [61.0, 45.0], [99.0, 47.0], [22.0, 37.0]]
    PB.suite_release = '2.1.0'
    PB.base_function = 'branin_hoo'
    PB.m = 151
    PB.pattern = 'K'
    PB.strategy = 'baseline'
    PB.scale_factors_C = []
    PB.lb_orig = [-5.0, 2.128654168571864]
    PB.ub_orig = [7.283648681640625, 15.0]
    PB.lb_C_work = [-5.0, 2.128654168571864]
    PB.ub_C_work = [7.283648681640625, 15.0]
    PB.ProblemData = _problem_data()
    PB.func_idx = _objective
    PB.func_F = lambda xr: _objective(_cell2idx(xr))
    PB.idx2cell = _idx2cell
    PB.cell2idx = _cell2idx
    PB.decode = lambda xr: _decode(_cell2idx(xr))
    return PB
