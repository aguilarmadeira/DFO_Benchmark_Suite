"""becker_lago_k_2D  DFO Benchmark Suite v2.1.0 core, Python version (self-contained; standard library only).
BECKER_LAGO_K_2D  DFO Benchmark Suite v2.1.0 core: becker_lago, dimension 2, pattern K, strategy baseline.

 Copyright (C) 2026 J. F. A. Madeira

 Objective F(z) = f(T(z)), where f is the original function (embedded below, translated from the MATLAB original) and T acts on the finite
 coordinates only: uniform grid of 151 levels on the interval of the domain, piecewise-affine map with |T - Id| <=
 half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
 Global optimum of F on this domain: f* = 0 (E: exact by structure); the optimal solution is listed in
 reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
 Variables: x1:K x2:K


Python version of the MATLAB wrapper becker_lago_k_2D.m: same instance, same values (level indices 1..m as in MATLAB).
Usage:  from becker_lago_k_2D import becker_lago_k_2D;  PB = becker_lago_k_2D();  f = PB.func_idx(PB.x0_idx)
"""
# Copyright (C) 2026 J. F. A. Madeira. LGPL-3.0-or-later.
import math
from types import SimpleNamespace

_NAME = 'becker_lago_k_2D'
_N = 2
_TYPES = ['K', 'K']
_CARD = [151.0, 151.0]
_LB = [-6.2031982421875, -6.85068333079307]
_UB = [10.0, 10.0]
_LBW = [-6.2031982421875, -6.85068333079307]
_UBW = [10.0, 10.0]
_TB = [-5.014963704427084, -5.05327710884181]
_TI = [-5.0, -5.0]
_CW = [math.nan, math.nan]
_CX = [math.nan, math.nan]
_M = 151
_PERM = {0: [51, 29, 150, 135, 120, 30, 20, 146, 49, 60, 92, 55, 66, 11, 107, 40, 95, 41, 86, 3, 110, 121, 34, 44, 12, 43, 100, 38, 129, 74, 99, 138, 4, 141, 46, 94, 15, 114, 109, 50, 56, 35, 103, 108, 104, 132, 93, 87, 106, 26, 57, 72, 136, 10, 23, 31, 79, 71, 76, 88, 22, 28, 32, 39, 53, 126, 6, 137, 21, 42, 47, 70, 78, 27, 48, 151, 101, 130, 36, 142, 122, 96, 24, 97, 2, 149, 82, 5, 25, 9, 67, 123, 105, 54, 148, 102, 62, 147, 1, 133, 77, 16, 37, 81, 116, 111, 98, 52, 83, 33, 117, 73, 127, 125, 89, 144, 134, 18, 84, 58, 91, 113, 80, 69, 61, 64, 145, 143, 115, 17, 85, 90, 19, 124, 118, 8, 112, 139, 140, 45, 7, 119, 75, 14, 131, 128, 13, 68, 63, 59, 65], 1: [125, 83, 70, 18, 89, 62, 45, 117, 132, 33, 109, 111, 53, 103, 147, 47, 108, 143, 112, 121, 4, 60, 63, 130, 13, 51, 114, 128, 123, 40, 56, 126, 81, 88, 46, 59, 16, 61, 100, 146, 142, 136, 79, 52, 122, 10, 134, 113, 116, 118, 119, 43, 135, 148, 25, 32, 37, 48, 54, 6, 69, 30, 34, 133, 55, 107, 7, 73, 23, 44, 49, 95, 101, 29, 50, 67, 2, 96, 38, 64, 140, 131, 26, 65, 127, 115, 85, 137, 28, 12, 71, 84, 98, 141, 72, 20, 66, 17, 5, 24, 124, 21, 42, 86, 19, 75, 104, 58, 90, 39, 76, 138, 93, 120, 144, 22, 151, 27, 94, 68, 102, 129, 92, 139, 74, 77, 99, 91, 3, 35, 149, 110, 41, 15, 9, 36, 8, 145, 82, 87, 57, 150, 1, 97, 80, 106, 14, 11, 105, 31, 78]}

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
    return (abs(x[0]) - 5) ** 2 + (abs(x[1]) - 5) ** 2

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

def becker_lago_k_2D():
    PB = SimpleNamespace()
    PB.name = _NAME
    PB.nvar = 2
    PB.var_types = ['K', 'K']
    PB.cardinality = [151.0, 151.0]
    PB.f_star = 0.0
    PB.f_star_justification = 'E'
    PB.n_minimizers_documented = 4
    PB.description = 'DFO Benchmark Suite v2.1.0 core (Python): K adaptation of becker_lago (F = f o T).'
    PB.x0_idx = [127.0, 70.0]
    PB.x0_alt_idx = [[50.0, 74.0], [89.0, 83.0], [13.0, 73.0], [136.0, 76.0], [60.0, 80.0], [99.0, 69.0], [23.0, 79.0], [118.0, 81.0], [41.0, 70.0]]
    PB.suite_release = '2.1.0'
    PB.base_function = 'becker_lago'
    PB.m = 151
    PB.pattern = 'K'
    PB.strategy = 'baseline'
    PB.scale_factors_C = []
    PB.lb_orig = [-6.2031982421875, -6.85068333079307]
    PB.ub_orig = [10.0, 10.0]
    PB.lb_C_work = [-6.2031982421875, -6.85068333079307]
    PB.ub_C_work = [10.0, 10.0]
    PB.ProblemData = _problem_data()
    PB.func_idx = _objective
    PB.func_F = lambda xr: _objective(_cell2idx(xr))
    PB.idx2cell = _idx2cell
    PB.cell2idx = _cell2idx
    PB.decode = lambda xr: _decode(_cell2idx(xr))
    return PB
