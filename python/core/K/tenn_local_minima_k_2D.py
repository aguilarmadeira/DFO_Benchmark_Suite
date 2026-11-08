"""tenn_local_minima_k_2D  DFO Benchmark Suite v2.1.0 core, Python version (self-contained; standard library only).
TENN_LOCAL_MINIMA_K_2D  DFO Benchmark Suite v2.1.0 core: tenn_local_minima, dimension 2, pattern K, strategy baseline.

 Copyright (C) 2026 J. F. A. Madeira

 Objective F(z) = f(T(z)), where f is the original function (embedded below, translated from the MATLAB original) and T acts on the finite
 coordinates only: uniform grid of 151 levels on the interval of the domain, piecewise-affine map with |T - Id| <=
 half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
 Global optimum of F on this domain: f* = 2.355817157679958e-31 (E: exact by structure); the optimal solution is listed in
 reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
 Variables: x1:K x2:K


Python version of the MATLAB wrapper tenn_local_minima_k_2D.m: same instance, same values (level indices 1..m as in MATLAB).
Usage:  from tenn_local_minima_k_2D import tenn_local_minima_k_2D;  PB = tenn_local_minima_k_2D();  f = PB.func_idx(PB.x0_idx)
"""
# Copyright (C) 2026 J. F. A. Madeira. LGPL-3.0-or-later.
import math
from types import SimpleNamespace

_NAME = 'tenn_local_minima_k_2D'
_N = 2
_TYPES = ['K', 'K']
_CARD = [151.0, 151.0]
_LB = [-10.0, -10.0]
_UB = [6.6461669921875, 6.877569476197734]
_LBW = [-10.0, -10.0]
_UBW = [6.6461669921875, 6.877569476197734]
_TB = [0.9864702148437505, 1.0266787244491873]
_TI = [1.0, 1.0]
_CW = [math.nan, math.nan]
_CX = [math.nan, math.nan]
_M = 151
_PERM = {0: [95, 58, 68, 37, 146, 65, 33, 106, 114, 25, 100, 138, 81, 111, 34, 42, 74, 80, 125, 77, 6, 66, 71, 49, 50, 132, 90, 46, 82, 20, 89, 69, 140, 10, 51, 115, 27, 97, 4, 75, 15, 91, 62, 61, 143, 85, 86, 131, 36, 92, 113, 110, 9, 13, 39, 107, 7, 119, 137, 134, 122, 48, 55, 93, 31, 53, 1, 116, 83, 102, 139, 136, 78, 41, 96, 40, 63, 148, 150, 29, 70, 17, 151, 105, 72, 22, 59, 79, 147, 88, 23, 47, 141, 127, 18, 8, 2, 98, 73, 11, 24, 3, 126, 28, 101, 103, 54, 142, 149, 26, 118, 104, 14, 76, 16, 94, 123, 99, 112, 128, 121, 35, 135, 64, 133, 45, 84, 130, 129, 56, 52, 60, 43, 144, 87, 120, 38, 19, 12, 124, 145, 67, 44, 108, 5, 117, 30, 109, 21, 57, 32], 1: [18, 67, 69, 94, 114, 28, 46, 7, 87, 146, 63, 19, 39, 147, 53, 68, 45, 76, 80, 137, 124, 126, 117, 131, 122, 38, 116, 128, 36, 135, 30, 90, 64, 142, 92, 129, 130, 79, 145, 57, 16, 115, 103, 35, 44, 120, 134, 52, 61, 144, 123, 109, 85, 14, 41, 141, 93, 101, 40, 2, 37, 139, 125, 97, 89, 138, 106, 86, 83, 5, 136, 58, 127, 43, 10, 42, 133, 23, 107, 32, 118, 20, 149, 91, 75, 25, 62, 82, 104, 143, 151, 148, 9, 88, 22, 12, 6, 4, 78, 15, 150, 8, 50, 110, 1, 119, 60, 73, 27, 33, 105, 111, 21, 84, 24, 102, 96, 108, 3, 66, 11, 47, 100, 77, 56, 59, 99, 71, 17, 74, 72, 81, 65, 132, 113, 13, 70, 54, 49, 26, 55, 121, 112, 48, 98, 95, 31, 51, 34, 29, 140]}

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
    d = len(x)
    s = 0.0
    for i in range(d - 1): s += ((x[i] - 1) ** 2) * (1 + 10 * math.sin(math.pi * x[i + 1]) ** 2)
    f = math.pi / d * s
    return f + math.pi / d * (10 * math.sin(math.pi * x[0]) ** 2 + (x[d - 1] - 1) ** 2)

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

def tenn_local_minima_k_2D():
    PB = SimpleNamespace()
    PB.name = _NAME
    PB.nvar = 2
    PB.var_types = ['K', 'K']
    PB.cardinality = [151.0, 151.0]
    PB.f_star = 2.3558171576799584e-31
    PB.f_star_justification = 'E'
    PB.n_minimizers_documented = 1
    PB.description = 'DFO Benchmark Suite v2.1.0 core (Python): K adaptation of tenn_local_minima (F = f o T).'
    PB.x0_idx = [78.0, 69.0]
    PB.x0_alt_idx = [[5.0, 78.0], [147.0, 82.0], [73.0, 70.0], [109.0, 74.0], [35.0, 84.0], [129.0, 73.0], [55.0, 75.0], [90.0, 80.0], [16.0, 69.0]]
    PB.suite_release = '2.1.0'
    PB.base_function = 'tenn_local_minima'
    PB.m = 151
    PB.pattern = 'K'
    PB.strategy = 'baseline'
    PB.scale_factors_C = []
    PB.lb_orig = [-10.0, -10.0]
    PB.ub_orig = [6.6461669921875, 6.877569476197734]
    PB.lb_C_work = [-10.0, -10.0]
    PB.ub_C_work = [6.6461669921875, 6.877569476197734]
    PB.ProblemData = _problem_data()
    PB.func_idx = _objective
    PB.func_F = lambda xr: _objective(_cell2idx(xr))
    PB.idx2cell = _idx2cell
    PB.cell2idx = _cell2idx
    PB.decode = lambda xr: _decode(_cell2idx(xr))
    return PB
