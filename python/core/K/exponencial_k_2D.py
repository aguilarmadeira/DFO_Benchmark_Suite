"""exponencial_k_2D  DFO Benchmark Suite v2.1.0 core, Python version (self-contained; standard library only).
EXPONENCIAL_K_2D  DFO Benchmark Suite v2.1.0 core: exponencial, dimension 2, pattern K, strategy baseline.

 Copyright (C) 2026 J. F. A. Madeira

 Objective F(z) = f(T(z)), where f is the original function (embedded below, translated from the MATLAB original) and T acts on the finite
 coordinates only: uniform grid of 151 levels on the interval of the domain, piecewise-affine map with |T - Id| <=
 half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
 Global optimum of F on this domain: f* = -1 (E: exact by structure); the optimal solution is listed in
 reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
 Variables: x1:K x2:K


Python version of the MATLAB wrapper exponencial_k_2D.m: same instance, same values (level indices 1..m as in MATLAB).
Usage:  from exponencial_k_2D import exponencial_k_2D;  PB = exponencial_k_2D();  f = PB.func_idx(PB.x0_idx)
"""
# Copyright (C) 2026 J. F. A. Madeira. LGPL-3.0-or-later.
import math
from types import SimpleNamespace

_NAME = 'exponencial_k_2D'
_N = 2
_TYPES = ['K', 'K']
_CARD = [151.0, 151.0]
_LB = [-0.64875732421875, -0.6332164812274552]
_UB = [1.0, 1.0]
_LBW = [-0.64875732421875, -0.6332164812274552]
_UBW = [1.0, 1.0]
_TB = [-0.0002461100260416904, -0.0017061084861724751]
_TI = [0.0, 0.0]
_CW = [math.nan, math.nan]
_CX = [math.nan, math.nan]
_M = 151
_PERM = {0: [135, 145, 151, 36, 23, 58, 87, 32, 128, 142, 57, 109, 117, 132, 131, 107, 111, 28, 103, 9, 66, 13, 88, 138, 22, 115, 108, 96, 79, 98, 16, 17, 67, 51, 56, 76, 35, 74, 49, 129, 50, 144, 21, 26, 59, 41, 124, 148, 133, 119, 24, 95, 53, 99, 18, 61, 86, 45, 149, 6, 27, 34, 71, 73, 70, 11, 20, 33, 89, 127, 146, 100, 136, 101, 10, 134, 43, 125, 8, 38, 116, 80, 7, 130, 81, 123, 4, 91, 147, 37, 118, 30, 60, 104, 31, 82, 102, 97, 3, 25, 93, 40, 65, 126, 77, 114, 90, 54, 92, 72, 5, 29, 150, 55, 110, 122, 12, 15, 94, 106, 141, 63, 68, 78, 113, 14, 85, 62, 140, 52, 121, 46, 83, 112, 44, 47, 64, 120, 139, 75, 48, 42, 105, 19, 143, 69, 39, 2, 84, 137, 1], 1: [59, 33, 121, 57, 23, 76, 39, 17, 24, 85, 139, 148, 136, 14, 134, 68, 5, 56, 87, 114, 30, 128, 46, 77, 43, 26, 93, 38, 83, 11, 16, 88, 67, 51, 149, 127, 122, 6, 73, 110, 50, 117, 22, 72, 52, 123, 53, 115, 108, 15, 118, 32, 105, 101, 91, 106, 129, 47, 54, 99, 28, 112, 135, 146, 119, 13, 21, 35, 151, 81, 107, 126, 49, 62, 111, 143, 144, 132, 131, 41, 150, 125, 9, 55, 84, 113, 7, 94, 42, 147, 133, 95, 64, 138, 34, 86, 44, 4, 8, 29, 109, 45, 130, 19, 82, 61, 137, 60, 141, 79, 12, 36, 120, 63, 3, 71, 20, 25, 103, 140, 58, 74, 80, 90, 1, 27, 100, 78, 145, 70, 10, 66, 104, 2, 69, 75, 142, 18, 37, 116, 98, 96, 31, 102, 92, 48, 40, 89, 97, 65, 124]}

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
    s = 0.0
    for v in x: s += v ** 2
    return -math.exp(-0.5 * s)

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

def exponencial_k_2D():
    PB = SimpleNamespace()
    PB.name = _NAME
    PB.nvar = 2
    PB.var_types = ['K', 'K']
    PB.cardinality = [151.0, 151.0]
    PB.f_star = -1.0
    PB.f_star_justification = 'E'
    PB.n_minimizers_documented = 1
    PB.description = 'DFO Benchmark Suite v2.1.0 core (Python): K adaptation of exponencial (F = f o T).'
    PB.x0_idx = [104.0, 125.0]
    PB.x0_alt_idx = [[27.0, 130.0], [122.0, 119.0], [46.0, 128.0], [85.0, 131.0], [8.0, 121.0], [146.0, 125.0], [70.0, 134.0], [109.0, 122.0], [32.0, 126.0]]
    PB.suite_release = '2.1.0'
    PB.base_function = 'exponencial'
    PB.m = 151
    PB.pattern = 'K'
    PB.strategy = 'baseline'
    PB.scale_factors_C = []
    PB.lb_orig = [-0.64875732421875, -0.6332164812274552]
    PB.ub_orig = [1.0, 1.0]
    PB.lb_C_work = [-0.64875732421875, -0.6332164812274552]
    PB.ub_C_work = [1.0, 1.0]
    PB.ProblemData = _problem_data()
    PB.func_idx = _objective
    PB.func_F = lambda xr: _objective(_cell2idx(xr))
    PB.idx2cell = _idx2cell
    PB.cell2idx = _cell2idx
    PB.decode = lambda xr: _decode(_cell2idx(xr))
    return PB
