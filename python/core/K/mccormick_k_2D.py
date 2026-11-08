"""mccormick_k_2D  DFO Benchmark Suite v2.1.0 core, Python version (self-contained; standard library only).
MCCORMICK_K_2D  DFO Benchmark Suite v2.1.0 core: mccormick, dimension 2, pattern K, strategy baseline.

 Copyright (C) 2026 J. F. A. Madeira

 Objective F(z) = f(T(z)), where f is the original function (embedded below, translated from the MATLAB original) and T acts on the finite
 coordinates only: uniform grid of 151 levels on the interval of the domain, piecewise-affine map with |T - Id| <=
 half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
 Global optimum of F on this domain: f* = -1.913222954981037 (L: literature value at a refined point); the optimal solution is listed in
 reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
 Variables: x1:K x2:K


Python version of the MATLAB wrapper mccormick_k_2D.m: same instance, same values (level indices 1..m as in MATLAB).
Usage:  from mccormick_k_2D import mccormick_k_2D;  PB = mccormick_k_2D();  f = PB.func_idx(PB.x0_idx)
"""
# Copyright (C) 2026 J. F. A. Madeira. LGPL-3.0-or-later.
import math
from types import SimpleNamespace

_NAME = 'mccormick_k_2D'
_N = 2
_TYPES = ['K', 'K']
_CARD = [151.0, 151.0]
_LB = [-1.5, -3.0]
_UB = [3.1183013916015625, 2.0500198140527357]
_LBW = [-1.5, -3.0]
_UBW = [3.1183013916015625, 2.0500198140527357]
_TB = [-0.545551045735677, -1.5523276533048824]
_TI = [-0.5471975537167932, -1.547197552675593]
_CW = [math.nan, math.nan]
_CX = [math.nan, math.nan]
_M = 151
_PERM = {0: [26, 78, 130, 131, 2, 79, 100, 48, 99, 137, 144, 123, 83, 138, 113, 51, 62, 121, 116, 14, 44, 92, 111, 133, 39, 103, 81, 17, 40, 24, 38, 128, 108, 9, 105, 122, 147, 73, 60, 37, 67, 134, 10, 102, 107, 132, 41, 45, 90, 143, 25, 34, 77, 11, 31, 141, 50, 8, 16, 33, 66, 7, 96, 94, 88, 59, 3, 118, 112, 82, 65, 87, 119, 71, 28, 109, 70, 19, 47, 6, 36, 58, 125, 95, 148, 110, 27, 42, 72, 4, 91, 18, 149, 142, 43, 12, 80, 104, 106, 145, 101, 150, 93, 126, 120, 86, 23, 124, 52, 127, 29, 89, 49, 69, 53, 54, 20, 117, 98, 129, 46, 13, 84, 5, 114, 146, 151, 1, 57, 135, 68, 76, 32, 115, 35, 64, 75, 139, 30, 15, 74, 22, 63, 136, 140, 55, 21, 97, 56, 61, 85], 1: [48, 37, 89, 121, 20, 41, 15, 122, 64, 145, 126, 10, 27, 32, 4, 6, 80, 51, 43, 36, 112, 127, 77, 148, 55, 81, 50, 141, 78, 88, 1, 86, 79, 95, 12, 14, 68, 31, 102, 151, 99, 53, 119, 70, 24, 129, 42, 97, 71, 114, 94, 111, 109, 105, 150, 65, 52, 116, 17, 146, 115, 9, 13, 23, 44, 85, 101, 147, 117, 59, 120, 137, 113, 73, 143, 130, 83, 106, 49, 118, 39, 93, 72, 46, 82, 2, 30, 45, 75, 7, 125, 21, 103, 26, 47, 128, 84, 104, 110, 19, 134, 138, 98, 60, 63, 136, 29, 11, 58, 91, 35, 123, 57, 142, 135, 62, 28, 8, 108, 139, 140, 25, 132, 18, 3, 149, 100, 16, 74, 133, 87, 96, 54, 5, 61, 92, 107, 66, 67, 56, 124, 144, 131, 76, 90, 34, 22, 40, 69, 33, 38]}

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
    return math.sin(x[0] + x[1]) + (x[0] - x[1]) ** 2 - 3 / 2 * x[0] + 5 / 2 * x[1] + 1

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

def mccormick_k_2D():
    PB = SimpleNamespace()
    PB.name = _NAME
    PB.nvar = 2
    PB.var_types = ['K', 'K']
    PB.cardinality = [151.0, 151.0]
    PB.f_star = -1.9132229549810367
    PB.f_star_justification = 'L'
    PB.n_minimizers_documented = 1
    PB.description = 'DFO Benchmark Suite v2.1.0 core (Python): K adaptation of mccormick (F = f o T).'
    PB.x0_idx = [57.0, 72.0]
    PB.x0_alt_idx = [[96.0, 76.0], [20.0, 80.0], [115.0, 68.0], [39.0, 77.0], [78.0, 82.0], [2.0, 71.0], [150.0, 73.0], [76.0, 83.0], [111.0, 72.0]]
    PB.suite_release = '2.1.0'
    PB.base_function = 'mccormick'
    PB.m = 151
    PB.pattern = 'K'
    PB.strategy = 'baseline'
    PB.scale_factors_C = []
    PB.lb_orig = [-1.5, -3.0]
    PB.ub_orig = [3.1183013916015625, 2.0500198140527357]
    PB.lb_C_work = [-1.5, -3.0]
    PB.ub_C_work = [3.1183013916015625, 2.0500198140527357]
    PB.ProblemData = _problem_data()
    PB.func_idx = _objective
    PB.func_F = lambda xr: _objective(_cell2idx(xr))
    PB.idx2cell = _idx2cell
    PB.cell2idx = _cell2idx
    PB.decode = lambda xr: _decode(_cell2idx(xr))
    return PB
