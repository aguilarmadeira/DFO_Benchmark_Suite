"""sixhumpcamel_k_2D  DFO Benchmark Suite v2.1.0 core, Python version (self-contained; standard library only).
SIXHUMPCAMEL_K_2D  DFO Benchmark Suite v2.1.0 core: sixhumpcamel, dimension 2, pattern K, strategy baseline.

 Copyright (C) 2026 J. F. A. Madeira

 Objective F(z) = f(T(z)), where f is the original function (embedded below, translated from the MATLAB original) and T acts on the finite
 coordinates only: uniform grid of 151 levels on the interval of the domain, piecewise-affine map with |T - Id| <=
 half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
 Global optimum of F on this domain: f* = -1.031628453488552 (A: analytic reduction, 1-D root to machine precision); the optimal solution is listed in
 reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
 Variables: x1:K x2:K


Python version of the MATLAB wrapper sixhumpcamel_k_2D.m: same instance, same values (level indices 1..m as in MATLAB).
Usage:  from sixhumpcamel_k_2D import sixhumpcamel_k_2D;  PB = sixhumpcamel_k_2D();  f = PB.func_idx(PB.x0_idx)
"""
# Copyright (C) 2026 J. F. A. Madeira. LGPL-3.0-or-later.
import math
from types import SimpleNamespace

_NAME = 'sixhumpcamel_k_2D'
_N = 2
_TYPES = ['K', 'K']
_CARD = [151.0, 151.0]
_LB = [-2.03322509765625, -1.3236620433876949]
_UB = [3.0, 2.0]
_LBW = [-2.03322509765625, -1.3236620433876949]
_UBW = [3.0, 2.0]
_TB = [0.08072944335937482, -0.7032451286219918]
_TI = [0.089842, -0.712656]
_CW = [math.nan, math.nan]
_CX = [math.nan, math.nan]
_M = 151
_PERM = {0: [145, 150, 39, 49, 14, 140, 75, 129, 124, 34, 26, 109, 70, 144, 77, 6, 55, 95, 80, 57, 19, 79, 42, 30, 108, 17, 127, 48, 113, 7, 103, 18, 23, 40, 71, 68, 62, 56, 126, 119, 122, 59, 72, 13, 52, 36, 83, 147, 137, 128, 96, 138, 32, 37, 8, 22, 118, 51, 81, 67, 76, 66, 15, 112, 110, 9, 24, 46, 106, 25, 115, 123, 41, 105, 148, 58, 134, 38, 74, 92, 107, 31, 45, 141, 98, 111, 125, 151, 35, 114, 20, 44, 1, 116, 10, 94, 121, 82, 139, 89, 104, 78, 16, 100, 61, 131, 3, 12, 73, 69, 120, 28, 43, 101, 33, 102, 93, 88, 133, 29, 86, 87, 5, 135, 65, 142, 84, 21, 63, 97, 54, 47, 130, 4, 99, 53, 11, 50, 136, 2, 143, 60, 90, 146, 132, 64, 117, 27, 149, 91, 85], 1: [103, 143, 53, 49, 141, 87, 27, 42, 142, 15, 71, 20, 58, 28, 64, 73, 55, 128, 149, 61, 122, 130, 3, 72, 12, 31, 65, 80, 59, 8, 126, 133, 91, 125, 34, 1, 62, 45, 132, 30, 54, 100, 77, 95, 113, 116, 46, 147, 105, 22, 39, 51, 136, 84, 119, 23, 98, 52, 5, 70, 109, 10, 123, 75, 69, 11, 25, 102, 56, 26, 135, 7, 137, 63, 134, 106, 94, 120, 110, 47, 6, 33, 92, 89, 85, 32, 68, 88, 38, 121, 151, 48, 117, 129, 14, 104, 90, 138, 57, 93, 144, 83, 21, 2, 66, 60, 9, 18, 79, 76, 36, 139, 50, 108, 41, 111, 101, 97, 107, 40, 96, 127, 17, 112, 78, 131, 99, 37, 145, 148, 74, 67, 19, 29, 124, 82, 146, 86, 35, 44, 140, 114, 16, 150, 81, 43, 118, 115, 13, 4, 24]}

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
    return (4 - 2.1 * x[0] ** 2 + x[0] ** 4 / 3) * x[0] ** 2 + x[0] * x[1] + (4 * x[1] ** 2 - 4) * x[1] ** 2

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

def sixhumpcamel_k_2D():
    PB = SimpleNamespace()
    PB.name = _NAME
    PB.nvar = 2
    PB.var_types = ['K', 'K']
    PB.cardinality = [151.0, 151.0]
    PB.f_star = -1.0316284534885518
    PB.f_star_justification = 'A'
    PB.n_minimizers_documented = 2
    PB.description = 'DFO Benchmark Suite v2.1.0 core (Python): K adaptation of sixhumpcamel (F = f o T).'
    PB.x0_idx = [71.0, 52.0]
    PB.x0_alt_idx = [[110.0, 61.0], [34.0, 65.0], [129.0, 53.0], [53.0, 58.0], [92.0, 67.0], [16.0, 56.0], [139.0, 59.0], [62.0, 63.0], [101.0, 52.0]]
    PB.suite_release = '2.1.0'
    PB.base_function = 'sixhumpcamel'
    PB.m = 151
    PB.pattern = 'K'
    PB.strategy = 'baseline'
    PB.scale_factors_C = []
    PB.lb_orig = [-2.03322509765625, -1.3236620433876949]
    PB.ub_orig = [3.0, 2.0]
    PB.lb_C_work = [-2.03322509765625, -1.3236620433876949]
    PB.ub_C_work = [3.0, 2.0]
    PB.ProblemData = _problem_data()
    PB.func_idx = _objective
    PB.func_F = lambda xr: _objective(_cell2idx(xr))
    PB.idx2cell = _idx2cell
    PB.cell2idx = _cell2idx
    PB.decode = lambda xr: _decode(_cell2idx(xr))
    return PB
