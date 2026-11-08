"""gulf_k_3D  DFO Benchmark Suite v2.1.0 core, Python version (self-contained; standard library only).
GULF_K_3D  DFO Benchmark Suite v2.1.0 core: gulf, dimension 3, pattern K, strategy baseline.

 Copyright (C) 2026 J. F. A. Madeira

 Objective F(z) = f(T(z)), where f is the original function (embedded below, translated from the MATLAB original) and T acts on the finite
 coordinates only: uniform grid of 41 levels on the interval of the domain, piecewise-affine map with |T - Id| <=
 half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
 Global optimum of F on this domain: f* = 9.709124424797534e-31 (E: exact by structure); the optimal solution is listed in
 reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
 Variables: x1:K x2:K x3:K


Python version of the MATLAB wrapper gulf_k_3D.m: same instance, same values (level indices 1..m as in MATLAB).
Usage:  from gulf_k_3D import gulf_k_3D;  PB = gulf_k_3D();  f = PB.func_idx(PB.x0_idx)
"""
# Copyright (C) 2026 J. F. A. Madeira. LGPL-3.0-or-later.
import math
from types import SimpleNamespace

_NAME = 'gulf_k_3D'
_N = 3
_TYPES = ['K', 'K', 'K']
_CARD = [41.0, 41.0, 41.0]
_LB = [0.1, 4.185989534115735, 0.0]
_UB = [86.60752209472656, 25.6, 4.0132832]
_LBW = [0.1, 4.185989534115735, 0.0]
_UBW = [86.60752209472656, 25.6, 4.0132832]
_TB = [49.84182520446778, 25.064649738352895, 1.5049812]
_TI = [50.0, 25.0, 1.5]
_CW = [math.nan, math.nan, math.nan]
_CX = [math.nan, math.nan, math.nan]
_M = 41
_PERM = {0: [31, 23, 3, 25, 16, 15, 2, 7, 6, 29, 14, 17, 27, 18, 32, 30, 39, 28, 40, 19, 12, 9, 41, 8, 34, 13, 33, 20, 1, 11, 24, 4, 5, 37, 26, 38, 10, 21, 36, 22, 35], 1: [26, 7, 31, 22, 29, 5, 3, 8, 32, 39, 23, 34, 21, 37, 6, 24, 35, 19, 28, 38, 25, 12, 17, 13, 15, 41, 11, 27, 9, 33, 4, 36, 20, 30, 14, 16, 40, 1, 2, 10, 18], 2: [9, 36, 25, 1, 32, 5, 33, 26, 38, 19, 11, 35, 7, 12, 8, 31, 41, 40, 14, 4, 18, 30, 21, 27, 20, 24, 17, 6, 16, 29, 34, 28, 3, 15, 2, 37, 10, 23, 13, 39, 22]}

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
    for i in range(1, 100):
        u = 25 + (-50 * math.log(0.01 * i)) ** (2 / 3)
        s += (math.exp(-(u - x[1]) ** x[2] / x[0]) - 0.01 * i) ** 2
    return s

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

def gulf_k_3D():
    PB = SimpleNamespace()
    PB.name = _NAME
    PB.nvar = 3
    PB.var_types = ['K', 'K', 'K']
    PB.cardinality = [41.0, 41.0, 41.0]
    PB.f_star = 9.709124424797534e-31
    PB.f_star_justification = 'E'
    PB.n_minimizers_documented = 1
    PB.description = 'DFO Benchmark Suite v2.1.0 core (Python): K adaptation of gulf (F = f o T).'
    PB.x0_idx = [1.0, 27.0, 40.0]
    PB.x0_alt_idx = [[40.0, 24.0, 21.0], [21.0, 25.0, 2.0], [31.0, 28.0, 29.0], [11.0, 24.0, 9.0], [35.0, 25.0, 37.0], [16.0, 27.0, 17.0], [26.0, 23.0, 7.0], [6.0, 26.0, 25.0], [38.0, 27.0, 15.0]]
    PB.suite_release = '2.1.0'
    PB.base_function = 'gulf'
    PB.m = 41
    PB.pattern = 'K'
    PB.strategy = 'baseline'
    PB.scale_factors_C = []
    PB.lb_orig = [0.1, 4.185989534115735, 0.0]
    PB.ub_orig = [86.60752209472656, 25.6, 4.0132832]
    PB.lb_C_work = [0.1, 4.185989534115735, 0.0]
    PB.ub_C_work = [86.60752209472656, 25.6, 4.0132832]
    PB.ProblemData = _problem_data()
    PB.func_idx = _objective
    PB.func_F = lambda xr: _objective(_cell2idx(xr))
    PB.idx2cell = _idx2cell
    PB.cell2idx = _cell2idx
    PB.decode = lambda xr: _decode(_cell2idx(xr))
    return PB
