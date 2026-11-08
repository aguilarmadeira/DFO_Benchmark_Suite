"""salomon_ck_5D  DFO Benchmark Suite v2.1.0 core, Python version (self-contained; standard library only).
SALOMON_CK_5D  DFO Benchmark Suite v2.1.0 core: salomon, dimension 5, pattern CK, strategy moderate.

 Copyright (C) 2026 J. F. A. Madeira

 Objective F(z) = f(T(z)), where f is the original function (embedded below, translated from the MATLAB original) and T acts on the finite
 coordinates only: uniform grid of 151 levels on the interval of the domain, piecewise-affine map with |T - Id| <=
 half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
 Global optimum of F on this domain: f* = 0 (E: exact by structure); the optimal solution is listed in
 reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
 Variables: x1:C x2:K x3:C x4:K x5:C


Python version of the MATLAB wrapper salomon_ck_5D.m: same instance, same values (level indices 1..m as in MATLAB).
Usage:  from salomon_ck_5D import salomon_ck_5D;  PB = salomon_ck_5D();  f = PB.func_idx(PB.x0_idx)
"""
# Copyright (C) 2026 J. F. A. Madeira. LGPL-3.0-or-later.
import math
from types import SimpleNamespace

_NAME = 'salomon_ck_5D'
_N = 5
_TYPES = ['C', 'K', 'C', 'K', 'C']
_CARD = [math.inf, 151.0, math.inf, 151.0, math.inf]
_LB = [-70.399169921875, -100.0, -66.579328, -100.0, -71.18680418004234]
_UB = [100.0, 66.2599197276838, 100.0, 72.41232819658475, 100.0]
_LBW = [-70.399169921875, -100.0, -6657.9328, -100.0, -0.7118680418004235]
_UBW = [100.0, 66.2599197276838, 10000.0, 72.41232819658475, 1.0]
_TB = [math.nan, -0.2440481633897349, math.nan, -0.0008496459808355894, math.nan]
_TI = [math.nan, 0.0, math.nan, 0.0, math.nan]
_CW = [0.0, math.nan, 0.0, math.nan, 0.0]
_CX = [0.0, math.nan, 0.0, math.nan, 0.0]
_M = 151
_PERM = {1: [116, 20, 59, 151, 101, 142, 94, 71, 88, 72, 79, 107, 114, 52, 45, 106, 73, 99, 124, 38, 150, 10, 105, 61, 113, 32, 89, 126, 111, 37, 98, 26, 64, 2, 42, 7, 16, 80, 75, 34, 122, 82, 67, 115, 23, 66, 44, 55, 109, 112, 33, 143, 21, 135, 119, 90, 14, 93, 125, 147, 30, 85, 47, 149, 132, 36, 49, 129, 103, 43, 51, 104, 50, 68, 146, 96, 1, 31, 131, 128, 63, 87, 18, 84, 41, 60, 127, 22, 56, 5, 145, 35, 95, 54, 148, 123, 39, 24, 133, 117, 4, 137, 8, 83, 12, 121, 130, 15, 69, 29, 91, 65, 118, 92, 144, 28, 6, 139, 48, 53, 76, 141, 57, 3, 40, 9, 100, 27, 138, 19, 81, 134, 102, 17, 25, 86, 74, 77, 140, 120, 62, 97, 70, 108, 13, 110, 136, 78, 58, 46, 11], 3: [141, 77, 104, 69, 15, 99, 40, 145, 103, 48, 89, 128, 96, 36, 106, 88, 79, 26, 60, 140, 92, 38, 13, 30, 74, 41, 136, 59, 101, 70, 25, 114, 112, 109, 7, 9, 44, 144, 107, 118, 137, 130, 49, 17, 53, 126, 46, 127, 115, 19, 78, 116, 117, 124, 51, 8, 100, 34, 76, 3, 135, 43, 50, 94, 91, 39, 52, 4, 82, 122, 55, 22, 54, 151, 98, 62, 102, 35, 119, 32, 149, 90, 129, 5, 121, 66, 133, 29, 142, 11, 20, 110, 95, 61, 86, 134, 47, 33, 83, 16, 14, 84, 18, 93, 23, 148, 31, 27, 123, 42, 105, 80, 68, 108, 138, 45, 24, 71, 125, 73, 97, 85, 81, 28, 146, 37, 2, 132, 111, 147, 120, 58, 12, 65, 75, 6, 1, 10, 56, 143, 21, 63, 150, 131, 57, 67, 139, 87, 113, 64, 72]}

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
    r = math.sqrt(s)
    return 1 - math.cos(2 * math.pi * r) + 0.1 * r

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

def salomon_ck_5D():
    PB = SimpleNamespace()
    PB.name = _NAME
    PB.nvar = 5
    PB.var_types = ['C', 'K', 'C', 'K', 'C']
    PB.cardinality = [math.inf, 151.0, math.inf, 151.0, math.inf]
    PB.f_star = 0.0
    PB.f_star_justification = 'E'
    PB.n_minimizers_documented = 1
    PB.description = 'DFO Benchmark Suite v2.1.0 core (Python): CK adaptation of salomon (F = f o T).'
    PB.x0_idx = [-22.39120066165924, 94.0, 2222.744334336, 133.0, 0.6283021306684282]
    PB.x0_alt_idx = [[21.539835333824158, 98.0, -5373.273022464, 139.0, 0.09197683132148837], [-64.32537138462067, 88.0, 5580.983586816, 146.0, -0.29001025238604417], [91.4301198720932, 90.0, -2015.0337699839993, 132.0, 0.8713848203004944], [4.56648051738739, 100.0, 8811.289915392, 138.0, 0.4893977365929616], [48.49751651287079, 89.0, 1215.2725585919989, 145.0, -0.06095830963028781], [-37.367690205574036, 93.0, -3822.086320128, 148.0, -0.44294539333782024], [70.46303451061249, 96.0, 4573.511811072002, 134.0, 0.7184496793487183], [-16.899821162223816, 85.0, -463.84706764800103, 140.0, 0.3364625956411855], [27.031214833259583, 94.0, 7931.751063552, 147.0, -0.1998627037057542]]
    PB.suite_release = '2.1.0'
    PB.base_function = 'salomon'
    PB.m = 151
    PB.pattern = 'CK'
    PB.strategy = 'moderate'
    PB.scale_factors_C = [1.0, 100.0, 0.01]
    PB.lb_orig = [-70.399169921875, -100.0, -66.579328, -100.0, -71.18680418004234]
    PB.ub_orig = [100.0, 66.2599197276838, 100.0, 72.41232819658475, 100.0]
    PB.lb_C_work = [-70.399169921875, -100.0, -6657.9328, -100.0, -0.7118680418004235]
    PB.ub_C_work = [100.0, 66.2599197276838, 10000.0, 72.41232819658475, 1.0]
    PB.ProblemData = _problem_data()
    PB.func_idx = _objective
    PB.func_F = lambda xr: _objective(_cell2idx(xr))
    PB.idx2cell = _idx2cell
    PB.cell2idx = _cell2idx
    PB.decode = lambda xr: _decode(_cell2idx(xr))
    return PB
