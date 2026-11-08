"""griewank_ck_5D  DFO Benchmark Suite v2.1.0 core, Python version (self-contained; standard library only).
GRIEWANK_CK_5D  DFO Benchmark Suite v2.1.0 core: griewank, dimension 5, pattern CK, strategy moderate.

 Copyright (C) 2026 J. F. A. Madeira

 Objective F(z) = f(T(z)), where f is the original function (embedded below, translated from the MATLAB original) and T acts on the finite
 coordinates only: uniform grid of 151 levels on the interval of the domain, piecewise-affine map with |T - Id| <=
 half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
 Global optimum of F on this domain: f* = 0 (E: exact by structure); the optimal solution is listed in
 reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
 Variables: x1:C x2:K x3:C x4:K x5:C


Python version of the MATLAB wrapper griewank_ck_5D.m: same instance, same values (level indices 1..m as in MATLAB).
Usage:  from griewank_ck_5D import griewank_ck_5D;  PB = griewank_ck_5D();  f = PB.func_idx(PB.x0_idx)
"""
# Copyright (C) 2026 J. F. A. Madeira. LGPL-3.0-or-later.
import math
from types import SimpleNamespace

_NAME = 'griewank_ck_5D'
_N = 5
_TYPES = ['C', 'K', 'C', 'K', 'C']
_CARD = [math.inf, 151.0, math.inf, 151.0, math.inf]
_LB = [-600.0, -600.0, -379.987968, -600.0, -600.0]
_UB = [363.00439453125, 431.78174058832496, 600.0, 415.6168263223657, 405.6001639232293]
_LBW = [-600.0, -600.0, -37998.7968, -600.0, -6.0]
_UBW = [363.00439453125, 431.78174058832496, 60000.0, 415.6168263223657, 4.056001639232292]
_TB = [math.nan, -1.5665904587714294, math.nan, 2.599316951270339, math.nan]
_TI = [math.nan, 0.0, math.nan, 0.0, math.nan]
_CW = [0.0, math.nan, 0.0, math.nan, 0.0]
_CX = [0.0, math.nan, 0.0, math.nan, 0.0]
_M = 151
_PERM = {1: [43, 32, 102, 37, 41, 46, 15, 90, 105, 60, 98, 130, 125, 79, 66, 73, 26, 83, 13, 70, 38, 74, 71, 62, 6, 147, 51, 93, 21, 64, 114, 1, 4, 8, 82, 12, 19, 28, 30, 59, 39, 3, 116, 77, 136, 27, 34, 36, 84, 2, 112, 97, 108, 129, 151, 122, 31, 81, 86, 110, 33, 42, 47, 89, 103, 25, 85, 52, 149, 87, 115, 53, 7, 20, 119, 94, 11, 135, 61, 17, 120, 128, 146, 44, 78, 76, 16, 29, 138, 134, 111, 75, 142, 72, 92, 56, 118, 109, 10, 133, 121, 69, 124, 150, 40, 131, 88, 63, 144, 9, 58, 126, 49, 57, 23, 140, 80, 132, 113, 137, 14, 24, 68, 99, 54, 104, 100, 65, 50, 106, 48, 148, 35, 123, 141, 117, 45, 96, 127, 5, 139, 18, 55, 107, 95, 145, 67, 22, 91, 143, 101], 3: [68, 101, 54, 69, 148, 33, 103, 27, 37, 111, 58, 107, 32, 140, 98, 8, 96, 78, 146, 105, 62, 116, 67, 3, 43, 137, 15, 141, 94, 74, 1, 151, 135, 9, 114, 119, 86, 145, 31, 47, 118, 95, 7, 106, 18, 102, 88, 134, 46, 4, 77, 24, 143, 90, 17, 71, 34, 122, 147, 144, 36, 45, 83, 59, 14, 150, 39, 55, 61, 20, 26, 57, 11, 123, 53, 108, 16, 115, 66, 131, 129, 70, 81, 50, 84, 92, 110, 121, 132, 75, 104, 82, 25, 125, 142, 113, 30, 117, 19, 44, 23, 133, 112, 120, 51, 60, 100, 76, 10, 22, 72, 6, 64, 73, 40, 21, 128, 41, 13, 65, 35, 48, 93, 124, 80, 5, 149, 97, 85, 12, 87, 139, 79, 136, 42, 38, 109, 29, 63, 89, 99, 126, 49, 130, 138, 91, 56, 52, 2, 28, 127]}

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
    n = len(x)
    s = 0.0
    for v in x: s += v ** 2
    p = 1.0
    for i in range(n): p *= math.cos(x[i] / math.sqrt(i + 1))
    return 1 + 1 / 4000 * s - p

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

def griewank_ck_5D():
    PB = SimpleNamespace()
    PB.name = _NAME
    PB.nvar = 5
    PB.var_types = ['C', 'K', 'C', 'K', 'C']
    PB.cardinality = [math.inf, 151.0, math.inf, 151.0, math.inf]
    PB.f_star = 0.0
    PB.f_star_justification = 'E'
    PB.n_minimizers_documented = 1
    PB.description = 'DFO Benchmark Suite v2.1.0 core (Python): CK adaptation of griewank (F = f o T).'
    PB.x0_idx = [362.0639605522156, 48.0, 36982.04260761601, 44.0, -5.524020959224917]
    PB.x0_alt_idx = [[-133.07452940940857, 37.0, -7705.408733183991, 26.0, 2.2049720362180842], [115.20004105567932, 41.0, 56738.60004249601, 33.0, -0.9455559003407936], [-370.0638921260834, 51.0, 12051.148701695995, 39.0, -4.0960838368996715], [239.33732628822327, 39.0, -37152.087195648, 25.0, 3.632909158543331], [-251.56921076774597, 43.0, 31055.07537715201, 31.0, 0.3999606088632257], [-3.294640302658081, 47.0, -17395.529760767997, 35.0, -1.8439405833621612], [-488.5585734844208, 34.0, 50811.63281203201, 41.0, -4.994468519921038], [301.40596890449524, 44.0, 2361.0276741120033, 27.0, 2.7345244755219635], [-192.32187008857727, 48.0, -23510.654681087995, 34.0, -0.4160034610369152]]
    PB.suite_release = '2.1.0'
    PB.base_function = 'griewank'
    PB.m = 151
    PB.pattern = 'CK'
    PB.strategy = 'moderate'
    PB.scale_factors_C = [1.0, 100.0, 0.01]
    PB.lb_orig = [-600.0, -600.0, -379.987968, -600.0, -600.0]
    PB.ub_orig = [363.00439453125, 431.78174058832496, 600.0, 415.6168263223657, 405.6001639232293]
    PB.lb_C_work = [-600.0, -600.0, -37998.7968, -600.0, -6.0]
    PB.ub_C_work = [363.00439453125, 431.78174058832496, 60000.0, 415.6168263223657, 4.056001639232292]
    PB.ProblemData = _problem_data()
    PB.func_idx = _objective
    PB.func_F = lambda xr: _objective(_cell2idx(xr))
    PB.idx2cell = _idx2cell
    PB.cell2idx = _cell2idx
    PB.decode = lambda xr: _decode(_cell2idx(xr))
    return PB
