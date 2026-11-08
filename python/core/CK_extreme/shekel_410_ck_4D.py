"""shekel_410_ck_4D  DFO Benchmark Suite v2.1.0 core, Python version (self-contained; standard library only).
SHEKEL_410_CK_4D  DFO Benchmark Suite v2.1.0 core: shekel_410, dimension 4, pattern CK, strategy extreme.

 Copyright (C) 2026 J. F. A. Madeira

 Objective F(z) = f(T(z)), where f is the original function (embedded below, translated from the MATLAB original) and T acts on the finite
 coordinates only: uniform grid of 151 levels on the interval of the domain, piecewise-affine map with |T - Id| <=
 half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
 Global optimum of F on this domain: f* = -10.53640981669204 (L: literature value at a refined point); the optimal solution is listed in
 reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
 Variables: x1:C x2:K x3:C x4:K


Python version of the MATLAB wrapper shekel_410_ck_4D.m: same instance, same values (level indices 1..m as in MATLAB).
Usage:  from shekel_410_ck_4D import shekel_410_ck_4D;  PB = shekel_410_ck_4D();  f = PB.func_idx(PB.x0_idx)
"""
# Copyright (C) 2026 J. F. A. Madeira. LGPL-3.0-or-later.
import math
from types import SimpleNamespace

_NAME = 'shekel_410_ck_4D'
_N = 4
_TYPES = ['C', 'K', 'C', 'K']
_CARD = [math.inf, 151.0, math.inf, 151.0]
_LB = [1.69879150390625, 1.7686226693085405, 1.9790336000000002, 1.3936693044564765]
_UB = [10.0, 10.0, 10.0, 10.0]
_LBW = [1.69879150390625, 1.7686226693085405, 197903360.00000006, 1.3936693044564765]
_UBW = [10.0, 10.0, 1000000000.0, 10.0]
_TB = [math.nan, 4.018532473030873, math.nan, 3.9755685131195335]
_TI = [math.nan, 4.000592937456021, math.nan, 3.999509802482246]
_CW = [4.000746531113532, math.nan, 399966339.77086973, math.nan]
_CX = [4.000746531113532, math.nan, 3.999663397708697, math.nan]
_M = 151
_PERM = {1: [147, 23, 126, 149, 3, 29, 119, 62, 123, 49, 116, 58, 26, 146, 89, 54, 102, 78, 82, 121, 41, 105, 11, 83, 109, 125, 77, 36, 122, 70, 127, 98, 32, 52, 151, 48, 80, 115, 90, 88, 59, 150, 14, 22, 30, 44, 6, 9, 128, 92, 45, 129, 37, 42, 12, 60, 34, 53, 66, 64, 81, 135, 106, 25, 79, 5, 18, 76, 148, 16, 24, 65, 140, 141, 10, 40, 113, 96, 51, 19, 112, 84, 15, 110, 124, 117, 21, 38, 131, 95, 50, 73, 27, 103, 33, 17, 143, 93, 56, 132, 130, 108, 7, 87, 134, 142, 67, 71, 20, 1, 68, 57, 137, 111, 31, 86, 72, 39, 43, 75, 118, 63, 85, 138, 114, 100, 99, 107, 101, 97, 133, 91, 120, 35, 94, 139, 136, 46, 61, 28, 47, 13, 4, 69, 8, 144, 104, 74, 2, 55, 145], 3: [72, 142, 42, 144, 130, 118, 151, 64, 67, 38, 70, 126, 91, 59, 99, 97, 94, 50, 109, 104, 41, 53, 36, 68, 149, 133, 34, 81, 147, 116, 26, 61, 127, 146, 139, 48, 49, 113, 93, 90, 52, 114, 137, 23, 73, 148, 77, 128, 108, 30, 115, 136, 39, 44, 75, 122, 37, 56, 84, 112, 60, 10, 125, 131, 54, 8, 22, 43, 74, 19, 28, 69, 32, 85, 14, 88, 66, 21, 55, 24, 119, 129, 20, 7, 145, 78, 27, 45, 143, 120, 57, 102, 35, 123, 40, 106, 12, 4, 65, 100, 5, 80, 17, 98, 46, 25, 79, 134, 33, 15, 82, 124, 76, 13, 47, 103, 89, 58, 62, 95, 18, 87, 110, 71, 16, 2, 1, 11, 141, 3, 29, 150, 31, 83, 9, 138, 63, 117, 135, 111, 6, 121, 140, 92, 51, 101, 107, 105, 86, 132, 96]}

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

_A = [[4, 4, 4, 4, 0.1], [1, 1, 1, 1, 0.2], [8, 8, 8, 8, 0.2], [6, 6, 6, 6, 0.4], [3, 7, 3, 7, 0.4],
      [2, 9, 2, 9, 0.6], [5, 5, 3, 3, 0.3], [8, 1, 8, 1, 0.7], [6, 2, 6, 2, 0.5], [7, 3.6, 7, 3.6, 0.5]]
def _orig(x):
    s = 0.0
    for r in _A[:10]:
        e = 0.0
        for j in range(4): e += (x[j] - r[j]) ** 2
        s += 1 / (e + r[4])
    return -s

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

def shekel_410_ck_4D():
    PB = SimpleNamespace()
    PB.name = _NAME
    PB.nvar = 4
    PB.var_types = ['C', 'K', 'C', 'K']
    PB.cardinality = [math.inf, 151.0, math.inf, 151.0]
    PB.f_star = -10.536409816692045
    PB.f_star_justification = 'L'
    PB.n_minimizers_documented = 1
    PB.description = 'DFO Benchmark Suite v2.1.0 core (Python): CK adaptation of shekel_410 (F = f o T).'
    PB.x0_idx = [6.631687372922897, 133.0, 978439642.3168001, 136.0]
    PB.x0_alt_idx = [[2.5459363162517548, 122.0, 612683574.4768001, 142.0], [9.25824162364006, 127.0, 216126995.66080007, 149.0], [5.172490566968918, 129.0, 774386257.1008, 135.0], [7.16672620177269, 119.0, 371669576.0896001, 141.0], [3.0809751451015472, 128.0, 929928837.5296001, 148.0], [8.236803859472275, 132.0, 533372258.71360004, 130.0], [4.151052802801132, 120.0, 321618745.75360006, 137.0], [6.1209684908390045, 124.0, 695074941.3376001, 143.0], [2.035217434167862, 133.0, 483321428.3776, 150.0]]
    PB.suite_release = '2.1.0'
    PB.base_function = 'shekel_410'
    PB.m = 151
    PB.pattern = 'CK'
    PB.strategy = 'extreme'
    PB.scale_factors_C = [1.0, 100000000.0]
    PB.lb_orig = [1.69879150390625, 1.7686226693085405, 1.9790336000000002, 1.3936693044564765]
    PB.ub_orig = [10.0, 10.0, 10.0, 10.0]
    PB.lb_C_work = [1.69879150390625, 1.7686226693085405, 197903360.00000006, 1.3936693044564765]
    PB.ub_C_work = [10.0, 10.0, 1000000000.0, 10.0]
    PB.ProblemData = _problem_data()
    PB.func_idx = _objective
    PB.func_F = lambda xr: _objective(_cell2idx(xr))
    PB.idx2cell = _idx2cell
    PB.cell2idx = _cell2idx
    PB.decode = lambda xr: _decode(_cell2idx(xr))
    return PB
