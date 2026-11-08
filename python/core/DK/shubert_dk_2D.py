"""shubert_dk_2D  DFO Benchmark Suite v2.1.0 core, Python version (self-contained; standard library only).
SHUBERT_DK_2D  DFO Benchmark Suite v2.1.0 core: shubert, dimension 2, pattern DK, strategy baseline.

 Copyright (C) 2026 J. F. A. Madeira

 Objective F(z) = f(T(z)), where f is the original function (embedded below, translated from the MATLAB original) and T acts on the finite
 coordinates only: uniform grid of 151 levels on the interval of the domain, piecewise-affine map with |T - Id| <=
 half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
 Global optimum of F on this domain: f* = -186.7309088310238 (L: literature value at a refined point); the optimal solution is listed in
 reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
 Variables: x1:D x2:K


Python version of the MATLAB wrapper shubert_dk_2D.m: same instance, same values (level indices 1..m as in MATLAB).
Usage:  from shubert_dk_2D import shubert_dk_2D;  PB = shubert_dk_2D();  f = PB.func_idx(PB.x0_idx)
"""
# Copyright (C) 2026 J. F. A. Madeira. LGPL-3.0-or-later.
import math
from types import SimpleNamespace

_NAME = 'shubert_dk_2D'
_N = 2
_TYPES = ['D', 'K']
_CARD = [151.0, 151.0]
_LB = [-10.0, -6.3071991058273635]
_UB = [7.1274169921874995, 10.0]
_LBW = [-10.0, -6.3071991058273635]
_UBW = [7.1274169921874995, 10.0]
_TB = [-7.031247721354166, 4.89041094684076]
_TI = [-7.0835064072, 4.8580568773]
_CW = [math.nan, math.nan]
_CX = [math.nan, math.nan]
_M = 151
_PERM = {1: [72, 132, 92, 39, 49, 86, 104, 84, 47, 3, 69, 54, 81, 100, 41, 82, 118, 6, 31, 120, 88, 59, 142, 130, 5, 38, 16, 106, 129, 29, 73, 140, 67, 19, 23, 26, 141, 131, 13, 110, 24, 91, 126, 55, 7, 149, 58, 61, 116, 80, 101, 90, 114, 18, 43, 60, 10, 30, 95, 83, 150, 46, 53, 89, 28, 48, 63, 97, 137, 64, 2, 151, 8, 96, 117, 22, 123, 127, 34, 134, 42, 70, 1, 74, 35, 65, 143, 144, 71, 125, 146, 76, 113, 115, 40, 111, 108, 107, 109, 12, 112, 98, 36, 138, 139, 62, 122, 147, 77, 68, 27, 21, 32, 87, 136, 78, 66, 50, 57, 102, 133, 124, 44, 148, 93, 135, 99, 14, 37, 56, 119, 105, 33, 4, 85, 128, 145, 15, 75, 9, 103, 20, 94, 121, 25, 52, 17, 11, 51, 79, 45]}

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
    s1 = 0.0
    for i in range(1, 6): s1 += i * math.cos((i + 1) * x[0] + i)
    s2 = 0.0
    for i in range(1, 6): s2 += i * math.cos((i + 1) * x[1] + i)
    return s1 * s2

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

def shubert_dk_2D():
    PB = SimpleNamespace()
    PB.name = _NAME
    PB.nvar = 2
    PB.var_types = ['D', 'K']
    PB.cardinality = [151.0, 151.0]
    PB.f_star = -186.73090883102384
    PB.f_star_justification = 'L'
    PB.n_minimizers_documented = 12
    PB.description = 'DFO Benchmark Suite v2.1.0 core (Python): DK adaptation of shubert (F = f o T).'
    PB.x0_idx = [34.0, 149.0]
    PB.x0_alt_idx = [[129.0, 139.0], [53.0, 143.0], [92.0, 146.0], [15.0, 135.0], [138.0, 145.0], [63.0, 149.0], [102.0, 136.0], [25.0, 141.0], [120.0, 150.0]]
    PB.suite_release = '2.1.0'
    PB.base_function = 'shubert'
    PB.m = 151
    PB.pattern = 'DK'
    PB.strategy = 'baseline'
    PB.scale_factors_C = []
    PB.lb_orig = [-10.0, -6.3071991058273635]
    PB.ub_orig = [7.1274169921874995, 10.0]
    PB.lb_C_work = [-10.0, -6.3071991058273635]
    PB.ub_C_work = [7.1274169921874995, 10.0]
    PB.ProblemData = _problem_data()
    PB.func_idx = _objective
    PB.func_F = lambda xr: _objective(_cell2idx(xr))
    PB.idx2cell = _idx2cell
    PB.cell2idx = _cell2idx
    PB.decode = lambda xr: _decode(_cell2idx(xr))
    return PB
