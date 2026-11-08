"""schaffer1_k_2D  DFO Benchmark Suite v2.1.0 core, Python version (self-contained; standard library only).
SCHAFFER1_K_2D  DFO Benchmark Suite v2.1.0 core: schaffer1, dimension 2, pattern K, strategy baseline.

 Copyright (C) 2026 J. F. A. Madeira

 Objective F(z) = f(T(z)), where f is the original function (embedded below, translated from the MATLAB original) and T acts on the finite
 coordinates only: uniform grid of 151 levels on the interval of the domain, piecewise-affine map with |T - Id| <=
 half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
 Global optimum of F on this domain: f* = 0 (E: exact by structure); the optimal solution is listed in
 reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
 Variables: x1:K x2:K


Python version of the MATLAB wrapper schaffer1_k_2D.m: same instance, same values (level indices 1..m as in MATLAB).
Usage:  from schaffer1_k_2D import schaffer1_k_2D;  PB = schaffer1_k_2D();  f = PB.func_idx(PB.x0_idx)
"""
# Copyright (C) 2026 J. F. A. Madeira. LGPL-3.0-or-later.
import math
from types import SimpleNamespace

_NAME = 'schaffer1_k_2D'
_N = 2
_TYPES = ['K', 'K']
_CARD = [151.0, 151.0]
_LB = [-72.149169921875, -100.0]
_UB = [100.0, 69.3710308387949]
_LBW = [-72.149169921875, -100.0]
_UBW = [100.0, 69.3710308387949]
_TB = [0.15348144531249375, 0.4934782976849732]
_TI = [0.0, 0.0]
_CW = [math.nan, math.nan]
_CX = [math.nan, math.nan]
_M = 151
_PERM = {0: [150, 42, 55, 129, 71, 80, 52, 35, 25, 119, 60, 58, 133, 118, 1, 31, 139, 5, 56, 114, 120, 98, 36, 43, 32, 51, 87, 123, 89, 82, 135, 100, 7, 41, 131, 76, 27, 77, 2, 6, 132, 19, 22, 29, 78, 137, 115, 146, 81, 117, 105, 8, 40, 75, 13, 26, 128, 116, 4, 79, 102, 70, 108, 20, 44, 63, 10, 142, 37, 85, 73, 99, 97, 96, 65, 18, 38, 69, 126, 68, 127, 141, 62, 45, 148, 24, 57, 74, 21, 130, 112, 3, 47, 23, 48, 28, 107, 95, 50, 83, 84, 54, 124, 109, 103, 67, 11, 145, 140, 30, 94, 72, 64, 106, 125, 53, 33, 147, 88, 101, 9, 15, 134, 91, 136, 111, 86, 34, 104, 61, 16, 93, 143, 17, 39, 12, 144, 110, 46, 113, 49, 138, 66, 59, 121, 151, 90, 92, 122, 14, 149], 1: [117, 81, 128, 63, 87, 43, 129, 35, 44, 28, 116, 85, 113, 29, 102, 68, 124, 46, 33, 75, 110, 119, 57, 23, 104, 127, 12, 107, 115, 151, 98, 150, 8, 53, 139, 18, 56, 1, 3, 138, 47, 137, 90, 74, 38, 59, 78, 16, 112, 148, 82, 9, 134, 66, 14, 86, 131, 54, 62, 120, 118, 109, 10, 22, 45, 83, 11, 31, 39, 5, 106, 52, 13, 123, 67, 20, 136, 80, 25, 71, 126, 93, 65, 48, 147, 94, 60, 77, 24, 143, 76, 6, 95, 122, 51, 32, 108, 2, 55, 133, 89, 58, 96, 69, 4, 73, 17, 140, 145, 37, 146, 130, 72, 132, 149, 61, 42, 141, 97, 111, 19, 27, 30, 103, 26, 125, 100, 144, 121, 79, 36, 114, 88, 41, 64, 40, 91, 7, 84, 15, 99, 50, 135, 142, 70, 101, 92, 34, 21, 105, 49]}

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
    r = math.sqrt(x[0] ** 2 + x[1] ** 2)
    return 0.5 + (math.sin(r) ** 2 - 0.5) / ((1 + 0.001 * (x[0] ** 2 + x[1] ** 2)) ** 2)

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

def schaffer1_k_2D():
    PB = SimpleNamespace()
    PB.name = _NAME
    PB.nvar = 2
    PB.var_types = ['K', 'K']
    PB.cardinality = [151.0, 151.0]
    PB.f_star = 0.0
    PB.f_star_justification = 'E'
    PB.n_minimizers_documented = 1
    PB.description = 'DFO Benchmark Suite v2.1.0 core (Python): K adaptation of schaffer1 (F = f o T).'
    PB.x0_idx = [25.0, 16.0]
    PB.x0_alt_idx = [[119.0, 6.0], [44.0, 10.0], [83.0, 13.0], [6.0, 2.0], [143.0, 11.0], [68.0, 16.0], [106.0, 3.0], [29.0, 7.0], [124.0, 17.0]]
    PB.suite_release = '2.1.0'
    PB.base_function = 'schaffer1'
    PB.m = 151
    PB.pattern = 'K'
    PB.strategy = 'baseline'
    PB.scale_factors_C = []
    PB.lb_orig = [-72.149169921875, -100.0]
    PB.ub_orig = [100.0, 69.3710308387949]
    PB.lb_C_work = [-72.149169921875, -100.0]
    PB.ub_C_work = [100.0, 69.3710308387949]
    PB.ProblemData = _problem_data()
    PB.func_idx = _objective
    PB.func_F = lambda xr: _objective(_cell2idx(xr))
    PB.idx2cell = _idx2cell
    PB.cell2idx = _cell2idx
    PB.decode = lambda xr: _decode(_cell2idx(xr))
    return PB
