"""schaffer2_k_2D  DFO Benchmark Suite v2.1.0 core, Python version (self-contained; standard library only).
SCHAFFER2_K_2D  DFO Benchmark Suite v2.1.0 core: schaffer2, dimension 2, pattern K, strategy baseline.

 Copyright (C) 2026 J. F. A. Madeira

 Objective F(z) = f(T(z)), where f is the original function (embedded below, translated from the MATLAB original) and T acts on the finite
 coordinates only: uniform grid of 151 levels on the interval of the domain, piecewise-affine map with |T - Id| <=
 half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
 Global optimum of F on this domain: f* = 0 (E: exact by structure); the optimal solution is listed in
 reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
 Variables: x1:K x2:K


Python version of the MATLAB wrapper schaffer2_k_2D.m: same instance, same values (level indices 1..m as in MATLAB).
Usage:  from schaffer2_k_2D import schaffer2_k_2D;  PB = schaffer2_k_2D();  f = PB.func_idx(PB.x0_idx)
"""
# Copyright (C) 2026 J. F. A. Madeira. LGPL-3.0-or-later.
import math
from types import SimpleNamespace

_NAME = 'schaffer2_k_2D'
_N = 2
_TYPES = ['K', 'K']
_CARD = [151.0, 151.0]
_LB = [-100.0, -64.70436417212824]
_UB = [65.149169921875, 100.0]
_LBW = [-100.0, -64.70436417212824]
_UBW = [65.149169921875, 100.0]
_TB = [0.19049641927082916, 0.07935240224219342]
_TI = [0.0, 0.0]
_CW = [math.nan, math.nan]
_CX = [math.nan, math.nan]
_M = 151
_PERM = {0: [55, 107, 66, 9, 88, 137, 110, 3, 95, 61, 134, 144, 41, 102, 83, 123, 111, 138, 69, 133, 131, 52, 149, 99, 31, 80, 62, 114, 121, 73, 38, 39, 146, 113, 21, 48, 103, 65, 51, 77, 87, 26, 29, 56, 93, 91, 20, 112, 33, 42, 75, 18, 127, 118, 72, 36, 74, 6, 15, 143, 104, 16, 22, 92, 58, 78, 100, 128, 119, 19, 27, 67, 85, 145, 8, 37, 94, 11, 44, 117, 105, 148, 2, 122, 28, 50, 115, 14, 49, 1, 116, 34, 136, 57, 82, 108, 47, 35, 130, 24, 25, 139, 32, 129, 46, 13, 64, 63, 132, 150, 43, 23, 17, 140, 84, 12, 109, 54, 53, 70, 101, 120, 151, 79, 5, 106, 86, 142, 125, 90, 59, 10, 135, 96, 126, 124, 141, 7, 97, 147, 68, 45, 98, 76, 81, 40, 60, 71, 30, 89, 4], 1: [89, 65, 148, 76, 84, 77, 151, 119, 27, 93, 85, 22, 73, 41, 3, 90, 55, 32, 75, 44, 133, 34, 127, 2, 6, 11, 91, 141, 138, 111, 56, 144, 98, 62, 57, 82, 49, 58, 142, 131, 116, 96, 71, 103, 45, 66, 107, 113, 33, 43, 9, 106, 51, 117, 25, 115, 114, 94, 88, 36, 7, 146, 23, 35, 59, 78, 26, 47, 54, 21, 29, 130, 101, 81, 104, 39, 60, 129, 46, 12, 48, 74, 5, 149, 145, 53, 87, 109, 97, 4, 123, 38, 83, 128, 86, 67, 105, 40, 110, 28, 121, 140, 37, 122, 52, 125, 70, 69, 136, 95, 50, 132, 143, 68, 124, 20, 1, 64, 63, 80, 112, 10, 139, 92, 19, 134, 137, 147, 13, 108, 79, 30, 15, 120, 150, 17, 100, 42, 135, 102, 118, 99, 24, 8, 31, 18, 61, 14, 72, 16, 126]}

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
    return (x[0] ** 2 + x[1] ** 2) ** 0.25 * (math.sin(50 * (x[0] ** 2 + x[1] ** 2) ** 0.1) ** 2 + 1)

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

def schaffer2_k_2D():
    PB = SimpleNamespace()
    PB.name = _NAME
    PB.nvar = 2
    PB.var_types = ['K', 'K']
    PB.cardinality = [151.0, 151.0]
    PB.f_star = 0.0
    PB.f_star_justification = 'E'
    PB.n_minimizers_documented = 1
    PB.description = 'DFO Benchmark Suite v2.1.0 core (Python): K adaptation of schaffer2 (F = f o T).'
    PB.x0_idx = [100.0, 66.0]
    PB.x0_alt_idx = [[26.0, 56.0], [119.0, 60.0], [45.0, 63.0], [81.0, 52.0], [7.0, 61.0], [143.0, 66.0], [69.0, 53.0], [104.0, 57.0], [31.0, 67.0]]
    PB.suite_release = '2.1.0'
    PB.base_function = 'schaffer2'
    PB.m = 151
    PB.pattern = 'K'
    PB.strategy = 'baseline'
    PB.scale_factors_C = []
    PB.lb_orig = [-100.0, -64.70436417212824]
    PB.ub_orig = [65.149169921875, 100.0]
    PB.lb_C_work = [-100.0, -64.70436417212824]
    PB.ub_C_work = [65.149169921875, 100.0]
    PB.ProblemData = _problem_data()
    PB.func_idx = _objective
    PB.func_F = lambda xr: _objective(_cell2idx(xr))
    PB.idx2cell = _idx2cell
    PB.cell2idx = _cell2idx
    PB.decode = lambda xr: _decode(_cell2idx(xr))
    return PB
