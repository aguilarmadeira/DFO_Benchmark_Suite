"""tenn_local_minima_k_6D  DFO Benchmark Suite v2.1.0 core, Python version (self-contained; standard library only).
TENN_LOCAL_MINIMA_K_6D  DFO Benchmark Suite v2.1.0 core: tenn_local_minima, dimension 6, pattern K, strategy baseline.

 Copyright (C) 2026 J. F. A. Madeira

 Objective F(z) = f(T(z)), where f is the original function (embedded below, translated from the MATLAB original) and T acts on the finite
 coordinates only: uniform grid of 17 levels on the interval of the domain, piecewise-affine map with |T - Id| <=
 half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
 Global optimum of F on this domain: f* = 7.852723858933194e-32 (E: exact by structure); the optimal solution is listed in
 reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
 Variables: x1:K x2:K x3:K x4:K x5:K x6:K


Python version of the MATLAB wrapper tenn_local_minima_k_6D.m: same instance, same values (level indices 1..m as in MATLAB).
Usage:  from tenn_local_minima_k_6D import tenn_local_minima_k_6D;  PB = tenn_local_minima_k_6D();  f = PB.func_idx(PB.x0_idx)
"""
# Copyright (C) 2026 J. F. A. Madeira. LGPL-3.0-or-later.
import math
from types import SimpleNamespace

_NAME = 'tenn_local_minima_k_6D'
_N = 6
_TYPES = ['K', 'K', 'K', 'K', 'K', 'K']
_CARD = [17.0, 17.0, 17.0, 17.0, 17.0, 17.0]
_LB = [-10.0, -10.0, -10.0, -10.0, -10.0, -10.0]
_UB = [6.2961669921875, 7.188680587308845, 6.1987328, 6.780008329862557, 6.343473806433986, 7.315738244459228]
_LBW = [-10.0, -10.0, -10.0, -10.0, -10.0, -10.0]
_UBW = [6.2961669921875, 7.188680587308845, 6.1987328, 6.780008329862557, 6.343473806433986, 7.315738244459228]
_TB = [1.2036148071289077, 0.7429253670680289, 1.1366288000000004, 0.48750520616409787, 1.2361382419233635, 0.8223364027870161]
_TI = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
_CW = [math.nan, math.nan, math.nan, math.nan, math.nan, math.nan]
_CX = [math.nan, math.nan, math.nan, math.nan, math.nan, math.nan]
_M = 17
_PERM = {0: [6, 4, 13, 3, 11, 1, 17, 2, 15, 10, 16, 9, 12, 7, 14, 5, 8], 1: [7, 10, 9, 4, 17, 15, 13, 5, 2, 6, 14, 8, 12, 1, 3, 16, 11], 2: [15, 13, 17, 1, 10, 5, 14, 11, 16, 2, 8, 6, 12, 9, 4, 7, 3], 3: [11, 7, 6, 9, 13, 1, 10, 3, 2, 15, 5, 14, 12, 17, 8, 4, 16], 4: [15, 16, 8, 3, 7, 14, 5, 9, 17, 4, 1, 2, 11, 12, 10, 13, 6], 5: [15, 3, 2, 6, 7, 4, 13, 16, 1, 9, 8, 12, 17, 5, 14, 10, 11]}

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
    d = len(x)
    s = 0.0
    for i in range(d - 1): s += ((x[i] - 1) ** 2) * (1 + 10 * math.sin(math.pi * x[i + 1]) ** 2)
    f = math.pi / d * s
    return f + math.pi / d * (10 * math.sin(math.pi * x[0]) ** 2 + (x[d - 1] - 1) ** 2)

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

def tenn_local_minima_k_6D():
    PB = SimpleNamespace()
    PB.name = _NAME
    PB.nvar = 6
    PB.var_types = ['K', 'K', 'K', 'K', 'K', 'K']
    PB.cardinality = [17.0, 17.0, 17.0, 17.0, 17.0, 17.0]
    PB.f_star = 7.852723858933194e-32
    PB.f_star_justification = 'E'
    PB.n_minimizers_documented = 1
    PB.description = 'DFO Benchmark Suite v2.1.0 core (Python): K adaptation of tenn_local_minima (F = f o T).'
    PB.x0_idx = [13.0, 5.0, 15.0, 7.0, 7.0, 9.0]
    PB.x0_alt_idx = [[5.0, 6.0, 11.0, 8.0, 2.0, 7.0], [9.0, 6.0, 3.0, 6.0, 14.0, 4.0], [1.0, 5.0, 14.0, 7.0, 9.0, 2.0], [17.0, 5.0, 6.0, 7.0, 5.0, 15.0], [9.0, 6.0, 17.0, 6.0, 16.0, 12.0], [13.0, 5.0, 9.0, 6.0, 12.0, 10.0], [5.0, 5.0, 1.0, 7.0, 7.0, 8.0], [15.0, 6.0, 12.0, 8.0, 4.0, 5.0], [7.0, 5.0, 4.0, 6.0, 15.0, 3.0]]
    PB.suite_release = '2.1.0'
    PB.base_function = 'tenn_local_minima'
    PB.m = 17
    PB.pattern = 'K'
    PB.strategy = 'baseline'
    PB.scale_factors_C = []
    PB.lb_orig = [-10.0, -10.0, -10.0, -10.0, -10.0, -10.0]
    PB.ub_orig = [6.2961669921875, 7.188680587308845, 6.1987328, 6.780008329862557, 6.343473806433986, 7.315738244459228]
    PB.lb_C_work = [-10.0, -10.0, -10.0, -10.0, -10.0, -10.0]
    PB.ub_C_work = [6.2961669921875, 7.188680587308845, 6.1987328, 6.780008329862557, 6.343473806433986, 7.315738244459228]
    PB.ProblemData = _problem_data()
    PB.func_idx = _objective
    PB.func_F = lambda xr: _objective(_cell2idx(xr))
    PB.idx2cell = _idx2cell
    PB.cell2idx = _cell2idx
    PB.decode = lambda xr: _decode(_cell2idx(xr))
    return PB
