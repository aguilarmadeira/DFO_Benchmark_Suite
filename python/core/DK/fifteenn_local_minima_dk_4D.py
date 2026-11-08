"""fifteenn_local_minima_dk_4D  DFO Benchmark Suite v2.1.0 core, Python version (self-contained; standard library only).
FIFTEENN_LOCAL_MINIMA_DK_4D  DFO Benchmark Suite v2.1.0 core: fifteenn_local_minima, dimension 4, pattern DK, strategy baseline.

 Copyright (C) 2026 J. F. A. Madeira

 Objective F(z) = f(T(z)), where f is the original function (embedded below, translated from the MATLAB original) and T acts on the finite
 coordinates only: uniform grid of 25 levels on the interval of the domain, piecewise-affine map with |T - Id| <=
 half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
 Global optimum of F on this domain: f* = 1.349783804395672e-32 (E: exact by structure); the optimal solution is listed in
 reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
 Variables: x1:D x2:K x3:D x4:K


Python version of the MATLAB wrapper fifteenn_local_minima_dk_4D.m: same instance, same values (level indices 1..m as in MATLAB).
Usage:  from fifteenn_local_minima_dk_4D import fifteenn_local_minima_dk_4D;  PB = fifteenn_local_minima_dk_4D();  f = PB.func_idx(PB.x0_idx)
"""
# Copyright (C) 2026 J. F. A. Madeira. LGPL-3.0-or-later.
import math
from types import SimpleNamespace

_NAME = 'fifteenn_local_minima_dk_4D'
_N = 4
_TYPES = ['D', 'K', 'D', 'K']
_CARD = [25.0, 25.0, 25.0, 25.0]
_LB = [-10.0, -10.0, -10.0, -10.0]
_UB = [7.275073242187499, 6.176609256718996, 6.949132799999999, 6.95551853394419]
_LBW = [-10.0, -10.0, -10.0, -10.0]
_UBW = [7.275073242187499, 6.176609256718996, 6.949132799999999, 6.95551853394419]
_TB = [0.7969207763671875, 0.7844061711459975, 1.299421866666668, 1.30367902262946]
_TI = [1.0, 1.0, 1.0, 1.0]
_CW = [math.nan, math.nan, math.nan, math.nan]
_CX = [math.nan, math.nan, math.nan, math.nan]
_M = 25
_PERM = {1: [3, 17, 6, 13, 8, 9, 20, 7, 22, 18, 14, 11, 5, 10, 2, 24, 19, 1, 21, 16, 15, 4, 23, 25, 12], 3: [3, 19, 7, 10, 24, 6, 11, 2, 16, 21, 18, 17, 12, 5, 13, 20, 14, 25, 1, 8, 15, 23, 9, 22, 4]}

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
    for i in range(d - 1): s += ((x[i] - 1) ** 2) * (1 + 10 * math.sin(3 * math.pi * x[i + 1]) ** 2)
    f = 0.1 * s
    return f + 0.1 * (math.sin(3 * math.pi * x[0]) ** 2 + (x[d - 1] - 1) ** 2 * (1 + math.sin(2 * math.pi * x[d - 1]) ** 2))

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

def fifteenn_local_minima_dk_4D():
    PB = SimpleNamespace()
    PB.name = _NAME
    PB.nvar = 4
    PB.var_types = ['D', 'K', 'D', 'K']
    PB.cardinality = [25.0, 25.0, 25.0, 25.0]
    PB.f_star = 1.3497838043956716e-32
    PB.f_star_justification = 'E'
    PB.n_minimizers_documented = 1
    PB.description = 'DFO Benchmark Suite v2.1.0 core (Python): DK adaptation of fifteenn_local_minima (F = f o T).'
    PB.x0_idx = [4.0, 24.0, 9.0, 7.0]
    PB.x0_alt_idx = [[23.0, 24.0, 21.0, 5.0], [11.0, 23.0, 14.0, 6.0], [17.0, 24.0, 3.0, 7.0], [5.0, 24.0, 18.0, 8.0], [20.0, 23.0, 7.0, 5.0], [8.0, 23.0, 23.0, 6.0], [14.0, 25.0, 12.0, 7.0], [2.0, 23.0, 5.0, 5.0], [24.0, 24.0, 17.0, 6.0]]
    PB.suite_release = '2.1.0'
    PB.base_function = 'fifteenn_local_minima'
    PB.m = 25
    PB.pattern = 'DK'
    PB.strategy = 'baseline'
    PB.scale_factors_C = []
    PB.lb_orig = [-10.0, -10.0, -10.0, -10.0]
    PB.ub_orig = [7.275073242187499, 6.176609256718996, 6.949132799999999, 6.95551853394419]
    PB.lb_C_work = [-10.0, -10.0, -10.0, -10.0]
    PB.ub_C_work = [7.275073242187499, 6.176609256718996, 6.949132799999999, 6.95551853394419]
    PB.ProblemData = _problem_data()
    PB.func_idx = _objective
    PB.func_F = lambda xr: _objective(_cell2idx(xr))
    PB.idx2cell = _idx2cell
    PB.cell2idx = _cell2idx
    PB.decode = lambda xr: _decode(_cell2idx(xr))
    return PB
