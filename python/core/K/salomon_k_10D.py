"""salomon_k_10D  DFO Benchmark Suite v2.1.0 core, Python version (self-contained; standard library only).
SALOMON_K_10D  DFO Benchmark Suite v2.1.0 core: salomon, dimension 10, pattern K, strategy baseline.

 Copyright (C) 2026 J. F. A. Madeira

 Objective F(z) = f(T(z)), where f is the original function (embedded below, translated from the MATLAB original) and T acts on the finite
 coordinates only: uniform grid of 9 levels on the interval of the domain, piecewise-affine map with |T - Id| <=
 half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
 Global optimum of F on this domain: f* = 0 (E: exact by structure); the optimal solution is listed in
 reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
 Variables: x1:K x2:K x3:K x4:K x5:K x6:K x7:K x8:K x9:K x10:K


Python version of the MATLAB wrapper salomon_k_10D.m: same instance, same values (level indices 1..m as in MATLAB).
Usage:  from salomon_k_10D import salomon_k_10D;  PB = salomon_k_10D();  f = PB.func_idx(PB.x0_idx)
"""
# Copyright (C) 2026 J. F. A. Madeira. LGPL-3.0-or-later.
import math
from types import SimpleNamespace

_NAME = 'salomon_k_10D'
_N = 10
_TYPES = ['K', 'K', 'K', 'K', 'K', 'K', 'K', 'K', 'K', 'K']
_CARD = [9.0, 9.0, 9.0, 9.0, 9.0, 9.0, 9.0, 9.0, 9.0, 9.0]
_LB = [-100.0, -100.0, -63.779328, -70.41232819658475, -69.91407690731508, -62.55383214873429, -100.0, -100.0, -100.0, -100.0]
_UB = [63.399169921875, 61.59325306101712, 100.0, 100.0, 100.0, 100.0, 69.5715808000383, 64.0719147336193, 70.75630804635489, 61.75997375866169]
_LBW = [-100.0, -100.0, -63.779328, -70.41232819658475, -69.91407690731508, -62.55383214873429, -100.0, -100.0, -100.0, -100.0]
_UBW = [63.399169921875, 61.59325306101712, 100.0, 100.0, 100.0, 100.0, 69.5715808000383, 64.0719147336193, 70.75630804635489, 61.75997375866169]
_TB = [2.124481201171875, 0.9957831631356981, -2.362079999999999, -6.507705122865474, -6.196298067071929, -1.5961450929589347, 5.982238000023941, 2.544946708512043, 6.722692528971805, 1.0999835991635507]
_TI = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
_CW = [math.nan, math.nan, math.nan, math.nan, math.nan, math.nan, math.nan, math.nan, math.nan, math.nan]
_CX = [math.nan, math.nan, math.nan, math.nan, math.nan, math.nan, math.nan, math.nan, math.nan, math.nan]
_M = 9
_PERM = {0: [7, 3, 1, 2, 8, 9, 4, 6, 5], 1: [4, 2, 7, 1, 8, 3, 9, 5, 6], 2: [8, 5, 9, 1, 4, 7, 6, 3, 2], 3: [7, 6, 1, 5, 4, 3, 2, 8, 9], 4: [6, 7, 5, 9, 8, 1, 2, 4, 3], 5: [5, 8, 2, 3, 9, 7, 4, 1, 6], 6: [3, 5, 9, 8, 4, 2, 6, 7, 1], 7: [3, 1, 5, 2, 4, 7, 9, 6, 8], 8: [3, 6, 5, 8, 9, 7, 2, 1, 4], 9: [1, 8, 3, 9, 4, 6, 2, 7, 5]}

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

def salomon_k_10D():
    PB = SimpleNamespace()
    PB.name = _NAME
    PB.nvar = 10
    PB.var_types = ['K', 'K', 'K', 'K', 'K', 'K', 'K', 'K', 'K', 'K']
    PB.cardinality = [9.0, 9.0, 9.0, 9.0, 9.0, 9.0, 9.0, 9.0, 9.0, 9.0]
    PB.f_star = 0.0
    PB.f_star_justification = 'E'
    PB.n_minimizers_documented = 1
    PB.description = 'DFO Benchmark Suite v2.1.0 core (Python): K adaptation of salomon (F = f o T).'
    PB.x0_idx = [7.0, 9.0, 7.0, 1.0, 8.0, 3.0, 4.0, 9.0, 4.0, 9.0]
    PB.x0_alt_idx = [[3.0, 9.0, 3.0, 2.0, 5.0, 2.0, 2.0, 3.0, 2.0, 2.0], [5.0, 8.0, 8.0, 2.0, 4.0, 9.0, 7.0, 5.0, 8.0, 4.0], [1.0, 8.0, 5.0, 1.0, 1.0, 7.0, 5.0, 8.0, 6.0, 5.0], [9.0, 9.0, 2.0, 2.0, 7.0, 6.0, 3.0, 3.0, 4.0, 6.0], [5.0, 8.0, 6.0, 2.0, 5.0, 5.0, 8.0, 5.0, 2.0, 8.0], [7.0, 9.0, 4.0, 1.0, 3.0, 4.0, 6.0, 8.0, 7.0, 1.0], [3.0, 9.0, 8.0, 1.0, 8.0, 3.0, 3.0, 2.0, 5.0, 3.0], [8.0, 8.0, 6.0, 2.0, 7.0, 1.0, 1.0, 4.0, 3.0, 4.0], [4.0, 9.0, 2.0, 2.0, 4.0, 8.0, 7.0, 7.0, 1.0, 5.0]]
    PB.suite_release = '2.1.0'
    PB.base_function = 'salomon'
    PB.m = 9
    PB.pattern = 'K'
    PB.strategy = 'baseline'
    PB.scale_factors_C = []
    PB.lb_orig = [-100.0, -100.0, -63.779328, -70.41232819658475, -69.91407690731508, -62.55383214873429, -100.0, -100.0, -100.0, -100.0]
    PB.ub_orig = [63.399169921875, 61.59325306101712, 100.0, 100.0, 100.0, 100.0, 69.5715808000383, 64.0719147336193, 70.75630804635489, 61.75997375866169]
    PB.lb_C_work = [-100.0, -100.0, -63.779328, -70.41232819658475, -69.91407690731508, -62.55383214873429, -100.0, -100.0, -100.0, -100.0]
    PB.ub_C_work = [63.399169921875, 61.59325306101712, 100.0, 100.0, 100.0, 100.0, 69.5715808000383, 64.0719147336193, 70.75630804635489, 61.75997375866169]
    PB.ProblemData = _problem_data()
    PB.func_idx = _objective
    PB.func_F = lambda xr: _objective(_cell2idx(xr))
    PB.idx2cell = _idx2cell
    PB.cell2idx = _cell2idx
    PB.decode = lambda xr: _decode(_cell2idx(xr))
    return PB
