"""neumaier3_d_10D  DFO Benchmark Suite v2.1.0 core, Python version (self-contained; standard library only).
NEUMAIER3_D_10D  DFO Benchmark Suite v2.1.0 core: neumaier3, dimension 10, pattern D, strategy baseline.

 Copyright (C) 2026 J. F. A. Madeira

 Objective F(z) = f(T(z)), where f is the original function (embedded below, translated from the MATLAB original) and T acts on the finite
 coordinates only: uniform grid of 9 levels on the interval of the domain, piecewise-affine map with |T - Id| <=
 half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
 Global optimum of F on this domain: f* = -210 (E: exact by structure); the optimal solution is listed in
 reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
 Variables: x1:D x2:D x3:D x4:D x5:D x6:D x7:D x8:D x9:D x10:D


Python version of the MATLAB wrapper neumaier3_d_10D.m: same instance, same values (level indices 1..m as in MATLAB).
Usage:  from neumaier3_d_10D import neumaier3_d_10D;  PB = neumaier3_d_10D();  f = PB.func_idx(PB.x0_idx)
"""
# Copyright (C) 2026 J. F. A. Madeira. LGPL-3.0-or-later.
import math
from types import SimpleNamespace

_NAME = 'neumaier3_d_10D'
_N = 10
_TYPES = ['D', 'D', 'D', 'D', 'D', 'D', 'D', 'D', 'D', 'D']
_CARD = [9.0, 9.0, 9.0, 9.0, 9.0, 9.0, 9.0, 9.0, 9.0, 9.0]
_LB = [-100.0, -100.0, -100.0, -72.69804248229904, -66.21159756847209, -71.1692167641189, -100.0, -100.0, -100.0, -100.0]
_UB = [61.813232421875, 69.88954935731343, 72.739328, 100.0, 100.0, 100.0, 62.20825900073035, 69.96665157572455, 61.65233829210158, 65.62204272417893]
_LBW = [-100.0, -100.0, -100.0, -72.69804248229904, -66.21159756847209, -71.1692167641189, -100.0, -100.0, -100.0, -100.0]
_UBW = [61.813232421875, 69.88954935731343, 72.739328, 100.0, 100.0, 100.0, 62.20825900073035, 69.96665157572455, 61.65233829210158, 65.62204272417893]
_TB = [1.133270263671875, 27.417162017985078, 29.554496, 35.238234069137874, 37.670650911822975, 35.81154371345542, 21.65619425054777, 27.474988681793405, 21.239253719076174, 3.5137767026118354]
_TI = [10.0, 18.0, 24.0, 28.0, 30.0, 30.0, 28.0, 24.0, 18.0, 10.0]
_CW = [math.nan, math.nan, math.nan, math.nan, math.nan, math.nan, math.nan, math.nan, math.nan, math.nan]
_CX = [math.nan, math.nan, math.nan, math.nan, math.nan, math.nan, math.nan, math.nan, math.nan, math.nan]
_M = 9
_PERM = {}

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
    s1 = 0.0
    for v in x: s1 += (v - 1) ** 2
    s2 = 0.0
    for i in range(1, d): s2 += x[i] * x[i - 1]
    return s1 - s2

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

def neumaier3_d_10D():
    PB = SimpleNamespace()
    PB.name = _NAME
    PB.nvar = 10
    PB.var_types = ['D', 'D', 'D', 'D', 'D', 'D', 'D', 'D', 'D', 'D']
    PB.cardinality = [9.0, 9.0, 9.0, 9.0, 9.0, 9.0, 9.0, 9.0, 9.0, 9.0]
    PB.f_star = -210.0
    PB.f_star_justification = 'E'
    PB.n_minimizers_documented = 1
    PB.description = 'DFO Benchmark Suite v2.1.0 core (Python): D adaptation of neumaier3 (F = f o T).'
    PB.x0_idx = [8.0, 2.0, 2.0, 8.0, 2.0, 6.0, 8.0, 5.0, 1.0, 7.0]
    PB.x0_alt_idx = [[4.0, 2.0, 7.0, 8.0, 8.0, 5.0, 6.0, 7.0, 7.0, 8.0], [6.0, 1.0, 3.0, 9.0, 6.0, 4.0, 3.0, 2.0, 5.0, 1.0], [2.0, 1.0, 9.0, 9.0, 3.0, 2.0, 1.0, 5.0, 3.0, 3.0], [8.0, 2.0, 5.0, 8.0, 2.0, 1.0, 7.0, 7.0, 9.0, 4.0], [4.0, 1.0, 1.0, 8.0, 7.0, 8.0, 5.0, 2.0, 7.0, 6.0], [6.0, 1.0, 7.0, 9.0, 5.0, 7.0, 2.0, 4.0, 5.0, 7.0], [2.0, 2.0, 3.0, 8.0, 3.0, 6.0, 8.0, 7.0, 3.0, 8.0], [7.0, 1.0, 8.0, 8.0, 9.0, 4.0, 5.0, 1.0, 9.0, 2.0], [3.0, 2.0, 4.0, 9.0, 6.0, 3.0, 3.0, 4.0, 6.0, 3.0]]
    PB.suite_release = '2.1.0'
    PB.base_function = 'neumaier3'
    PB.m = 9
    PB.pattern = 'D'
    PB.strategy = 'baseline'
    PB.scale_factors_C = []
    PB.lb_orig = [-100.0, -100.0, -100.0, -72.69804248229904, -66.21159756847209, -71.1692167641189, -100.0, -100.0, -100.0, -100.0]
    PB.ub_orig = [61.813232421875, 69.88954935731343, 72.739328, 100.0, 100.0, 100.0, 62.20825900073035, 69.96665157572455, 61.65233829210158, 65.62204272417893]
    PB.lb_C_work = [-100.0, -100.0, -100.0, -72.69804248229904, -66.21159756847209, -71.1692167641189, -100.0, -100.0, -100.0, -100.0]
    PB.ub_C_work = [61.813232421875, 69.88954935731343, 72.739328, 100.0, 100.0, 100.0, 62.20825900073035, 69.96665157572455, 61.65233829210158, 65.62204272417893]
    PB.ProblemData = _problem_data()
    PB.func_idx = _objective
    PB.func_F = lambda xr: _objective(_cell2idx(xr))
    PB.idx2cell = _idx2cell
    PB.cell2idx = _cell2idx
    PB.decode = lambda xr: _decode(_cell2idx(xr))
    return PB
