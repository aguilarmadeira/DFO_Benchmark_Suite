"""shekel_45_cd_4D  DFO Benchmark Suite v2.1.0 core, Python version (self-contained; standard library only).
SHEKEL_45_CD_4D  DFO Benchmark Suite v2.1.0 core: shekel_45, dimension 4, pattern CD, strategy moderate.

 Copyright (C) 2026 J. F. A. Madeira

 Objective F(z) = f(T(z)), where f is the original function (embedded below, translated from the MATLAB original) and T acts on the finite
 coordinates only: uniform grid of 151 levels on the interval of the domain, piecewise-affine map with |T - Id| <=
 half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
 Global optimum of F on this domain: f* = -10.15319967905823 (L: literature value at a refined point); the optimal solution is listed in
 reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
 Variables: x1:C x2:D x3:C x4:D


Python version of the MATLAB wrapper shekel_45_cd_4D.m: same instance, same values (level indices 1..m as in MATLAB).
Usage:  from shekel_45_cd_4D import shekel_45_cd_4D;  PB = shekel_45_cd_4D();  f = PB.func_idx(PB.x0_idx)
"""
# Copyright (C) 2026 J. F. A. Madeira. LGPL-3.0-or-later.
import math
from types import SimpleNamespace

_NAME = 'shekel_45_cd_4D'
_N = 4
_TYPES = ['C', 'D', 'C', 'D']
_CARD = [math.inf, 151.0, math.inf, 151.0]
_LB = [1.9175415039062502, 1.3019560026418735, 1.6990336, 1.8793835901707623]
_UB = [10.0, 10.0, 10.0, 10.0]
_LBW = [1.91754150390625, 1.3019560026418735, 169.90336000000002, 1.8793835901707623]
_UBW = [10.0, 10.0, 1000.0, 10.0]
_TB = [math.nan, 4.027343121814086, math.nan, 3.990743856726364]
_TI = [math.nan, 4.000133276236526, math.nan, 4.00013327718805]
_CW = [4.000037153719125, math.nan, 400.0037153472134, math.nan]
_CX = [4.000037153719125, math.nan, 4.000037153472134, math.nan]
_M = 151
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

_A = [[4, 4, 4, 4, 0.1], [1, 1, 1, 1, 0.2], [8, 8, 8, 8, 0.2], [6, 6, 6, 6, 0.4], [3, 7, 3, 7, 0.4],
      [2, 9, 2, 9, 0.6], [5, 5, 3, 3, 0.3], [8, 1, 8, 1, 0.7], [6, 2, 6, 2, 0.5], [7, 3.6, 7, 3.6, 0.5]]
def _orig(x):
    s = 0.0
    for r in _A[:5]:
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

def shekel_45_cd_4D():
    PB = SimpleNamespace()
    PB.name = _NAME
    PB.nvar = 4
    PB.var_types = ['C', 'D', 'C', 'D']
    PB.cardinality = [math.inf, 151.0, math.inf, 151.0]
    PB.f_star = -10.153199679058229
    PB.f_star_justification = 'L'
    PB.n_minimizers_documented = 1
    PB.description = 'DFO Benchmark Suite v2.1.0 core (Python): CD adaptation of shekel_45 (F = f o T).'
    PB.x0_idx = [9.246216028928757, 33.0, 645.6483463168001, 93.0]
    PB.x0_alt_idx = [[5.268130987882614, 22.0, 267.1242784768, 100.0], [7.257173508405685, 27.0, 812.9958289408, 106.0], [3.279088467359543, 29.0, 434.47176110080005, 92.0], [8.251694768667221, 19.0, 973.9681693696, 98.0], [4.2736097276210785, 28.0, 595.4441015296, 105.0], [6.26265224814415, 32.0, 185.04432271360002, 87.0], [2.284567207098007, 20.0, 762.7915841536001, 94.0], [9.506685882806778, 24.0, 352.39180533760003, 101.0], [5.528600841760635, 33.0, 930.1390667776, 107.0]]
    PB.suite_release = '2.1.0'
    PB.base_function = 'shekel_45'
    PB.m = 151
    PB.pattern = 'CD'
    PB.strategy = 'moderate'
    PB.scale_factors_C = [1.0, 100.0]
    PB.lb_orig = [1.9175415039062502, 1.3019560026418735, 1.6990336, 1.8793835901707623]
    PB.ub_orig = [10.0, 10.0, 10.0, 10.0]
    PB.lb_C_work = [1.91754150390625, 1.3019560026418735, 169.90336000000002, 1.8793835901707623]
    PB.ub_C_work = [10.0, 10.0, 1000.0, 10.0]
    PB.ProblemData = _problem_data()
    PB.func_idx = _objective
    PB.func_F = lambda xr: _objective(_cell2idx(xr))
    PB.idx2cell = _idx2cell
    PB.cell2idx = _cell2idx
    PB.decode = lambda xr: _decode(_cell2idx(xr))
    return PB
