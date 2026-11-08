"""griewank_cd_5D  DFO Benchmark Suite v2.1.0 core, Python version (self-contained; standard library only).
GRIEWANK_CD_5D  DFO Benchmark Suite v2.1.0 core: griewank, dimension 5, pattern CD, strategy extreme.

 Copyright (C) 2026 J. F. A. Madeira

 Objective F(z) = f(T(z)), where f is the original function (embedded below, translated from the MATLAB original) and T acts on the finite
 coordinates only: uniform grid of 151 levels on the interval of the domain, piecewise-affine map with |T - Id| <=
 half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
 Global optimum of F on this domain: f* = 0 (E: exact by structure); the optimal solution is listed in
 reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
 Variables: x1:C x2:D x3:C x4:D x5:C


Python version of the MATLAB wrapper griewank_cd_5D.m: same instance, same values (level indices 1..m as in MATLAB).
Usage:  from griewank_cd_5D import griewank_cd_5D;  PB = griewank_cd_5D();  f = PB.func_idx(PB.x0_idx)
"""
# Copyright (C) 2026 J. F. A. Madeira. LGPL-3.0-or-later.
import math
from types import SimpleNamespace

_NAME = 'griewank_cd_5D'
_N = 5
_TYPES = ['C', 'D', 'C', 'D', 'C']
_CARD = [math.inf, 151.0, math.inf, 151.0, math.inf]
_LB = [-600.0, -600.0, -379.987968, -600.0, -600.0]
_UB = [363.00439453125, 431.78174058832496, 600.0, 415.6168263223657, 405.6001639232293]
_LBW = [-600.0, -600.0, -37998796800.0, -600.0, -60000000000.0]
_UBW = [363.00439453125, 431.78174058832496, 60000000000.0, 415.6168263223657, 40560016392.32292]
_TB = [math.nan, -1.5665904587714294, math.nan, 2.599316951270339, math.nan]
_TI = [math.nan, 0.0, math.nan, 0.0, math.nan]
_CW = [0.0, math.nan, 0.0, math.nan, 0.0]
_CX = [0.0, math.nan, 0.0, math.nan, 0.0]
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

def griewank_cd_5D():
    PB = SimpleNamespace()
    PB.name = _NAME
    PB.nvar = 5
    PB.var_types = ['C', 'D', 'C', 'D', 'C']
    PB.cardinality = [math.inf, 151.0, math.inf, 151.0, math.inf]
    PB.f_star = 0.0
    PB.f_star_justification = 'E'
    PB.n_minimizers_documented = 1
    PB.description = 'DFO Benchmark Suite v2.1.0 core (Python): CD adaptation of griewank (F = f o T).'
    PB.x0_idx = [362.0639605522156, 48.0, 36982042607.61601, 44.0, -55240209592.249176]
    PB.x0_alt_idx = [[-133.07452940940857, 37.0, -7705408733.183994, 26.0, 22049720362.180832], [115.20004105567932, 41.0, 56738600042.496, 33.0, -9455559003.407936], [-370.0638921260834, 51.0, 12051148701.695992, 39.0, -40960838368.99671], [239.33732628822327, 39.0, -37152087195.648, 25.0, 36329091585.433304], [-251.56921076774597, 43.0, 31055075377.152008, 31.0, 3999606088.632263], [-3.294640302658081, 47.0, -17395529760.767998, 35.0, -18439405833.621613], [-488.5585734844208, 34.0, 50811632812.032, 41.0, -49944685199.21039], [301.40596890449524, 44.0, 2361027674.1119995, 27.0, 27345244755.219635], [-192.32187008857727, 48.0, -23510654681.087997, 34.0, -4160034610.369156]]
    PB.suite_release = '2.1.0'
    PB.base_function = 'griewank'
    PB.m = 151
    PB.pattern = 'CD'
    PB.strategy = 'extreme'
    PB.scale_factors_C = [1.0, 100000000.0, 100000000.0]
    PB.lb_orig = [-600.0, -600.0, -379.987968, -600.0, -600.0]
    PB.ub_orig = [363.00439453125, 431.78174058832496, 600.0, 415.6168263223657, 405.6001639232293]
    PB.lb_C_work = [-600.0, -600.0, -37998796800.0, -600.0, -60000000000.0]
    PB.ub_C_work = [363.00439453125, 431.78174058832496, 60000000000.0, 415.6168263223657, 40560016392.32292]
    PB.ProblemData = _problem_data()
    PB.func_idx = _objective
    PB.func_F = lambda xr: _objective(_cell2idx(xr))
    PB.idx2cell = _idx2cell
    PB.cell2idx = _cell2idx
    PB.decode = lambda xr: _decode(_cell2idx(xr))
    return PB
