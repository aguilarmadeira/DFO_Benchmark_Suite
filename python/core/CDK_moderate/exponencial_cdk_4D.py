"""exponencial_cdk_4D  DFO Benchmark Suite v2.1.0 core, Python version (self-contained; standard library only).
EXPONENCIAL_CDK_4D  DFO Benchmark Suite v2.1.0 core: exponencial, dimension 4, pattern CDK, strategy moderate.

 Copyright (C) 2026 J. F. A. Madeira

 Objective F(z) = f(T(z)), where f is the original function (embedded below, translated from the MATLAB original) and T acts on the finite
 coordinates only: uniform grid of 151 levels on the interval of the domain, piecewise-affine map with |T - Id| <=
 half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
 Global optimum of F on this domain: f* = -1 (E: exact by structure); the optimal solution is listed in
 reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
 Variables: x1:C x2:K x3:D x4:C


Python version of the MATLAB wrapper exponencial_cdk_4D.m: same instance, same values (level indices 1..m as in MATLAB).
Usage:  from exponencial_cdk_4D import exponencial_cdk_4D;  PB = exponencial_cdk_4D();  f = PB.func_idx(PB.x0_idx)
"""
# Copyright (C) 2026 J. F. A. Madeira. LGPL-3.0-or-later.
import math
from types import SimpleNamespace

_NAME = 'exponencial_cdk_4D'
_N = 4
_TYPES = ['C', 'K', 'D', 'C']
_CARD = [math.inf, 151.0, 151.0, math.inf]
_LB = [-0.68375732421875, -0.7109942590052329, -0.6165132799999999, -1.0]
_UB = [1.0, 1.0, 1.0, 0.735551853394419]
_LBW = [-0.68375732421875, -0.7109942590052329, -0.6165132799999999, -100.0]
_UBW = [1.0, 1.0, 1.0, 73.55518533944189]
_TB = [math.nan, -0.0037832986164033633, -0.0022382335999999725, math.nan]
_TI = [math.nan, 0.0, 0.0, math.nan]
_CW = [0.0, math.nan, math.nan, 0.0]
_CX = [0.0, math.nan, math.nan, 0.0]
_M = 151
_PERM = {1: [112, 146, 60, 82, 73, 111, 32, 59, 122, 68, 19, 48, 1, 2, 93, 37, 66, 8, 43, 38, 95, 16, 18, 75, 96, 102, 106, 10, 13, 15, 61, 126, 94, 55, 83, 35, 5, 12, 54, 104, 135, 26, 110, 76, 40, 142, 78, 63, 71, 23, 87, 98, 128, 77, 116, 101, 109, 58, 115, 20, 41, 49, 131, 80, 14, 27, 36, 91, 56, 86, 129, 145, 123, 141, 84, 51, 92, 9, 133, 62, 6, 120, 151, 124, 113, 103, 34, 46, 90, 127, 50, 65, 4, 52, 125, 137, 134, 118, 47, 70, 39, 139, 147, 140, 25, 108, 42, 149, 45, 28, 72, 130, 24, 17, 74, 30, 97, 107, 148, 85, 29, 143, 67, 79, 144, 22, 100, 89, 132, 99, 136, 114, 31, 81, 21, 150, 105, 53, 88, 64, 117, 3, 57, 119, 44, 121, 11, 69, 7, 138, 33]}

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
    return -math.exp(-0.5 * s)

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

def exponencial_cdk_4D():
    PB = SimpleNamespace()
    PB.name = _NAME
    PB.nvar = 4
    PB.var_types = ['C', 'K', 'D', 'C']
    PB.cardinality = [math.inf, 151.0, 151.0, math.inf]
    PB.f_star = -1.0
    PB.f_star_justification = 'E'
    PB.n_minimizers_documented = 1
    PB.description = 'DFO Benchmark Suite v2.1.0 core (Python): CDK adaptation of exponencial (F = f o T).'
    PB.x0_idx = [0.05124221086502079, 42.0, 134.0, 66.10987751355162]
    PB.x0_alt_idx = [[0.4483392751216888, 46.0, 65.0, 70.23009640749092], [-0.3803850328922272, 36.0, 20.0, 53.53236720573702], [0.6653861176967621, 45.0, 96.0, 61.122244115625136], [-0.1633381903171539, 48.0, 49.0, 68.71212102551326], [0.24115819811820982, 37.0, 125.0, 52.01439182375938], [-0.5875661098957062, 42.0, 79.0, 56.13461071769865], [0.9317617881298066, 51.0, 11.0, 63.72448762758677], [0.10303748011589053, 38.0, 110.0, 71.31436453747492], [0.5026009857654572, 42.0, 41.0, 54.616635335721014]]
    PB.suite_release = '2.1.0'
    PB.base_function = 'exponencial'
    PB.m = 151
    PB.pattern = 'CDK'
    PB.strategy = 'moderate'
    PB.scale_factors_C = [1.0, 100.0]
    PB.lb_orig = [-0.68375732421875, -0.7109942590052329, -0.6165132799999999, -1.0]
    PB.ub_orig = [1.0, 1.0, 1.0, 0.735551853394419]
    PB.lb_C_work = [-0.68375732421875, -0.7109942590052329, -0.6165132799999999, -100.0]
    PB.ub_C_work = [1.0, 1.0, 1.0, 73.55518533944189]
    PB.ProblemData = _problem_data()
    PB.func_idx = _objective
    PB.func_F = lambda xr: _objective(_cell2idx(xr))
    PB.idx2cell = _idx2cell
    PB.cell2idx = _cell2idx
    PB.decode = lambda xr: _decode(_cell2idx(xr))
    return PB
