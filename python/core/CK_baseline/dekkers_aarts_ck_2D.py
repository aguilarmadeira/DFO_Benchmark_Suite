"""dekkers_aarts_ck_2D  DFO Benchmark Suite v2.1.0 core, Python version (self-contained; standard library only).
DEKKERS_AARTS_CK_2D  DFO Benchmark Suite v2.1.0 core: dekkers_aarts, dimension 2, pattern CK, strategy baseline.

 Copyright (C) 2026 J. F. A. Madeira

 Objective F(z) = f(T(z)), where f is the original function (embedded below, translated from the MATLAB original) and T acts on the finite
 coordinates only: uniform grid of 65 levels on the interval of the domain, piecewise-affine map with |T - Id| <=
 half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
 Global optimum of F on this domain: f* = -24776.51834231769 (A: analytic reduction, 1-D root to machine precision); the optimal solution is listed in
 reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
 Variables: x1:C x2:K


Python version of the MATLAB wrapper dekkers_aarts_ck_2D.m: same instance, same values (level indices 1..m as in MATLAB).
Usage:  from dekkers_aarts_ck_2D import dekkers_aarts_ck_2D;  PB = dekkers_aarts_ck_2D();  f = PB.func_idx(PB.x0_idx)
"""
# Copyright (C) 2026 J. F. A. Madeira. LGPL-3.0-or-later.
import math
from types import SimpleNamespace

_NAME = 'dekkers_aarts_ck_2D'
_N = 2
_TYPES = ['C', 'K']
_CARD = [math.inf, 65.0]
_LB = [-20.0, -12.145811106030584]
_UB = [14.025146484375, 20.0]
_LBW = [-20.0, -12.145811106030584]
_UBW = [14.025146484375, 20.0]
_TB = [math.nan, 14.977217014682724]
_TI = [math.nan, 14.94511215764898]
_CW = [0.0, math.nan]
_CX = [0.0, math.nan]
_M = 65
_PERM = {1: [30, 34, 13, 65, 10, 40, 29, 44, 18, 46, 41, 42, 33, 49, 11, 20, 57, 16, 39, 27, 23, 63, 64, 50, 61, 55, 60, 22, 5, 47, 56, 45, 25, 51, 7, 2, 37, 1, 8, 14, 32, 15, 53, 38, 62, 54, 21, 12, 28, 3, 26, 19, 9, 52, 17, 35, 59, 24, 31, 58, 48, 4, 36, 43, 6]}

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
    return 10 ** 5 * x[0] ** 2 + x[1] ** 2 - (x[0] ** 2 + x[1] ** 2) ** 2 + 1 / 10 ** 5 * (x[0] ** 2 + x[1] ** 2) ** 4

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

def dekkers_aarts_ck_2D():
    PB = SimpleNamespace()
    PB.name = _NAME
    PB.nvar = 2
    PB.var_types = ['C', 'K']
    PB.cardinality = [math.inf, 65.0]
    PB.f_star = -24776.51834231769
    PB.f_star_justification = 'A'
    PB.n_minimizers_documented = 1
    PB.description = 'DFO Benchmark Suite v2.1.0 core (Python): CK adaptation of dekkers_aarts (F = f o T).'
    PB.x0_idx = [-9.40036940574646, 47.0]
    PB.x0_alt_idx = [[-1.3758841753005981, 49.0], [-18.12263596057892, 44.0], [12.579742312431335, 48.0], [-4.167009472846985, 49.0], [4.0070003271102905, 45.0], [-12.73975145816803, 47.0], [8.393054366111755, 51.0], [-8.353697419166565, 45.0], [-0.27937066555023193, 47.0]]
    PB.suite_release = '2.1.0'
    PB.base_function = 'dekkers_aarts'
    PB.m = 65
    PB.pattern = 'CK'
    PB.strategy = 'baseline'
    PB.scale_factors_C = [1.0]
    PB.lb_orig = [-20.0, -12.145811106030584]
    PB.ub_orig = [14.025146484375, 20.0]
    PB.lb_C_work = [-20.0, -12.145811106030584]
    PB.ub_C_work = [14.025146484375, 20.0]
    PB.ProblemData = _problem_data()
    PB.func_idx = _objective
    PB.func_F = lambda xr: _objective(_cell2idx(xr))
    PB.idx2cell = _idx2cell
    PB.cell2idx = _cell2idx
    PB.decode = lambda xr: _decode(_cell2idx(xr))
    return PB
