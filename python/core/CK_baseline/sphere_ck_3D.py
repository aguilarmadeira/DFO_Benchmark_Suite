"""sphere_ck_3D  DFO Benchmark Suite v2.1.0 core, Python version (self-contained; standard library only).
SPHERE_CK_3D  DFO Benchmark Suite v2.1.0 core: sphere, dimension 3, pattern CK, strategy baseline.

 Copyright (C) 2026 J. F. A. Madeira

 Objective F(z) = f(T(z)), where f is the original function (embedded below, translated from the MATLAB original) and T acts on the finite
 coordinates only: uniform grid of 65 levels on the interval of the domain, piecewise-affine map with |T - Id| <=
 half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
 Global optimum of F on this domain: f* = 0 (E: exact by structure); the optimal solution is listed in
 reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
 Variables: x1:C x2:K x3:C


Python version of the MATLAB wrapper sphere_ck_3D.m: same instance, same values (level indices 1..m as in MATLAB).
Usage:  from sphere_ck_3D import sphere_ck_3D;  PB = sphere_ck_3D();  f = PB.func_idx(PB.x0_idx)
"""
# Copyright (C) 2026 J. F. A. Madeira. LGPL-3.0-or-later.
import math
from types import SimpleNamespace

_NAME = 'sphere_ck_3D'
_N = 3
_TYPES = ['C', 'K', 'C']
_CARD = [math.inf, 65.0, math.inf]
_LB = [-3.1116375, -5.12, -3.7471911936]
_UB = [5.12, 3.149641497739166, 5.12]
_LBW = [-3.1116375, -5.12, -3.7471911936]
_UBW = [5.12, 3.149641497739166, 5.12]
_TB = [math.nan, 0.048525936086978305, math.nan]
_TI = [math.nan, 0.0, math.nan]
_CW = [0.0, math.nan, 0.0]
_CX = [0.0, math.nan, 0.0]
_M = 65
_PERM = {1: [10, 53, 63, 8, 25, 26, 60, 57, 14, 65, 12, 47, 21, 44, 55, 4, 11, 56, 2, 27, 13, 15, 36, 6, 48, 23, 61, 38, 29, 18, 42, 59, 35, 3, 54, 24, 1, 31, 64, 50, 32, 16, 34, 5, 52, 41, 49, 62, 28, 22, 20, 46, 33, 51, 7, 58, 19, 9, 17, 43, 37, 39, 30, 40, 45]}

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
    return s

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

def sphere_ck_3D():
    PB = SimpleNamespace()
    PB.name = _NAME
    PB.nvar = 3
    PB.var_types = ['C', 'K', 'C']
    PB.cardinality = [math.inf, 65.0, math.inf]
    PB.f_star = 0.0
    PB.f_star_justification = 'E'
    PB.n_minimizers_documented = 1
    PB.description = 'DFO Benchmark Suite v2.1.0 core (Python): CK adaptation of sphere (F = f o T).'
    PB.x0_idx = [4.866780682373047, 44.0, -3.205228467847168]
    PB.x0_alt_idx = [[0.8152716003417968, 48.0, 2.9663366028984335], [2.8410261413574216, 50.0, -1.4176027232174078], [-1.2104829406738282, 45.0, 4.753962347528192], [3.8539034118652338, 47.0, 0.30192299304550385], [-0.19760567016601582, 50.0, -2.0390154820648956], [1.8281488708496094, 46.0, 2.0895487376752655], [-2.2233602111816406, 47.0, -0.25138973743513615], [4.36034204711914, 49.0, 3.8771744823050245], [0.30883296508789027, 44.0, 1.5362360071946242]]
    PB.suite_release = '2.1.0'
    PB.base_function = 'sphere'
    PB.m = 65
    PB.pattern = 'CK'
    PB.strategy = 'baseline'
    PB.scale_factors_C = [1.0, 1.0]
    PB.lb_orig = [-3.1116375, -5.12, -3.7471911936]
    PB.ub_orig = [5.12, 3.149641497739166, 5.12]
    PB.lb_C_work = [-3.1116375, -5.12, -3.7471911936]
    PB.ub_C_work = [5.12, 3.149641497739166, 5.12]
    PB.ProblemData = _problem_data()
    PB.func_idx = _objective
    PB.func_F = lambda xr: _objective(_cell2idx(xr))
    PB.idx2cell = _idx2cell
    PB.cell2idx = _cell2idx
    PB.decode = lambda xr: _decode(_cell2idx(xr))
    return PB
