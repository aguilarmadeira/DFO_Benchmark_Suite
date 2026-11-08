"""sphere_cdk_3D  DFO Benchmark Suite v2.1.0 core, Python version (self-contained; standard library only).
SPHERE_CDK_3D  DFO Benchmark Suite v2.1.0 core: sphere, dimension 3, pattern CDK, strategy baseline.

 Copyright (C) 2026 J. F. A. Madeira

 Objective F(z) = f(T(z)), where f is the original function (embedded below, translated from the MATLAB original) and T acts on the finite
 coordinates only: uniform grid of 151 levels on the interval of the domain, piecewise-affine map with |T - Id| <=
 half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
 Global optimum of F on this domain: f* = 0 (E: exact by structure); the optimal solution is listed in
 reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
 Variables: x1:C x2:K x3:D


Python version of the MATLAB wrapper sphere_cdk_3D.m: same instance, same values (level indices 1..m as in MATLAB).
Usage:  from sphere_cdk_3D import sphere_cdk_3D;  PB = sphere_cdk_3D();  f = PB.func_idx(PB.x0_idx)
"""
# Copyright (C) 2026 J. F. A. Madeira. LGPL-3.0-or-later.
import math
from types import SimpleNamespace

_NAME = 'sphere_cdk_3D'
_N = 3
_TYPES = ['C', 'K', 'D']
_CARD = [math.inf, 151.0, 151.0]
_LB = [-3.1116375, -5.12, -3.7471911936]
_UB = [5.12, 3.149641497739166, 5.12]
_LBW = [-3.1116375, -5.12, -3.7471911936]
_UBW = [5.12, 3.149641497739166, 5.12]
_TB = [math.nan, 0.0071777285982825845, -0.022970892287999867]
_TI = [math.nan, 0.0, 0.0]
_CW = [0.0, math.nan, math.nan]
_CX = [0.0, math.nan, math.nan]
_M = 151
_PERM = {1: [2, 17, 50, 49, 147, 14, 35, 78, 130, 75, 33, 134, 148, 142, 90, 83, 73, 53, 85, 19, 82, 41, 67, 145, 93, 39, 88, 151, 9, 12, 101, 84, 45, 1, 135, 139, 136, 122, 31, 36, 132, 51, 13, 150, 29, 57, 74, 146, 91, 80, 54, 87, 64, 95, 20, 124, 44, 6, 16, 38, 11, 23, 127, 70, 106, 77, 40, 63, 113, 43, 137, 129, 140, 7, 104, 123, 115, 60, 141, 138, 28, 56, 125, 107, 109, 58, 120, 27, 66, 32, 55, 133, 111, 144, 128, 37, 30, 26, 97, 34, 46, 24, 65, 47, 10, 5, 114, 119, 110, 25, 96, 98, 3, 62, 112, 126, 118, 61, 72, 8, 68, 81, 4, 100, 149, 69, 99, 143, 102, 15, 117, 116, 79, 108, 76, 59, 52, 121, 86, 131, 18, 22, 48, 21, 42, 105, 89, 92, 71, 94, 103]}

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

def sphere_cdk_3D():
    PB = SimpleNamespace()
    PB.name = _NAME
    PB.nvar = 3
    PB.var_types = ['C', 'K', 'D']
    PB.cardinality = [math.inf, 151.0, 151.0]
    PB.f_star = 0.0
    PB.f_star_justification = 'E'
    PB.n_minimizers_documented = 1
    PB.description = 'DFO Benchmark Suite v2.1.0 core (Python): CDK adaptation of sphere (F = f o T).'
    PB.x0_idx = [4.866780682373047, 102.0, 10.0]
    PB.x0_alt_idx = [[0.8152716003417968, 111.0, 115.0], [2.8410261413574216, 115.0, 40.0], [-1.2104829406738282, 103.0, 145.0], [3.8539034118652338, 108.0, 69.0], [-0.19760567016601582, 117.0, 30.0], [1.8281488708496094, 106.0, 100.0], [-2.2233602111816406, 109.0, 60.0], [4.36034204711914, 113.0, 130.0], [0.30883296508789027, 102.0, 90.0]]
    PB.suite_release = '2.1.0'
    PB.base_function = 'sphere'
    PB.m = 151
    PB.pattern = 'CDK'
    PB.strategy = 'baseline'
    PB.scale_factors_C = [1.0]
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
