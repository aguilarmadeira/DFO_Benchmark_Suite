"""cosine_mixture_cdk_4D  DFO Benchmark Suite v2.1.0 core, Python version (self-contained; standard library only).
COSINE_MIXTURE_CDK_4D  DFO Benchmark Suite v2.1.0 core: cosine_mixture, dimension 4, pattern CDK, strategy baseline.

 Copyright (C) 2026 J. F. A. Madeira

 Objective F(z) = f(T(z)), where f is the original function (embedded below, translated from the MATLAB original) and T acts on the finite
 coordinates only: uniform grid of 151 levels on the interval of the domain, piecewise-affine map with |T - Id| <=
 half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
 Global optimum of F on this domain: f* = -4.4 (E: exact by structure); the optimal solution is listed in
 reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
 Variables: x1:C x2:K x3:D x4:C


Python version of the MATLAB wrapper cosine_mixture_cdk_4D.m: same instance, same values (level indices 1..m as in MATLAB).
Usage:  from cosine_mixture_cdk_4D import cosine_mixture_cdk_4D;  PB = cosine_mixture_cdk_4D();  f = PB.func_idx(PB.x0_idx)
"""
# Copyright (C) 2026 J. F. A. Madeira. LGPL-3.0-or-later.
import math
from types import SimpleNamespace

_NAME = 'cosine_mixture_cdk_4D'
_N = 4
_TYPES = ['C', 'K', 'D', 'C']
_CARD = [math.inf, 151.0, 151.0, math.inf]
_LB = [-0.66625732421875, -0.653957221968196, -0.62211328, -0.6984089962515618]
_UB = [1.0, 1.0, 1.0, 1.0]
_LBW = [-0.66625732421875, -0.653957221968196, -0.62211328, -0.6984089962515618]
_UBW = [1.0, 1.0, 1.0, 1.0]
_TB = [math.nan, 1.0, 1.0, math.nan]
_TI = [math.nan, 1.0, 1.0, math.nan]
_CW = [1.0, math.nan, math.nan, 1.0]
_CX = [1.0, math.nan, math.nan, 1.0]
_M = 151
_PERM = {1: [38, 121, 109, 75, 118, 138, 137, 65, 148, 147, 56, 44, 35, 108, 62, 26, 77, 93, 17, 60, 55, 33, 3, 94, 36, 90, 53, 98, 66, 23, 103, 101, 115, 119, 6, 8, 124, 145, 122, 149, 82, 80, 37, 68, 2, 150, 16, 54, 141, 72, 96, 81, 15, 89, 128, 144, 88, 70, 7, 19, 39, 105, 50, 58, 112, 21, 29, 84, 135, 120, 79, 92, 86, 146, 12, 125, 43, 110, 116, 113, 57, 74, 130, 42, 69, 52, 73, 85, 18, 123, 102, 97, 27, 71, 134, 142, 91, 48, 139, 106, 31, 76, 127, 46, 107, 117, 131, 59, 95, 63, 104, 9, 34, 20, 64, 129, 67, 49, 4, 1, 51, 22, 11, 13, 30, 140, 100, 136, 151, 83, 99, 5, 132, 87, 126, 32, 61, 14, 111, 47, 24, 40, 45, 41, 114, 78, 143, 28, 25, 10, 133]}

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
    s1 = 0.0
    for v in x: s1 += math.cos(5 * math.pi * v)
    s2 = 0.0
    for v in x: s2 += v ** 2
    return 0.1 * s1 - s2

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

def cosine_mixture_cdk_4D():
    PB = SimpleNamespace()
    PB.name = _NAME
    PB.nvar = 4
    PB.var_types = ['C', 'K', 'D', 'C']
    PB.cardinality = [math.inf, 151.0, 151.0, math.inf]
    PB.f_star = -4.4
    PB.f_star_justification = 'E'
    PB.n_minimizers_documented = 1
    PB.description = 'DFO Benchmark Suite v2.1.0 core (Python): CDK adaptation of cosine_mixture (F = f o T).'
    PB.x0_idx = [0.26938521623611456, 59.0, 128.0, -0.28600897758656374]
    PB.x0_alt_idx = [[-0.5873379099369048, 63.0, 88.0, -0.24568856409787615], [0.9357254254817962, 53.0, 14.0, -0.4090923450783471], [0.09608794569969181, 62.0, 118.0, -0.334817899178133], [0.5256699120998383, 64.0, 43.0, -0.2605434532779189], [-0.3237307941913604, 54.0, 147.0, -0.42394723425838987], [0.7306976687908173, 58.0, 73.0, -0.3836268207697022], [-0.10893981099128724, 68.0, 5.0, -0.3093523748694882], [0.3206421554088593, 55.0, 104.0, -0.23507792896927415], [-0.5336401641368865, 59.0, 35.0, -0.398481709949745]]
    PB.suite_release = '2.1.0'
    PB.base_function = 'cosine_mixture'
    PB.m = 151
    PB.pattern = 'CDK'
    PB.strategy = 'baseline'
    PB.scale_factors_C = [1.0, 1.0]
    PB.lb_orig = [-0.66625732421875, -0.653957221968196, -0.62211328, -0.6984089962515618]
    PB.ub_orig = [1.0, 1.0, 1.0, 1.0]
    PB.lb_C_work = [-0.66625732421875, -0.653957221968196, -0.62211328, -0.6984089962515618]
    PB.ub_C_work = [1.0, 1.0, 1.0, 1.0]
    PB.ProblemData = _problem_data()
    PB.func_idx = _objective
    PB.func_F = lambda xr: _objective(_cell2idx(xr))
    PB.idx2cell = _idx2cell
    PB.cell2idx = _cell2idx
    PB.decode = lambda xr: _decode(_cell2idx(xr))
    return PB
