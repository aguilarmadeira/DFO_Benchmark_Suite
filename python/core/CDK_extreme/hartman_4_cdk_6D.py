"""hartman_4_cdk_6D  DFO Benchmark Suite v2.1.0 core, Python version (self-contained; standard library only).
HARTMAN_4_CDK_6D  DFO Benchmark Suite v2.1.0 core: hartman_4, dimension 6, pattern CDK, strategy extreme.

 Copyright (C) 2026 J. F. A. Madeira

 Objective F(z) = f(T(z)), where f is the original function (embedded below, translated from the MATLAB original) and T acts on the finite
 coordinates only: uniform grid of 25 levels on the interval of the domain, piecewise-affine map with |T - Id| <=
 half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
 Global optimum of F on this domain: f* = -3.322368011391339 (L: literature value at a refined point); the optimal solution is listed in
 reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
 Variables: x1:C x2:K x3:D x4:C x5:D x6:K


Python version of the MATLAB wrapper hartman_4_cdk_6D.m: same instance, same values (level indices 1..m as in MATLAB).
Usage:  from hartman_4_cdk_6D import hartman_4_cdk_6D;  PB = hartman_4_cdk_6D();  f = PB.func_idx(PB.x0_idx)
"""
# Copyright (C) 2026 J. F. A. Madeira. LGPL-3.0-or-later.
import math
from types import SimpleNamespace

_NAME = 'hartman_4_cdk_6D'
_N = 6
_TYPES = ['C', 'K', 'D', 'C', 'D', 'K']
_CARD = [math.inf, 25.0, 25.0, math.inf, 25.0, 25.0]
_LB = [0.0, 0.14795966062084032, 0.14470336, 0.0, 0.18109077248821803, 0.16527817653443508]
_UB = [0.848441162109375, 1.0, 1.0, 0.816347355268638, 1.0, 1.0]
_LBW = [0.0, 0.14795966062084032, 0.14470336, 0.0, 0.18109077248821803, 0.16527817653443508]
_UBW = [0.848441162109375, 1.0, 1.0, 81634735.5268638, 1.0, 1.0]
_TB = [math.nan, 0.14795966062084032, 0.46543959999999995, math.nan, 0.31757564374018166, 0.6521992402226813]
_TI = [math.nan, 0.150011, 0.476874, math.nan, 0.311652, 0.6573]
_CW = [0.20169, math.nan, math.nan, 27533200.000000004, math.nan, math.nan]
_CX = [0.20169, math.nan, math.nan, 0.275332, math.nan, math.nan]
_M = 25
_PERM = {1: [14, 10, 16, 7, 9, 21, 5, 11, 20, 4, 12, 23, 24, 13, 6, 25, 1, 2, 15, 8, 18, 19, 22, 17, 3], 5: [22, 21, 18, 1, 5, 24, 10, 25, 23, 15, 12, 6, 14, 2, 13, 11, 20, 7, 8, 16, 19, 3, 17, 9, 4]}

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

_A3 = [[3, 10, 30, 1, 0.3689, 0.1170, 0.2673],
       [0.1, 10, 35, 1.2, 0.4699, 0.4387, 0.7470],
       [3, 10, 30, 3, 0.1091, 0.8732, 0.5547],
       [0.1, 10, 35, 3.2, 0.03815, 0.5743, 0.8828]]
_A6 = [[10, 3, 17, 3.5, 1.7, 8, 1, 0.1312, 0.1696, 0.5569, 0.0124, 0.8283, 0.5886],
       [0.05, 10, 17, 0.1, 8, 14, 1.2, 0.2329, 0.4135, 0.8307, 0.3736, 0.1004, 0.9991],
       [3, 3.5, 1.7, 10, 17, 8, 3, 0.2348, 0.1451, 0.3522, 0.2883, 0.3047, 0.6650],
       [17, 8, 0.05, 10, 0.1, 14, 3.2, 0.4047, 0.8828, 0.8732, 0.5743, 0.1091, 0.0381]]
def _orig(x):
    n = len(x)
    A = _A3 if n == 3 else _A6
    f = 0.0
    for r in A:
        e = 0.0
        for j in range(n): e += r[j] * (x[j] - r[n + 1 + j]) ** 2
        f += r[n] * math.exp(-e)
    return -f

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

def hartman_4_cdk_6D():
    PB = SimpleNamespace()
    PB.name = _NAME
    PB.nvar = 6
    PB.var_types = ['C', 'K', 'D', 'C', 'D', 'K']
    PB.cardinality = [math.inf, 25.0, 25.0, math.inf, 25.0, 25.0]
    PB.f_star = -3.322368011391339
    PB.f_star_justification = 'L'
    PB.n_minimizers_documented = 1
    PB.description = 'DFO Benchmark Suite v2.1.0 core (Python): CDK adaptation of hartman_4 (F = f o T).'
    PB.x0_idx = [0.21252456843853, 2.0, 6.0, 58140523.84462187, 9.0, 23.0]
    PB.x0_alt_idx = [[0.4312633055448532, 2.0, 23.0, 48654438.37523619, 1.0, 20.0], [0.0037285012006759645, 3.0, 11.0, 52224470.54113402, 20.0, 17.0], [0.838912770152092, 1.0, 5.0, 55794502.70703185, 12.0, 12.0], [0.421320635676384, 3.0, 16.0, 47940431.942056626, 7.0, 9.0], [0.6251453679800034, 3.0, 9.0, 51510464.107954465, 23.0, 6.0], [0.20755323350429536, 2.0, 20.0, 53448481.56944185, 17.0, 2.0], [0.734514736533165, 2.0, 14.0, 57018513.73533968, 10.0, 23.0], [0.316922602057457, 3.0, 2.0, 49164442.970364444, 5.0, 18.0], [0.518261666893959, 2.0, 19.0, 52734475.136262275, 21.0, 15.0]]
    PB.suite_release = '2.1.0'
    PB.base_function = 'hartman_4'
    PB.m = 25
    PB.pattern = 'CDK'
    PB.strategy = 'extreme'
    PB.scale_factors_C = [1.0, 100000000.0]
    PB.lb_orig = [0.0, 0.14795966062084032, 0.14470336, 0.0, 0.18109077248821803, 0.16527817653443508]
    PB.ub_orig = [0.848441162109375, 1.0, 1.0, 0.816347355268638, 1.0, 1.0]
    PB.lb_C_work = [0.0, 0.14795966062084032, 0.14470336, 0.0, 0.18109077248821803, 0.16527817653443508]
    PB.ub_C_work = [0.848441162109375, 1.0, 1.0, 81634735.5268638, 1.0, 1.0]
    PB.ProblemData = _problem_data()
    PB.func_idx = _objective
    PB.func_F = lambda xr: _objective(_cell2idx(xr))
    PB.idx2cell = _idx2cell
    PB.cell2idx = _cell2idx
    PB.decode = lambda xr: _decode(_cell2idx(xr))
    return PB
