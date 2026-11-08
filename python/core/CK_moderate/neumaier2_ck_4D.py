"""neumaier2_ck_4D  DFO Benchmark Suite v2.1.0 core, Python version (self-contained; standard library only).
NEUMAIER2_CK_4D  DFO Benchmark Suite v2.1.0 core: neumaier2, dimension 4, pattern CK, strategy moderate.

 Copyright (C) 2026 J. F. A. Madeira

 Objective F(z) = f(T(z)), where f is the original function (embedded below, translated from the MATLAB original) and T acts on the finite
 coordinates only: uniform grid of 151 levels on the interval of the domain, piecewise-affine map with |T - Id| <=
 half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
 Global optimum of F on this domain: f* = 0 (E: exact by structure); the optimal solution is listed in
 reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
 Variables: x1:C x2:K x3:C x4:K


Python version of the MATLAB wrapper neumaier2_ck_4D.m: same instance, same values (level indices 1..m as in MATLAB).
Usage:  from neumaier2_ck_4D import neumaier2_ck_4D;  PB = neumaier2_ck_4D();  f = PB.func_idx(PB.x0_idx)
"""
# Copyright (C) 2026 J. F. A. Madeira. LGPL-3.0-or-later.
import math
from types import SimpleNamespace

_NAME = 'neumaier2_ck_4D'
_N = 4
_TYPES = ['C', 'K', 'C', 'K']
_CARD = [math.inf, 151.0, math.inf, 151.0]
_LB = [0.0, 0.0, 0.75801344, 0.0]
_UB = [3.3762646484375, 3.242235431590713, 4.0, 3.2196751353602666]
_LBW = [0.0, 0.0, 75.80134400000003, 0.0]
_UBW = [3.3762646484375, 3.242235431590713, 400.0, 3.2196751353602666]
_TB = [math.nan, 2.010185967586242, math.nan, 3.0050301263362487]
_TI = [math.nan, 2.0, math.nan, 3.0]
_CW = [1.0, math.nan, 200.0, math.nan]
_CX = [1.0, math.nan, 2.0, math.nan]
_M = 151
_PERM = {1: [91, 124, 59, 110, 103, 13, 106, 38, 151, 65, 3, 105, 131, 41, 85, 140, 16, 18, 96, 121, 28, 45, 118, 119, 56, 40, 44, 114, 9, 87, 43, 83, 148, 48, 104, 123, 7, 17, 19, 60, 29, 90, 70, 2, 8, 94, 149, 75, 92, 81, 112, 84, 126, 79, 20, 31, 39, 73, 5, 23, 50, 77, 113, 14, 35, 117, 147, 80, 21, 143, 61, 25, 76, 150, 30, 55, 74, 24, 53, 98, 47, 138, 127, 141, 108, 22, 51, 145, 10, 36, 97, 57, 135, 66, 89, 109, 37, 15, 52, 82, 72, 34, 68, 33, 46, 86, 139, 116, 120, 107, 134, 11, 93, 137, 115, 12, 125, 102, 88, 128, 63, 78, 130, 95, 99, 49, 100, 62, 32, 42, 144, 133, 26, 54, 146, 136, 6, 69, 142, 122, 58, 132, 111, 129, 67, 4, 64, 101, 27, 71, 1], 3: [99, 89, 90, 71, 23, 86, 88, 85, 47, 27, 56, 83, 39, 48, 30, 14, 124, 36, 100, 5, 63, 42, 93, 15, 101, 146, 151, 6, 97, 12, 20, 46, 26, 50, 35, 135, 144, 111, 114, 126, 31, 131, 40, 118, 145, 81, 72, 32, 87, 49, 13, 129, 51, 123, 22, 33, 41, 96, 8, 69, 53, 95, 7, 141, 119, 55, 128, 80, 113, 79, 140, 103, 59, 77, 34, 60, 1, 117, 58, 19, 52, 75, 4, 92, 18, 28, 143, 73, 115, 104, 122, 64, 11, 109, 127, 133, 45, 24, 61, 91, 120, 44, 130, 43, 136, 142, 38, 112, 74, 10, 139, 25, 108, 21, 16, 29, 116, 3, 107, 150, 84, 134, 132, 121, 125, 76, 138, 94, 147, 78, 37, 67, 137, 102, 54, 57, 70, 2, 68, 65, 17, 98, 105, 9, 148, 106, 66, 110, 82, 149, 62]}

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

_b = [8, 18, 44, 114]
def _orig(x):
    s = 0.0
    for k in range(4):
        t = 0.0
        for v in x: t += v ** (k + 1)
        s += (_b[k] - t) ** 2
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

def neumaier2_ck_4D():
    PB = SimpleNamespace()
    PB.name = _NAME
    PB.nvar = 4
    PB.var_types = ['C', 'K', 'C', 'K']
    PB.cardinality = [math.inf, 151.0, math.inf, 151.0]
    PB.f_star = 0.0
    PB.f_star_justification = 'E'
    PB.n_minimizers_documented = 1
    PB.description = 'DFO Benchmark Suite v2.1.0 core (Python): CK adaptation of neumaier2 (F = f o T).'
    PB.x0_idx = [1.0567312693595887, 138.0, 352.38170140672, 108.0]
    PB.x0_alt_idx = [[1.9271744990348816, 143.0, 266.79325622272006, 115.0], [0.22585364103317263, 147.0, 106.50944069632003, 121.0], [2.9954457354545596, 135.0, 332.15170527232004, 128.0], [1.274342076778412, 144.0, 169.37804406784002, 114.0], [2.144785306453705, 148.0, 395.02030864384, 120.0], [0.44346444845199584, 138.0, 234.73649311744, 124.0], [2.5800069212913512, 140.0, 86.90190598144002, 109.0], [0.8490118622779846, 149.0, 300.09494216704, 116.0], [1.7194550919532776, 139.0, 152.26035503104004, 122.0]]
    PB.suite_release = '2.1.0'
    PB.base_function = 'neumaier2'
    PB.m = 151
    PB.pattern = 'CK'
    PB.strategy = 'moderate'
    PB.scale_factors_C = [1.0, 100.0]
    PB.lb_orig = [0.0, 0.0, 0.75801344, 0.0]
    PB.ub_orig = [3.3762646484375, 3.242235431590713, 4.0, 3.2196751353602666]
    PB.lb_C_work = [0.0, 0.0, 75.80134400000003, 0.0]
    PB.ub_C_work = [3.3762646484375, 3.242235431590713, 400.0, 3.2196751353602666]
    PB.ProblemData = _problem_data()
    PB.func_idx = _objective
    PB.func_F = lambda xr: _objective(_cell2idx(xr))
    PB.idx2cell = _idx2cell
    PB.cell2idx = _cell2idx
    PB.decode = lambda xr: _decode(_cell2idx(xr))
    return PB
