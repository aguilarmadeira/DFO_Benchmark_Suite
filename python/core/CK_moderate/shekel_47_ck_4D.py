"""shekel_47_ck_4D  DFO Benchmark Suite v2.1.0 core, Python version (self-contained; standard library only).
SHEKEL_47_CK_4D  DFO Benchmark Suite v2.1.0 core: shekel_47, dimension 4, pattern CK, strategy moderate.

 Copyright (C) 2026 J. F. A. Madeira

 Objective F(z) = f(T(z)), where f is the original function (embedded below, translated from the MATLAB original) and T acts on the finite
 coordinates only: uniform grid of 151 levels on the interval of the domain, piecewise-affine map with |T - Id| <=
 half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
 Global optimum of F on this domain: f* = -10.40294056681866 (L: literature value at a refined point); the optimal solution is listed in
 reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
 Variables: x1:C x2:K x3:C x4:K


Python version of the MATLAB wrapper shekel_47_ck_4D.m: same instance, same values (level indices 1..m as in MATLAB).
Usage:  from shekel_47_ck_4D import shekel_47_ck_4D;  PB = shekel_47_ck_4D();  f = PB.func_idx(PB.x0_idx)
"""
# Copyright (C) 2026 J. F. A. Madeira. LGPL-3.0-or-later.
import math
from types import SimpleNamespace

_NAME = 'shekel_47_ck_4D'
_N = 4
_TYPES = ['C', 'K', 'C', 'K']
_CARD = [math.inf, 151.0, math.inf, 151.0]
_LB = [1.3487915039062501, 1.5352893359752071, 1.8390336, 1.9793835901707624]
_UB = [10.0, 10.0, 10.0, 10.0]
_LBW = [1.3487915039062504, 1.5352893359752071, 183.90336000000002, 1.9793835901707624]
_UBW = [10.0, 10.0, 1000.0, 10.0]
_TB = [math.nan, 4.018271130755813, math.nan, 4.011273080660835]
_TI = [math.nan, 4.000689366380763, math.nan, 3.999606160106975]
_CW = [4.000572916244606, math.nan, 399.94897114568784, math.nan]
_CX = [4.000572916244606, math.nan, 3.999489711456878, math.nan]
_M = 151
_PERM = {1: [137, 121, 31, 29, 100, 86, 62, 6, 143, 18, 148, 117, 35, 40, 124, 49, 65, 44, 38, 126, 134, 87, 111, 129, 17, 76, 72, 58, 9, 97, 20, 59, 78, 33, 2, 94, 82, 25, 53, 54, 39, 5, 146, 15, 22, 51, 136, 123, 73, 91, 36, 89, 102, 151, 1, 107, 23, 125, 110, 13, 45, 115, 63, 116, 37, 56, 4, 24, 32, 70, 135, 145, 139, 141, 103, 85, 104, 83, 30, 128, 34, 147, 109, 61, 21, 46, 81, 138, 47, 3, 19, 41, 88, 130, 149, 119, 64, 101, 113, 50, 55, 27, 66, 42, 92, 68, 14, 16, 74, 52, 8, 120, 108, 132, 118, 12, 112, 75, 131, 142, 150, 140, 98, 57, 122, 105, 99, 96, 79, 67, 77, 43, 144, 90, 7, 84, 127, 11, 10, 80, 28, 106, 114, 60, 48, 71, 133, 95, 93, 26, 69], 3: [58, 88, 78, 67, 20, 62, 102, 77, 148, 19, 94, 32, 31, 36, 43, 65, 30, 98, 105, 27, 69, 68, 112, 142, 18, 111, 124, 42, 90, 131, 50, 45, 138, 150, 61, 149, 97, 135, 53, 74, 41, 95, 9, 56, 140, 51, 121, 63, 12, 21, 38, 49, 143, 34, 134, 82, 119, 44, 55, 16, 92, 100, 70, 15, 40, 60, 120, 116, 106, 3, 11, 109, 14, 33, 71, 99, 46, 141, 35, 118, 39, 132, 80, 110, 115, 147, 87, 17, 137, 129, 107, 48, 2, 75, 6, 85, 73, 101, 25, 114, 64, 37, 76, 52, 130, 79, 26, 29, 86, 66, 22, 8, 10, 59, 145, 146, 13, 93, 96, 126, 47, 133, 123, 144, 23, 7, 1, 127, 113, 103, 117, 139, 104, 136, 57, 4, 54, 81, 84, 151, 128, 72, 108, 83, 91, 28, 122, 125, 5, 89, 24]}

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
    for r in _A[:7]:
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

def shekel_47_ck_4D():
    PB = SimpleNamespace()
    PB.name = _NAME
    PB.nvar = 4
    PB.var_types = ['C', 'K', 'C', 'K']
    PB.cardinality = [math.inf, 151.0, math.inf, 151.0]
    PB.f_star = -10.402940566818662
    PB.f_star_justification = 'L'
    PB.n_minimizers_documented = 1
    PB.description = 'DFO Benchmark Suite v2.1.0 core (Python): CK adaptation of shekel_47 (F = f o T).'
    PB.x0_idx = [2.164066523313523, 83.0, 814.8439943168001, 114.0]
    PB.x0_alt_idx = [[9.260760992765427, 72.0, 442.70392647680006, 121.0], [4.901362961530685, 77.0, 979.3690769408, 128.0], [7.1317526519298555, 79.0, 607.2290091008, 113.0], [2.721663945913315, 69.0, 197.48320808960003, 120.0], [8.19625682234764, 78.0, 765.4864695296001, 126.0], [3.8368587911129, 82.0, 362.0082907136, 109.0], [6.0672484815120695, 70.0, 930.0115521536001, 115.0], [1.6318144381046298, 74.0, 526.5333733376, 122.0], [9.539559704065322, 83.0, 311.0838603776, 129.0]]
    PB.suite_release = '2.1.0'
    PB.base_function = 'shekel_47'
    PB.m = 151
    PB.pattern = 'CK'
    PB.strategy = 'moderate'
    PB.scale_factors_C = [1.0, 100.0]
    PB.lb_orig = [1.3487915039062501, 1.5352893359752071, 1.8390336, 1.9793835901707624]
    PB.ub_orig = [10.0, 10.0, 10.0, 10.0]
    PB.lb_C_work = [1.3487915039062504, 1.5352893359752071, 183.90336000000002, 1.9793835901707624]
    PB.ub_C_work = [10.0, 10.0, 1000.0, 10.0]
    PB.ProblemData = _problem_data()
    PB.func_idx = _objective
    PB.func_F = lambda xr: _objective(_cell2idx(xr))
    PB.idx2cell = _idx2cell
    PB.cell2idx = _cell2idx
    PB.decode = lambda xr: _decode(_cell2idx(xr))
    return PB
