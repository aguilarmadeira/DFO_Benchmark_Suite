"""miele_cantrell_ck_4D  DFO Benchmark Suite v2.1.0 core, Python version (self-contained; standard library only).
MIELE_CANTRELL_CK_4D  DFO Benchmark Suite v2.1.0 core: miele_cantrell, dimension 4, pattern CK, strategy moderate.

 Copyright (C) 2026 J. F. A. Madeira

 Objective F(z) = f(T(z)), where f is the original function (embedded below, translated from the MATLAB original) and T acts on the finite
 coordinates only: uniform grid of 151 levels on the interval of the domain, piecewise-affine map with |T - Id| <=
 half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
 Global optimum of F on this domain: f* = 0 (E: exact by structure); the optimal solution is listed in
 reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
 Variables: x1:C x2:K x3:C x4:K


Python version of the MATLAB wrapper miele_cantrell_ck_4D.m: same instance, same values (level indices 1..m as in MATLAB).
Usage:  from miele_cantrell_ck_4D import miele_cantrell_ck_4D;  PB = miele_cantrell_ck_4D();  f = PB.func_idx(PB.x0_idx)
"""
# Copyright (C) 2026 J. F. A. Madeira. LGPL-3.0-or-later.
import math
from types import SimpleNamespace

_NAME = 'miele_cantrell_ck_4D'
_N = 4
_TYPES = ['C', 'K', 'C', 'K']
_CARD = [math.inf, 151.0, math.inf, 151.0]
_LB = [-10.0, -10.0, -10.0, -10.0]
_UB = [7.2313232421875, 7.144510491286898, 6.769932799999999, 6.4983756768013325]
_LBW = [-10.0, -10.0, -1000.0, -10.0]
_UBW = [7.2313232421875, 7.144510491286898, 676.9932799999999, 6.4983756768013325]
_TB = [math.nan, 0.9724867144236153, math.nan, 0.9989171178675544]
_TI = [math.nan, 1.0, math.nan, 1.0]
_CW = [0.0, math.nan, 100.0, math.nan]
_CX = [0.0, math.nan, 1.0, math.nan]
_M = 151
_PERM = {1: [77, 122, 104, 101, 150, 18, 6, 128, 114, 56, 57, 68, 102, 78, 82, 129, 36, 48, 71, 19, 107, 52, 111, 8, 79, 123, 141, 26, 127, 3, 96, 110, 15, 73, 66, 131, 86, 121, 140, 92, 16, 49, 97, 31, 119, 4, 11, 14, 70, 139, 44, 81, 134, 142, 105, 9, 76, 94, 41, 58, 25, 137, 133, 80, 7, 23, 74, 53, 126, 22, 28, 63, 116, 148, 90, 17, 34, 62, 115, 51, 145, 24, 35, 10, 47, 124, 84, 108, 42, 67, 65, 85, 32, 93, 20, 87, 132, 37, 72, 1, 106, 50, 83, 138, 54, 89, 27, 12, 60, 99, 45, 147, 151, 98, 88, 95, 64, 55, 38, 5, 120, 112, 69, 118, 130, 59, 135, 39, 117, 113, 43, 61, 40, 29, 100, 143, 125, 136, 149, 91, 21, 33, 46, 103, 2, 13, 109, 144, 30, 146, 75], 3: [30, 85, 131, 64, 146, 41, 51, 62, 84, 59, 50, 37, 144, 129, 118, 139, 70, 122, 36, 142, 55, 87, 5, 8, 47, 17, 22, 27, 110, 137, 11, 145, 48, 76, 42, 97, 140, 106, 148, 68, 88, 89, 134, 33, 133, 6, 13, 79, 149, 75, 98, 4, 150, 147, 132, 116, 130, 35, 44, 1, 96, 72, 43, 54, 100, 151, 38, 80, 63, 120, 32, 67, 108, 77, 128, 21, 136, 66, 16, 56, 123, 29, 143, 125, 53, 61, 3, 18, 49, 91, 109, 92, 107, 7, 28, 95, 138, 46, 117, 10, 99, 60, 93, 119, 65, 101, 39, 24, 73, 127, 58, 12, 90, 113, 104, 112, 81, 74, 57, 25, 124, 14, 94, 20, 135, 86, 2, 71, 23, 19, 82, 102, 83, 78, 15, 103, 52, 69, 105, 34, 121, 141, 40, 126, 45, 115, 111, 9, 114, 31, 26]}

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
    return (math.exp(x[0]) - x[1]) ** 4 + 100 * (x[1] - x[2]) ** 6 + math.tan(x[2] - x[3]) ** 4 + x[0] ** 8

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

def miele_cantrell_ck_4D():
    PB = SimpleNamespace()
    PB.name = _NAME
    PB.nvar = 4
    PB.var_types = ['C', 'K', 'C', 'K']
    PB.cardinality = [math.inf, 151.0, math.inf, 151.0]
    PB.f_star = 0.0
    PB.f_star_justification = 'E'
    PB.n_minimizers_documented = 1
    PB.description = 'DFO Benchmark Suite v2.1.0 core (Python): CK adaptation of miele_cantrell (F = f o T).'
    PB.x0_idx = [-8.914628565311432, 38.0, -240.1208049664, 65.0]
    PB.x0_alt_idx = [[5.220441281795502, 43.0, -682.8470308864, 72.0], [-3.462530195713043, 47.0, 97.96104028160016, 79.0], [0.9799203276634216, 35.0, -344.7651856384, 85.0], [-7.804015934467316, 44.0, 423.1635771392, 71.0], [3.1001808047294617, 48.0, -19.56264878079992, 77.0], [-5.582790672779083, 38.0, -848.6681264128, 81.0], [-1.1403401494026184, 40.0, 318.5191964672001, 66.0], [-9.974758803844452, 49.0, -510.5862811648, 73.0], [6.987325012683868, 39.0, 656.6010417151999, 80.0]]
    PB.suite_release = '2.1.0'
    PB.base_function = 'miele_cantrell'
    PB.m = 151
    PB.pattern = 'CK'
    PB.strategy = 'moderate'
    PB.scale_factors_C = [1.0, 100.0]
    PB.lb_orig = [-10.0, -10.0, -10.0, -10.0]
    PB.ub_orig = [7.2313232421875, 7.144510491286898, 6.769932799999999, 6.4983756768013325]
    PB.lb_C_work = [-10.0, -10.0, -1000.0, -10.0]
    PB.ub_C_work = [7.2313232421875, 7.144510491286898, 676.9932799999999, 6.4983756768013325]
    PB.ProblemData = _problem_data()
    PB.func_idx = _objective
    PB.func_F = lambda xr: _objective(_cell2idx(xr))
    PB.idx2cell = _idx2cell
    PB.cell2idx = _cell2idx
    PB.decode = lambda xr: _decode(_cell2idx(xr))
    return PB
