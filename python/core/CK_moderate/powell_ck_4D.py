"""powell_ck_4D  DFO Benchmark Suite v2.1.0 core, Python version (self-contained; standard library only).
POWELL_CK_4D  DFO Benchmark Suite v2.1.0 core: powell, dimension 4, pattern CK, strategy moderate.

 Copyright (C) 2026 J. F. A. Madeira

 Objective F(z) = f(T(z)), where f is the original function (embedded below, translated from the MATLAB original) and T acts on the finite
 coordinates only: uniform grid of 151 levels on the interval of the domain, piecewise-affine map with |T - Id| <=
 half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
 Global optimum of F on this domain: f* = 0 (E: exact by structure); the optimal solution is listed in
 reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
 Variables: x1:C x2:K x3:C x4:K


Python version of the MATLAB wrapper powell_ck_4D.m: same instance, same values (level indices 1..m as in MATLAB).
Usage:  from powell_ck_4D import powell_ck_4D;  PB = powell_ck_4D();  f = PB.func_idx(PB.x0_idx)
"""
# Copyright (C) 2026 J. F. A. Madeira. LGPL-3.0-or-later.
import math
from types import SimpleNamespace

_NAME = 'powell_ck_4D'
_N = 4
_TYPES = ['C', 'K', 'C', 'K']
_CARD = [math.inf, 151.0, math.inf, 151.0]
_LB = [-6.0063232421875, -6.781547528323934, -10.0, -10.0]
_UB = [10.0, 10.0, 6.1539328, 6.469804248229904]
_LBW = [-6.0063232421874995, -6.781547528323934, -1000.0, -10.0]
_UBW = [10.0, 10.0, 615.39328, 6.469804248229904]
_TB = [math.nan, 0.04294846652779949, math.nan, -0.008318756073856903]
_TI = [math.nan, 0.0, math.nan, 0.0]
_CW = [8.881784197001252e-16, math.nan, 0.0, math.nan]
_CX = [0.0, math.nan, 0.0, math.nan]
_M = 151
_PERM = {1: [86, 57, 13, 44, 125, 103, 104, 143, 63, 119, 66, 118, 111, 107, 126, 98, 112, 151, 93, 144, 94, 56, 69, 145, 110, 73, 135, 35, 89, 122, 115, 102, 16, 137, 25, 139, 36, 149, 84, 127, 150, 88, 92, 38, 1, 87, 20, 130, 113, 121, 59, 15, 45, 49, 17, 29, 116, 55, 90, 108, 53, 71, 7, 81, 40, 133, 132, 23, 30, 65, 2, 39, 95, 14, 147, 76, 18, 48, 79, 43, 105, 136, 32, 142, 52, 67, 117, 28, 61, 5, 140, 31, 74, 96, 70, 134, 26, 8, 54, 85, 80, 146, 131, 51, 77, 128, 75, 68, 11, 72, 22, 101, 83, 3, 10, 124, 129, 42, 33, 27, 41, 114, 97, 37, 62, 24, 100, 120, 12, 60, 82, 4, 58, 47, 34, 19, 50, 109, 106, 138, 6, 91, 123, 64, 148, 46, 78, 141, 9, 21, 99], 3: [70, 112, 50, 9, 72, 129, 66, 148, 104, 36, 141, 136, 49, 138, 140, 59, 133, 108, 123, 33, 20, 86, 30, 100, 44, 58, 149, 111, 151, 76, 11, 13, 74, 53, 118, 28, 83, 150, 101, 16, 77, 143, 32, 103, 127, 14, 22, 25, 69, 65, 7, 135, 120, 126, 19, 31, 145, 110, 125, 115, 134, 71, 10, 21, 43, 121, 8, 27, 105, 68, 79, 139, 93, 107, 95, 4, 23, 119, 5, 48, 1, 109, 37, 98, 142, 73, 91, 113, 67, 12, 137, 39, 81, 116, 78, 55, 34, 17, 124, 94, 90, 56, 92, 61, 88, 38, 144, 147, 128, 85, 35, 146, 130, 18, 26, 54, 29, 60, 52, 47, 62, 15, 131, 63, 89, 51, 2, 24, 46, 97, 122, 45, 102, 96, 84, 75, 114, 41, 42, 82, 106, 57, 117, 87, 40, 3, 80, 6, 64, 99, 132]}

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
    return (x[0] + 10 * x[1]) ** 2 + 5 * (x[2] - x[3]) ** 2 + (x[1] - 2 * x[2]) ** 4 + 10 * (x[0] - x[3]) ** 4

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

def powell_ck_4D():
    PB = SimpleNamespace()
    PB.name = _NAME
    PB.nvar = 4
    PB.var_types = ['C', 'K', 'C', 'K']
    PB.cardinality = [math.inf, 151.0, math.inf, 151.0]
    PB.f_star = 0.0
    PB.f_star_justification = 'E'
    PB.n_minimizers_documented = 1
    PB.description = 'DFO Benchmark Suite v2.1.0 core (Python): CK adaptation of powell (F = f o T).'
    PB.x0_idx = [9.00742038488388, 77.0, 442.74004623359997, 69.0]
    PB.x0_alt_idx = [[1.1293081641197205, 82.0, -293.87928944639987, 75.0], [5.068364274501801, 71.0, -782.3742173184, 82.0], [-2.8097479462623594, 74.0, 31.783995801600213, 67.0], [7.03789232969284, 83.0, -469.11715246080007, 74.0], [-0.8402198910713192, 72.0, 345.0410606592002, 80.0], [3.098836219310761, 77.0, -143.45386721279988, 84.0], [-4.7792760014533995, 79.0, -880.0732028927999, 69.0], [8.022656357288358, 68.0, 182.2094180352001, 76.0], [0.14454413652420062, 78.0, -554.4099176448001, 83.0]]
    PB.suite_release = '2.1.0'
    PB.base_function = 'powell'
    PB.m = 151
    PB.pattern = 'CK'
    PB.strategy = 'moderate'
    PB.scale_factors_C = [1.0, 100.0]
    PB.lb_orig = [-6.0063232421875, -6.781547528323934, -10.0, -10.0]
    PB.ub_orig = [10.0, 10.0, 6.1539328, 6.469804248229904]
    PB.lb_C_work = [-6.0063232421874995, -6.781547528323934, -1000.0, -10.0]
    PB.ub_C_work = [10.0, 10.0, 615.39328, 6.469804248229904]
    PB.ProblemData = _problem_data()
    PB.func_idx = _objective
    PB.func_F = lambda xr: _objective(_cell2idx(xr))
    PB.idx2cell = _idx2cell
    PB.cell2idx = _cell2idx
    PB.decode = lambda xr: _decode(_cell2idx(xr))
    return PB
