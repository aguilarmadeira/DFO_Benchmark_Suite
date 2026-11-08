"""wood_ck_4D  DFO Benchmark Suite v2.1.0 core, Python version (self-contained; standard library only).
WOOD_CK_4D  DFO Benchmark Suite v2.1.0 core: wood, dimension 4, pattern CK, strategy baseline.

 Copyright (C) 2026 J. F. A. Madeira

 Objective F(z) = f(T(z)), where f is the original function (embedded below, translated from the MATLAB original) and T acts on the finite
 coordinates only: uniform grid of 151 levels on the interval of the domain, piecewise-affine map with |T - Id| <=
 half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
 Global optimum of F on this domain: f* = 0 (E: exact by structure); the optimal solution is listed in
 reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
 Variables: x1:C x2:K x3:C x4:K


Python version of the MATLAB wrapper wood_ck_4D.m: same instance, same values (level indices 1..m as in MATLAB).
Usage:  from wood_ck_4D import wood_ck_4D;  PB = wood_ck_4D();  f = PB.func_idx(PB.x0_idx)
"""
# Copyright (C) 2026 J. F. A. Madeira. LGPL-3.0-or-later.
import math
from types import SimpleNamespace

_NAME = 'wood_ck_4D'
_N = 4
_TYPES = ['C', 'K', 'C', 'K']
_CARD = [math.inf, 151.0, math.inf, 151.0]
_LB = [-10.0, -10.0, -10.0, -10.0]
_UB = [6.1211669921875, 6.566458365086623, 6.4227327999999995, 7.3514369012911285]
_LBW = [-10.0, -10.0, -10.0, -10.0]
_UBW = [6.1211669921875, 6.566458365086623, 6.4227327999999995, 7.3514369012911285]
_TB = [math.nan, 1.0443055767244172, math.nan, 0.9892433708177162]
_TI = [math.nan, 1.0, math.nan, 1.0]
_CW = [1.0, math.nan, 1.0, math.nan]
_CX = [1.0, math.nan, 1.0, math.nan]
_M = 151
_PERM = {1: [30, 135, 36, 60, 80, 141, 142, 61, 111, 151, 55, 102, 78, 62, 148, 43, 45, 70, 53, 116, 33, 103, 34, 118, 100, 42, 82, 146, 128, 17, 115, 139, 39, 73, 86, 88, 108, 68, 2, 83, 46, 145, 26, 98, 105, 147, 64, 28, 38, 87, 110, 91, 72, 127, 48, 94, 93, 112, 51, 15, 130, 3, 10, 89, 52, 150, 25, 49, 92, 96, 41, 16, 132, 84, 47, 6, 138, 131, 27, 1, 44, 74, 124, 4, 54, 125, 50, 134, 24, 90, 31, 58, 123, 8, 37, 104, 109, 114, 20, 59, 79, 63, 5, 149, 121, 101, 56, 75, 137, 57, 22, 40, 76, 35, 106, 117, 85, 12, 32, 119, 81, 65, 18, 144, 136, 29, 99, 13, 19, 122, 67, 107, 23, 71, 95, 11, 120, 140, 77, 126, 21, 129, 7, 143, 133, 66, 14, 69, 97, 9, 113], 3: [98, 25, 82, 11, 150, 140, 59, 123, 22, 62, 41, 118, 37, 100, 90, 88, 66, 148, 86, 113, 114, 5, 81, 110, 151, 138, 39, 108, 14, 129, 131, 48, 104, 84, 63, 103, 137, 78, 4, 74, 46, 24, 64, 111, 127, 94, 93, 101, 96, 2, 126, 124, 116, 144, 75, 76, 21, 42, 112, 18, 105, 80, 83, 27, 55, 146, 130, 52, 61, 35, 142, 20, 85, 8, 134, 135, 34, 71, 32, 6, 49, 79, 139, 149, 60, 147, 56, 77, 30, 106, 125, 65, 26, 132, 109, 97, 40, 128, 29, 68, 122, 73, 15, 3, 92, 115, 133, 87, 45, 70, 36, 54, 91, 50, 7, 117, 141, 145, 51, 19, 102, 89, 43, 33, 38, 57, 143, 44, 53, 28, 107, 16, 136, 119, 10, 67, 47, 58, 13, 69, 121, 95, 1, 23, 31, 12, 17, 99, 9, 72, 120]}

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
    return 100 * (x[0] ** 2 - x[1]) ** 2 + (1 - x[0]) ** 2 + 90 * (x[2] ** 2 - x[3]) ** 2 + (1 - x[2]) ** 2 + 10.1 * ((1 - x[1]) ** 2 + (1 - x[3]) ** 2) + 19.8 * (1 - x[1]) * (1 - x[3])

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

def wood_ck_4D():
    PB = SimpleNamespace()
    PB.name = _NAME
    PB.nvar = 4
    PB.var_types = ['C', 'K', 'C', 'K']
    PB.cardinality = [math.inf, 151.0, math.inf, 151.0]
    PB.f_star = 0.0
    PB.f_star_justification = 'E'
    PB.n_minimizers_documented = 1
    PB.description = 'DFO Benchmark Suite v2.1.0 core (Python): CK adaptation of wood (F = f o T).'
    PB.x0_idx = [4.365786015987396, 57.0, 1.5143064207360002, 142.0]
    PB.x0_alt_idx = [[-3.5688508629798887, 67.0, -6.605092675584, 149.0], [0.3984675765037533, 56.0, 4.825129353215997, 134.0], [-7.536169302463531, 59.0, -3.420396331008001, 141.0], [2.3821267962455757, 63.0, -7.755997790208, 147.0], [-5.55251008272171, 53.0, -0.10957339852799919, 130.0], [-1.5851916432380673, 62.0, -4.445174857728, 136.0], [-9.519828522205353, 64.0, 3.201249533951998, 143.0], [5.640995514392854, 54.0, -1.1343519252480014, 150.0], [-2.293641364574432, 58.0, -9.253751021568, 135.0]]
    PB.suite_release = '2.1.0'
    PB.base_function = 'wood'
    PB.m = 151
    PB.pattern = 'CK'
    PB.strategy = 'baseline'
    PB.scale_factors_C = [1.0, 1.0]
    PB.lb_orig = [-10.0, -10.0, -10.0, -10.0]
    PB.ub_orig = [6.1211669921875, 6.566458365086623, 6.4227327999999995, 7.3514369012911285]
    PB.lb_C_work = [-10.0, -10.0, -10.0, -10.0]
    PB.ub_C_work = [6.1211669921875, 6.566458365086623, 6.4227327999999995, 7.3514369012911285]
    PB.ProblemData = _problem_data()
    PB.func_idx = _objective
    PB.func_F = lambda xr: _objective(_cell2idx(xr))
    PB.idx2cell = _idx2cell
    PB.cell2idx = _cell2idx
    PB.decode = lambda xr: _decode(_cell2idx(xr))
    return PB
