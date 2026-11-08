"""fifteenn_local_minima_ck_4D  DFO Benchmark Suite v2.1.0 core, Python version (self-contained; standard library only).
FIFTEENN_LOCAL_MINIMA_CK_4D  DFO Benchmark Suite v2.1.0 core: fifteenn_local_minima, dimension 4, pattern CK, strategy baseline.

 Copyright (C) 2026 J. F. A. Madeira

 Objective F(z) = f(T(z)), where f is the original function (embedded below, translated from the MATLAB original) and T acts on the finite
 coordinates only: uniform grid of 151 levels on the interval of the domain, piecewise-affine map with |T - Id| <=
 half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
 Global optimum of F on this domain: f* = 1.349783804395672e-32 (E: exact by structure); the optimal solution is listed in
 reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
 Variables: x1:C x2:K x3:C x4:K


Python version of the MATLAB wrapper fifteenn_local_minima_ck_4D.m: same instance, same values (level indices 1..m as in MATLAB).
Usage:  from fifteenn_local_minima_ck_4D import fifteenn_local_minima_ck_4D;  PB = fifteenn_local_minima_ck_4D();  f = PB.func_idx(PB.x0_idx)
"""
# Copyright (C) 2026 J. F. A. Madeira. LGPL-3.0-or-later.
import math
from types import SimpleNamespace

_NAME = 'fifteenn_local_minima_ck_4D'
_N = 4
_TYPES = ['C', 'K', 'C', 'K']
_CARD = [math.inf, 151.0, math.inf, 151.0]
_LB = [-10.0, -10.0, -10.0, -10.0]
_UB = [7.275073242187499, 6.176609256718996, 6.949132799999999, 6.95551853394419]
_LBW = [-10.0, -10.0, -10.0, -10.0]
_UBW = [7.275073242187499, 6.176609256718996, 6.949132799999999, 6.95551853394419]
_TB = [math.nan, 1.0000942945689175, math.nan, 0.964568651950577]
_TI = [math.nan, 1.0, math.nan, 1.0]
_CW = [1.0, math.nan, 1.0, math.nan]
_CX = [1.0, math.nan, 1.0, math.nan]
_M = 151
_PERM = {1: [137, 48, 20, 33, 103, 71, 47, 84, 89, 58, 12, 144, 53, 61, 125, 45, 11, 62, 117, 116, 145, 115, 55, 135, 108, 110, 14, 82, 91, 114, 59, 32, 3, 72, 52, 97, 16, 24, 26, 83, 35, 90, 65, 67, 9, 140, 86, 77, 54, 41, 87, 131, 31, 34, 111, 106, 10, 120, 30, 151, 5, 75, 134, 80, 85, 138, 66, 76, 18, 44, 78, 6, 130, 95, 68, 88, 27, 49, 142, 94, 127, 143, 132, 43, 74, 64, 128, 126, 92, 107, 21, 37, 69, 148, 113, 96, 57, 19, 129, 139, 22, 150, 147, 50, 133, 4, 40, 136, 46, 119, 81, 1, 42, 38, 101, 63, 15, 36, 118, 122, 93, 23, 39, 56, 104, 141, 100, 149, 28, 29, 8, 79, 17, 105, 60, 7, 123, 124, 51, 98, 13, 102, 121, 2, 109, 112, 99, 25, 146, 73, 70], 3: [9, 100, 24, 51, 69, 14, 34, 123, 62, 85, 99, 23, 55, 89, 116, 35, 12, 113, 39, 119, 82, 131, 67, 5, 104, 42, 112, 88, 22, 38, 139, 111, 4, 124, 107, 97, 93, 25, 27, 31, 84, 140, 102, 103, 11, 129, 149, 91, 78, 132, 3, 75, 66, 36, 71, 87, 13, 144, 73, 86, 96, 142, 108, 98, 45, 109, 143, 106, 151, 90, 53, 10, 33, 122, 72, 17, 94, 121, 77, 30, 148, 74, 127, 135, 79, 70, 7, 19, 92, 47, 29, 44, 76, 147, 137, 145, 130, 105, 37, 60, 117, 83, 125, 133, 28, 128, 52, 134, 59, 49, 95, 15, 138, 118, 2, 80, 32, 54, 18, 43, 141, 46, 64, 81, 6, 40, 1, 16, 150, 65, 48, 120, 61, 20, 110, 63, 50, 146, 126, 41, 114, 68, 115, 26, 8, 57, 101, 56, 58, 136, 21]}

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
    d = len(x)
    s = 0.0
    for i in range(d - 1): s += ((x[i] - 1) ** 2) * (1 + 10 * math.sin(3 * math.pi * x[i + 1]) ** 2)
    f = 0.1 * s
    return f + 0.1 * (math.sin(3 * math.pi * x[0]) ** 2 + (x[d - 1] - 1) ** 2 * (1 + math.sin(2 * math.pi * x[d - 1]) ** 2))

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

def fifteenn_local_minima_ck_4D():
    PB = SimpleNamespace()
    PB.name = _NAME
    PB.nvar = 4
    PB.var_types = ['C', 'K', 'C', 'K']
    PB.cardinality = [math.inf, 151.0, math.inf, 151.0]
    PB.f_star = 1.3497838043956716e-32
    PB.f_star_justification = 'E'
    PB.n_minimizers_documented = 1
    PB.description = 'DFO Benchmark Suite v2.1.0 core (Python): CK adaptation of fifteenn_local_minima (F = f o T).'
    PB.x0_idx = [-7.857486033439637, 142.0, -4.4894979440639995, 40.0]
    PB.x0_alt_idx = [[5.933893239498138, 146.0, 4.052864987136003, 23.0], [-2.5686818718910223, 136.0, -1.0725527715840002, 30.0], [1.6826056838035583, 145.0, -8.801357328384, 36.0], [-6.819969427585602, 148.0, 2.2142230609920013, 43.0], [3.8082494616508473, 137.0, -5.514581495808, 28.0], [-4.694325649738312, 142.0, 5.6311682334720015, 32.0], [-0.4430380940437324, 151.0, -2.0976363233279995, 38.0], [-8.945613205432892, 138.0, -7.223054082048, 24.0], [6.4906094670295715, 142.0, 1.3193088491520015, 31.0]]
    PB.suite_release = '2.1.0'
    PB.base_function = 'fifteenn_local_minima'
    PB.m = 151
    PB.pattern = 'CK'
    PB.strategy = 'baseline'
    PB.scale_factors_C = [1.0, 1.0]
    PB.lb_orig = [-10.0, -10.0, -10.0, -10.0]
    PB.ub_orig = [7.275073242187499, 6.176609256718996, 6.949132799999999, 6.95551853394419]
    PB.lb_C_work = [-10.0, -10.0, -10.0, -10.0]
    PB.ub_C_work = [7.275073242187499, 6.176609256718996, 6.949132799999999, 6.95551853394419]
    PB.ProblemData = _problem_data()
    PB.func_idx = _objective
    PB.func_F = lambda xr: _objective(_cell2idx(xr))
    PB.idx2cell = _idx2cell
    PB.cell2idx = _cell2idx
    PB.decode = lambda xr: _decode(_cell2idx(xr))
    return PB
