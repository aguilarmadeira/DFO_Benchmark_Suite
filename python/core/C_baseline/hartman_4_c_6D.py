"""hartman_4_c_6D  DFO Benchmark Suite v2.1.0 core, Python version (self-contained; standard library only).
HARTMAN_4_C_6D  DFO Benchmark Suite v2.1.0 core: hartman_4, dimension 6, pattern C, strategy baseline.

 Copyright (C) 2026 J. F. A. Madeira

 Objective F(z) = f(T(z)), where f is the original function (embedded below, translated from the MATLAB original) and T acts on the finite
 coordinates only: uniform grid of 0 levels on the interval of the domain, piecewise-affine map with |T - Id| <=
 half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
 Global optimum of F on this domain: f* = -3.322368011391339 (L: literature value at a refined point); the optimal solution is listed in
 reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
 Variables: x1:C x2:C x3:C x4:C x5:C x6:C


Python version of the MATLAB wrapper hartman_4_c_6D.m: same instance, same values (level indices 1..m as in MATLAB).
Usage:  from hartman_4_c_6D import hartman_4_c_6D;  PB = hartman_4_c_6D();  f = PB.func_idx(PB.x0_idx)
"""
# Copyright (C) 2026 J. F. A. Madeira. LGPL-3.0-or-later.
import math
from types import SimpleNamespace

_NAME = 'hartman_4_c_6D'
_N = 6
_TYPES = ['C', 'C', 'C', 'C', 'C', 'C']
_CARD = [math.inf, math.inf, math.inf, math.inf, math.inf, math.inf]
_LB = [0.0, 0.14795966062084032, 0.14470336, 0.0, 0.18109077248821803, 0.16527817653443508]
_UB = [0.848441162109375, 1.0, 1.0, 0.816347355268638, 1.0, 1.0]
_LBW = [0.0, 0.14795966062084032, 0.14470336, 0.0, 0.18109077248821803, 0.16527817653443508]
_UBW = [0.848441162109375, 1.0, 1.0, 0.816347355268638, 1.0, 1.0]
_TB = [math.nan, math.nan, math.nan, math.nan, math.nan, math.nan]
_TI = [math.nan, math.nan, math.nan, math.nan, math.nan, math.nan]
_CW = [0.20169, 0.150011, 0.476874, 0.275332, 0.311652, 0.6573]
_CX = [0.20169, 0.150011, 0.476874, 0.275332, 0.311652, 0.6573]
_M = 0
_PERM = {}

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

def hartman_4_c_6D():
    PB = SimpleNamespace()
    PB.name = _NAME
    PB.nvar = 6
    PB.var_types = ['C', 'C', 'C', 'C', 'C', 'C']
    PB.cardinality = [math.inf, math.inf, math.inf, math.inf, math.inf, math.inf]
    PB.f_star = -3.322368011391339
    PB.f_star_justification = 'L'
    PB.n_minimizers_documented = 1
    PB.description = 'DFO Benchmark Suite v2.1.0 core (Python): C adaptation of hartman_4 (F = f o T).'
    PB.x0_idx = [0.21252456843853, 0.17055606651201144, 0.32698417991679996, 0.5814052384462187, 0.44319094598184616, 0.9380702516045121]
    PB.x0_alt_idx = [[0.4312633055448532, 0.19510043842828353, 0.9222706413568, 0.4865443837523619, 0.18662810009723835, 0.8343474214697377], [0.0037285012006759645, 0.21964481034455563, 0.4994119825408, 0.5222447054113402, 0.8160376716558807, 0.7306245913349634], [0.838912770152092, 0.15068681305598167, 0.27361366958079997, 0.5579450270703186, 0.5594748257712728, 0.5630723272710971], [0.421320635676384, 0.20328189573370756, 0.6652711069696, 0.4794043194205663, 0.3767430146736025, 0.4593494971363228], [0.6251453679800034, 0.22782626764997965, 0.4394727940096, 0.5151046410795447, 0.9256096391919494, 0.35562666700154844], [0.20755323350429536, 0.1682185072818903, 0.8376989095936, 0.5344848156944184, 0.7428778280942789, 0.18807440293768218], [0.734514736533165, 0.18029589663751625, 0.6119005966336, 0.5701851373533968, 0.48631498220967107, 0.9141342138811026], [0.316922602057457, 0.23289097931524216, 0.1890419378176, 0.49164442970364447, 0.3035831711120007, 0.7465819498172364], [0.518261666893959, 0.1732832189471528, 0.7843283992576, 0.5273447513626227, 0.8591617078837055, 0.6428591196824621]]
    PB.suite_release = '2.1.0'
    PB.base_function = 'hartman_4'
    PB.m = 0
    PB.pattern = 'C'
    PB.strategy = 'baseline'
    PB.scale_factors_C = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
    PB.lb_orig = [0.0, 0.14795966062084032, 0.14470336, 0.0, 0.18109077248821803, 0.16527817653443508]
    PB.ub_orig = [0.848441162109375, 1.0, 1.0, 0.816347355268638, 1.0, 1.0]
    PB.lb_C_work = [0.0, 0.14795966062084032, 0.14470336, 0.0, 0.18109077248821803, 0.16527817653443508]
    PB.ub_C_work = [0.848441162109375, 1.0, 1.0, 0.816347355268638, 1.0, 1.0]
    PB.ProblemData = _problem_data()
    PB.func_idx = _objective
    PB.func_F = lambda xr: _objective(_cell2idx(xr))
    PB.idx2cell = _idx2cell
    PB.cell2idx = _cell2idx
    PB.decode = lambda xr: _decode(_cell2idx(xr))
    return PB
