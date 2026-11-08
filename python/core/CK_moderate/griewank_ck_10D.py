"""griewank_ck_10D  DFO Benchmark Suite v2.1.0 core, Python version (self-contained; standard library only).
GRIEWANK_CK_10D  DFO Benchmark Suite v2.1.0 core: griewank, dimension 10, pattern CK, strategy moderate.

 Copyright (C) 2026 J. F. A. Madeira

 Objective F(z) = f(T(z)), where f is the original function (embedded below, translated from the MATLAB original) and T acts on the finite
 coordinates only: uniform grid of 17 levels on the interval of the domain, piecewise-affine map with |T - Id| <=
 half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
 Global optimum of F on this domain: f* = 0 (E: exact by structure); the optimal solution is listed in
 reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
 Variables: x1:C x2:K x3:C x4:K x5:C x6:K x7:C x8:K x9:C x10:K


Python version of the MATLAB wrapper griewank_ck_10D.m: same instance, same values (level indices 1..m as in MATLAB).
Usage:  from griewank_ck_10D import griewank_ck_10D;  PB = griewank_ck_10D();  f = PB.func_idx(PB.x0_idx)
"""
# Copyright (C) 2026 J. F. A. Madeira. LGPL-3.0-or-later.
import math
from types import SimpleNamespace

_NAME = 'griewank_ck_10D'
_N = 10
_TYPES = ['C', 'K', 'C', 'K', 'C', 'K', 'C', 'K', 'C', 'K']
_CARD = [math.inf, 17.0, math.inf, 17.0, math.inf, 17.0, math.inf, 17.0, math.inf, 17.0]
_LB = [-270.0029296875, -256.7433826144388, -400.0, -400.0, -400.0, -400.0, -400.0, -400.0, -278.2615270814498, -287.5916191725778]
_UB = [400.0, 400.0, 264.525312, 285.07788421491045, 275.49101837306193, 285.0082280032212, 291.6565654146861, 262.337520430322, 400.0, 400.0]
_LBW = [-270.0029296875, -256.7433826144388, -40000.0, -400.0, -4.0, -400.0, -400.0, -400.0, -27826.15270814498, -287.5916191725778]
_UBW = [400.0, 400.0, 26452.531199999998, 285.07788421491045, 2.75491018373062, 285.0082280032212, 291.65656541468616, 262.337520430322, 40000.0, 400.0]
_TB = [math.nan, -10.46461413402426, math.nan, -14.643690129112883, math.nan, -14.68287174818812, math.nan, 13.960950268951251, math.nan, 13.22971421542502]
_TI = [math.nan, 0.0, math.nan, 0.0, math.nan, 0.0, math.nan, 0.0, math.nan, 0.0]
_CW = [0.0, math.nan, 0.0, math.nan, 8.881784197001252e-16, math.nan, 5.684341886080802e-14, math.nan, 0.0, math.nan]
_CX = [0.0, math.nan, 0.0, math.nan, 0.0, math.nan, 0.0, math.nan, 0.0, math.nan]
_M = 17
_PERM = {1: [15, 14, 1, 9, 13, 3, 12, 10, 2, 7, 11, 16, 8, 5, 17, 6, 4], 3: [15, 13, 12, 7, 17, 1, 11, 5, 2, 9, 4, 3, 16, 8, 10, 14, 6], 5: [13, 5, 9, 12, 6, 4, 16, 3, 15, 14, 7, 11, 8, 10, 17, 2, 1], 7: [3, 5, 6, 4, 12, 14, 7, 15, 1, 2, 10, 13, 16, 11, 9, 17, 8], 9: [10, 6, 1, 11, 12, 8, 5, 7, 13, 3, 2, 17, 9, 14, 16, 15, 4]}

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
    n = len(x)
    s = 0.0
    for v in x: s += v ** 2
    p = 1.0
    for i in range(n): p *= math.cos(x[i] / math.sqrt(i + 1))
    return 1 + 1 / 4000 * s - p

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

def griewank_ck_10D():
    PB = SimpleNamespace()
    PB.name = _NAME
    PB.nvar = 10
    PB.var_types = ['C', 'K', 'C', 'K', 'C', 'K', 'C', 'K', 'C', 'K']
    PB.cardinality = [math.inf, 17.0, math.inf, 17.0, math.inf, 17.0, math.inf, 17.0, math.inf, 17.0]
    PB.f_star = 0.0
    PB.f_star_justification = 'E'
    PB.n_minimizers_documented = 1
    PB.description = 'DFO Benchmark Suite v2.1.0 core (Python): CK adaptation of griewank (F = f o T).'
    PB.x0_idx = [64.34423542022705, 15.0, -2446.3455682559943, 3.0, 2.404730767744586, 11.0, -296.66681884502754, 14.0, 4361.585085627441, 5.0]
    PB.x0_alt_idx = [[222.3576216697693, 14.0, -32748.699795456, 1.0, 0.8974367598046946, 9.0, 198.7411709156646, 3.0, -15890.935074534525, 8.0], [-107.40944528579712, 14.0, 10950.484721664005, 2.0, -1.2188649281109094, 7.0, 4.885870574524176, 8.0, 34497.870245505124, 11.0], [308.7251868247986, 15.0, -19351.869505535997, 3.0, -3.335166616026513, 4.0, -229.51432918437106, 13.0, 17188.73864701822, 13.0], [-21.041880130767822, 14.0, 23836.959571968, 1.0, 1.8566238557664434, 2.0, 265.89366057632105, 2.0, -120.39295146868244, 16.0], [139.91585493087769, 14.0, -6465.394655232005, 2.0, -0.3150421666924079, 16.0, 72.03836023518073, 7.0, -17429.524549955582, 3.0], [-189.85121202468872, 15.0, -39319.526080512, 2.0, -1.8223361746322984, 12.0, -121.81694010595976, 13.0, 32959.280770084064, 6.0], [351.90896940231323, 13.0, 6931.435634688001, 3.0, -3.9386378625479015, 10.0, -356.21713986485497, 1.0, 15650.14917159716, 9.0], [22.141902446746826, 14.0, -25922.695790592, 1.0, 1.2531526092450544, 7.0, 139.19084989583712, 6.0, -1658.9824268897355, 11.0], [181.1367383003235, 15.0, 20328.265924607993, 2.0, -0.8631490786705496, 5.0, -54.66445044530326, 11.0, -21911.50258705171, 14.0]]
    PB.suite_release = '2.1.0'
    PB.base_function = 'griewank'
    PB.m = 17
    PB.pattern = 'CK'
    PB.strategy = 'moderate'
    PB.scale_factors_C = [1.0, 100.0, 0.01, 1.0, 100.0]
    PB.lb_orig = [-270.0029296875, -256.7433826144388, -400.0, -400.0, -400.0, -400.0, -400.0, -400.0, -278.2615270814498, -287.5916191725778]
    PB.ub_orig = [400.0, 400.0, 264.525312, 285.07788421491045, 275.49101837306193, 285.0082280032212, 291.6565654146861, 262.337520430322, 400.0, 400.0]
    PB.lb_C_work = [-270.0029296875, -256.7433826144388, -40000.0, -400.0, -4.0, -400.0, -400.0, -400.0, -27826.15270814498, -287.5916191725778]
    PB.ub_C_work = [400.0, 400.0, 26452.531199999998, 285.07788421491045, 2.75491018373062, 285.0082280032212, 291.65656541468616, 262.337520430322, 40000.0, 400.0]
    PB.ProblemData = _problem_data()
    PB.func_idx = _objective
    PB.func_F = lambda xr: _objective(_cell2idx(xr))
    PB.idx2cell = _idx2cell
    PB.cell2idx = _cell2idx
    PB.decode = lambda xr: _decode(_cell2idx(xr))
    return PB
