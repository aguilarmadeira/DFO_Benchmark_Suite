"""tenn_local_minima_ck_8D  DFO Benchmark Suite v2.1.0 core, Python version (self-contained; standard library only).
TENN_LOCAL_MINIMA_CK_8D  DFO Benchmark Suite v2.1.0 core: tenn_local_minima, dimension 8, pattern CK, strategy extreme.

 Copyright (C) 2026 J. F. A. Madeira

 Objective F(z) = f(T(z)), where f is the original function (embedded below, translated from the MATLAB original) and T acts on the finite
 coordinates only: uniform grid of 25 levels on the interval of the domain, piecewise-affine map with |T - Id| <=
 half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
 Global optimum of F on this domain: f* = 5.889542894199896e-32 (E: exact by structure); the optimal solution is listed in
 reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
 Variables: x1:C x2:K x3:C x4:K x5:C x6:K x7:C x8:K


Python version of the MATLAB wrapper tenn_local_minima_ck_8D.m: same instance, same values (level indices 1..m as in MATLAB).
Usage:  from tenn_local_minima_ck_8D import tenn_local_minima_ck_8D;  PB = tenn_local_minima_ck_8D();  f = PB.func_idx(PB.x0_idx)
"""
# Copyright (C) 2026 J. F. A. Madeira. LGPL-3.0-or-later.
import math
from types import SimpleNamespace

_NAME = 'tenn_local_minima_ck_8D'
_N = 8
_TYPES = ['C', 'K', 'C', 'K', 'C', 'K', 'C', 'K']
_CARD = [math.inf, 25.0, math.inf, 25.0, math.inf, 25.0, math.inf, 25.0]
_LB = [-10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0]
_UB = [7.1711669921875, 6.722013920642178, 7.2627328, 6.580008329862557, 6.216201079161259, 7.20804593676692, 6.952313789346391, 6.55068177807107]
_LBW = [-10.0, -10.0, -10.0, -10.0, -1000000000.0, -10.0, -1000000000.0, -10.0]
_UBW = [7.1711669921875005, 6.722013920642178, 7.262732800000001, 6.580008329862557, 621620107.9161258, 7.20804593676692, 695231378.9346392, 6.55068177807107]
_TB = [math.nan, 1.1480092804281181, math.nan, 1.0533388865750393, math.nan, 0.7550287104793281, math.nan, 1.0337878520473804]
_TI = [math.nan, 1.0, math.nan, 1.0, math.nan, 1.0, math.nan, 1.0]
_CW = [1.0000000000000018, math.nan, 1.0, math.nan, 100000000.0, math.nan, 100000000.0, math.nan]
_CX = [1.0, math.nan, 1.0, math.nan, 1.0, math.nan, 1.0, math.nan]
_M = 25
_PERM = {1: [4, 11, 19, 14, 3, 5, 10, 25, 7, 23, 17, 2, 6, 22, 18, 21, 24, 12, 16, 1, 15, 20, 13, 8, 9], 3: [18, 6, 4, 23, 5, 8, 12, 1, 16, 14, 19, 9, 13, 24, 17, 11, 7, 20, 25, 21, 15, 3, 22, 10, 2], 5: [17, 12, 24, 13, 1, 2, 11, 14, 25, 5, 19, 15, 21, 22, 3, 23, 6, 18, 8, 9, 16, 4, 7, 20, 10], 7: [10, 17, 11, 13, 16, 4, 25, 7, 8, 23, 2, 22, 24, 12, 19, 18, 5, 20, 21, 3, 14, 9, 15, 1, 6]}

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
    for i in range(d - 1): s += ((x[i] - 1) ** 2) * (1 + 10 * math.sin(math.pi * x[i + 1]) ** 2)
    f = math.pi / d * s
    return f + math.pi / d * (10 * math.sin(math.pi * x[0]) ** 2 + (x[d - 1] - 1) ** 2)

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

def tenn_local_minima_ck_8D():
    PB = SimpleNamespace()
    PB.name = _NAME
    PB.nvar = 8
    PB.var_types = ['C', 'K', 'C', 'K', 'C', 'K', 'C', 'K']
    PB.cardinality = [math.inf, 25.0, math.inf, 25.0, math.inf, 25.0, math.inf, 25.0]
    PB.f_star = 5.889542894199896e-32
    PB.f_star_justification = 'E'
    PB.n_minimizers_documented = 1
    PB.description = 'DFO Benchmark Suite v2.1.0 core (Python): CK adaptation of tenn_local_minima (F = f o T).'
    PB.x0_idx = [-7.576920282840728, 14.0, -8.254392459264, 13.0, -293358630.6601405, 14.0, -336123514.5389283, 21.0]
    PB.x0_alt_idx = [[3.2892088294029254, 16.0, 3.7604695695360046, 14.0, -801409408.2717291, 11.0, -811257568.8424084, 4.0], [-5.363449537754058, 17.0, -4.774225526783999, 12.0, 444959765.581161, 7.0, 303599460.53634596, 11.0], [-0.9365080475807179, 15.0, 7.240636502016002, 13.0, -63091012.030427694, 4.0, -171534593.76713443, 19.0], [-9.68977872133255, 15.0, -1.4266363822079988, 14.0, -424940127.0199765, 25.0, -646668648.0706146, 3.0], [6.810639560222626, 17.0, -9.961331478528, 11.0, -946281964.8699323, 20.0, 567562824.0382793, 11.0], [-1.8420188069343553, 15.0, 2.0535305502720025, 12.0, 300087208.9829576, 17.0, 92428769.73479915, 19.0], [2.584922683238984, 16.0, -6.481164546047999, 13.0, -207963568.62863088, 14.0, -482079727.2988208, 1.0], [-6.218654143810272, 16.0, 5.533697482752002, 14.0, -569812683.6181796, 9.0, -957213781.602301, 9.0], [4.6474749684333805, 15.0, -3.000997613567999, 12.0, 530354827.61267066, 6.0, 257017690.50659275, 17.0]]
    PB.suite_release = '2.1.0'
    PB.base_function = 'tenn_local_minima'
    PB.m = 25
    PB.pattern = 'CK'
    PB.strategy = 'extreme'
    PB.scale_factors_C = [1.0, 1.0, 100000000.0, 100000000.0]
    PB.lb_orig = [-10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0]
    PB.ub_orig = [7.1711669921875, 6.722013920642178, 7.2627328, 6.580008329862557, 6.216201079161259, 7.20804593676692, 6.952313789346391, 6.55068177807107]
    PB.lb_C_work = [-10.0, -10.0, -10.0, -10.0, -1000000000.0, -10.0, -1000000000.0, -10.0]
    PB.ub_C_work = [7.1711669921875005, 6.722013920642178, 7.262732800000001, 6.580008329862557, 621620107.9161258, 7.20804593676692, 695231378.9346392, 6.55068177807107]
    PB.ProblemData = _problem_data()
    PB.func_idx = _objective
    PB.func_F = lambda xr: _objective(_cell2idx(xr))
    PB.idx2cell = _idx2cell
    PB.cell2idx = _cell2idx
    PB.decode = lambda xr: _decode(_cell2idx(xr))
    return PB
