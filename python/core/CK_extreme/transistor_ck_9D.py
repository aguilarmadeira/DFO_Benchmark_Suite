"""transistor_ck_9D  DFO Benchmark Suite v2.1.0 core, Python version (self-contained; standard library only).
TRANSISTOR_CK_9D  DFO Benchmark Suite v2.1.0 core: transistor, dimension 9, pattern CK, strategy extreme.

 Copyright (C) 2026 J. F. A. Madeira

 Objective F(z) = f(T(z)), where f is the original function (embedded below, translated from the MATLAB original) and T acts on the finite
 coordinates only: uniform grid of 25 levels on the interval of the domain, piecewise-affine map with |T - Id| <=
 half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
 Global optimum of F on this domain: f* = 1.60049326654736e-24 (E: exact by structure); the optimal solution is listed in
 reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
 Variables: x1:C x2:K x3:C x4:K x5:C x6:K x7:C x8:K x9:C


Python version of the MATLAB wrapper transistor_ck_9D.m: same instance, same values (level indices 1..m as in MATLAB).
Usage:  from transistor_ck_9D import transistor_ck_9D;  PB = transistor_ck_9D();  f = PB.func_idx(PB.x0_idx)
"""
# Copyright (C) 2026 J. F. A. Madeira. LGPL-3.0-or-later.
import math
from types import SimpleNamespace

_NAME = 'transistor_ck_9D'
_N = 9
_TYPES = ['C', 'K', 'C', 'K', 'C', 'K', 'C', 'K', 'C']
_CARD = [math.inf, 25.0, math.inf, 25.0, math.inf, 25.0, math.inf, 25.0, math.inf]
_LB = [-10.0, -10.0, -10.0, -10.0, -7.350085376681921, -6.992661321382304, -10.0, -10.0, -10.0]
_UB = [6.8211669921875, 7.03312503175329, 6.7027328, 6.180008329862558, 10.0, 10.0, 6.7876079069934505, 6.403313357018439, 7.316462562669516]
_LBW = [-10.0, -10.0, -10.0, -10.0, -735008537.6681921, -6.992661321382304, -1000000000.0, -10.0, -999999999.9999999]
_UBW = [6.8211669921875, 7.03312503175329, 6.7027328, 6.180008329862558, 1000000000.0, 10.0, 678760790.6993451, 6.403313357018439, 731646256.2669514]
_TB = [math.nan, 0.6457031448458057, math.nan, 2.1350062473969214, math.nan, 7.8759173348272125, math.nan, 0.9355422380122942, math.nan]
_TI = [math.nan, 0.4499874719815851, math.nan, 2.000068541623683, math.nan, 7.999692684218543, math.nan, 0.9999877234567898, math.nan]
_CW = [0.8999999526168558, math.nan, 1.0000064824653876, math.nan, 799997144.0507853, math.nan, 500003127.5929673, math.nan, 200005248.3486191]
_CX = [0.8999999526168565, math.nan, 1.000006482465387, math.nan, 7.999971440507853, math.nan, 5.000031275929672, math.nan, 2.00005248348619]
_M = 25
_PERM = {1: [2, 5, 22, 21, 12, 11, 25, 8, 17, 13, 10, 24, 3, 20, 1, 9, 15, 16, 4, 19, 23, 14, 7, 18, 6], 3: [2, 23, 17, 9, 25, 5, 22, 12, 15, 21, 4, 3, 19, 20, 11, 6, 14, 1, 10, 13, 24, 18, 16, 7, 8], 5: [17, 11, 5, 16, 13, 24, 3, 6, 8, 1, 12, 9, 18, 19, 22, 23, 14, 4, 15, 21, 7, 20, 25, 10, 2], 7: [18, 9, 4, 6, 21, 17, 10, 16, 24, 5, 25, 15, 23, 1, 3, 14, 13, 7, 20, 22, 8, 2, 11, 19, 12]}

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

_G = [[0.485, 0.369, 5.2095, 23.3037, 28.5132], [0.752, 1.254, 10.0677, 101.779, 111.8467],
      [0.869, 0.703, 22.9274, 111.461, 134.3884], [0.982, 1.455, 20.2153, 191.267, 211.4823]]
def _orig(x):
    s = 0.0
    for g in _G:
        alfa = (1 - x[0] * x[1]) * x[2] * (math.exp(x[4] * (g[0] - g[2] * x[6] * 0.001 - g[4] * x[7] * 0.001)) - 1) - g[4] + g[3] * x[1]
        beta = (1 - x[0] * x[1]) * x[3] * (math.exp(x[5] * (g[0] - g[1] - g[2] * x[6] * 0.001 + g[3] * x[8] * 0.001)) - 1) - g[4] * x[0] + g[3]
        s += alfa ** 2 + beta ** 2
    return (x[0] * x[2] - x[1] * x[3]) ** 2 + s

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

def transistor_ck_9D():
    PB = SimpleNamespace()
    PB.name = _NAME
    PB.nvar = 9
    PB.var_types = ['C', 'K', 'C', 'K', 'C', 'K', 'C', 'K', 'C']
    PB.cardinality = [math.inf, 25.0, math.inf, 25.0, math.inf, 25.0, math.inf, 25.0, math.inf]
    PB.f_star = 1.6004932665473596e-24
    PB.f_star_justification = 'E'
    PB.n_minimizers_documented = 1
    PB.description = 'DFO Benchmark Suite v2.1.0 core (Python): CK adaptation of transistor (F = f o T).'
    PB.x0_idx = [-3.4210181832313538, 2.0, -1.6299265392639999, 20.0, 336499364.63327587, 18.0, -145072359.38738823, 24.0, -698559482.9642966]
    PB.x0_alt_idx = [[0.9156889319419861, 4.0, -9.887757635584, 21.0, -207075812.0817024, 15.0, -615590089.6526027, 6.0, 587899834.0733439], [-7.560602247714996, 2.0, 1.7373443932160004, 19.0, -594226477.5117949, 10.0, 488435172.8651228, 14.0, 145986328.22071946], [3.084042489528656, 2.0, -6.648763691008, 20.0, 582867969.9069712, 7.0, 17917442.599908352, 22.0, -295927177.6319053], [-5.490810215473175, 3.0, 4.976338337791999, 21.0, 195717304.47687876, 4.0, -452600287.6653062, 6.0, -737840683.4845297], [-1.1541031002998352, 1.0, -3.2814927585280005, 18.0, -362078282.35490745, 23.0, -923118017.9305205, 13.0, 548618633.5531107], [-9.630394279956818, 3.0, -7.691014217728, 19.0, 971440676.3487444, 20.0, 279316181.6361382, 21.0, 106705127.70048606], [6.533695876598358, 3.0, 0.08577817395200071, 20.0, 427865499.63376594, 17.0, -289610485.67800975, 4.0, -335208378.1521385], [-2.090437591075897, 1.0, -4.323743285248001, 21.0, 40714834.20367348, 13.0, -760128215.9432241, 12.0, -777121884.0047631], [2.2462695240974426, 2.0, 3.453049106432001, 19.0, -502860342.51130474, 10.0, 442305983.6234348, 19.0, 434190788.5593878]]
    PB.suite_release = '2.1.0'
    PB.base_function = 'transistor'
    PB.m = 25
    PB.pattern = 'CK'
    PB.strategy = 'extreme'
    PB.scale_factors_C = [1.0, 1.0, 100000000.0, 100000000.0, 100000000.0]
    PB.lb_orig = [-10.0, -10.0, -10.0, -10.0, -7.350085376681921, -6.992661321382304, -10.0, -10.0, -10.0]
    PB.ub_orig = [6.8211669921875, 7.03312503175329, 6.7027328, 6.180008329862558, 10.0, 10.0, 6.7876079069934505, 6.403313357018439, 7.316462562669516]
    PB.lb_C_work = [-10.0, -10.0, -10.0, -10.0, -735008537.6681921, -6.992661321382304, -1000000000.0, -10.0, -999999999.9999999]
    PB.ub_C_work = [6.8211669921875, 7.03312503175329, 6.7027328, 6.180008329862558, 1000000000.0, 10.0, 678760790.6993451, 6.403313357018439, 731646256.2669514]
    PB.ProblemData = _problem_data()
    PB.func_idx = _objective
    PB.func_F = lambda xr: _objective(_cell2idx(xr))
    PB.idx2cell = _idx2cell
    PB.cell2idx = _cell2idx
    PB.decode = lambda xr: _decode(_cell2idx(xr))
    return PB
