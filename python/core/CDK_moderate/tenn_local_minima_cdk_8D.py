"""tenn_local_minima_cdk_8D  DFO Benchmark Suite v2.1.0 core, Python version (self-contained; standard library only).
TENN_LOCAL_MINIMA_CDK_8D  DFO Benchmark Suite v2.1.0 core: tenn_local_minima, dimension 8, pattern CDK, strategy moderate.

 Copyright (C) 2026 J. F. A. Madeira

 Objective F(z) = f(T(z)), where f is the original function (embedded below, translated from the MATLAB original) and T acts on the finite
 coordinates only: uniform grid of 17 levels on the interval of the domain, piecewise-affine map with |T - Id| <=
 half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
 Global optimum of F on this domain: f* = 5.889542894199896e-32 (E: exact by structure); the optimal solution is listed in
 reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
 Variables: x1:C x2:K x3:D x4:C x5:D x6:K x7:C x8:K


Python version of the MATLAB wrapper tenn_local_minima_cdk_8D.m: same instance, same values (level indices 1..m as in MATLAB).
Usage:  from tenn_local_minima_cdk_8D import tenn_local_minima_cdk_8D;  PB = tenn_local_minima_cdk_8D();  f = PB.func_idx(PB.x0_idx)
"""
# Copyright (C) 2026 J. F. A. Madeira. LGPL-3.0-or-later.
import math
from types import SimpleNamespace

_NAME = 'tenn_local_minima_cdk_8D'
_N = 8
_TYPES = ['C', 'K', 'D', 'C', 'D', 'K', 'C', 'K']
_CARD = [math.inf, 17.0, 17.0, math.inf, 17.0, 17.0, math.inf, 17.0]
_LB = [-10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0]
_UB = [7.1711669921875, 6.722013920642178, 7.2627328, 6.580008329862557, 6.216201079161259, 7.20804593676692, 6.952313789346391, 6.55068177807107]
_LBW = [-10.0, -10.0, -10.0, -1000.0, -10.0, -10.0, -0.09999999999999999, -10.0]
_UBW = [7.1711669921875005, 6.722013920642178, 7.2627328, 658.0008329862558, 6.216201079161259, 7.20804593676692, 0.06952313789346391, 6.55068177807107]
_TB = [math.nan, 1.4963845704414958, 0.7892080000000021, math.nan, 1.148638241923365, 0.7550287104793263, math.nan, 1.3785937224238598]
_TI = [math.nan, 1.0, 1.0, math.nan, 1.0, 1.0, math.nan, 1.0]
_CW = [1.0000000000000018, math.nan, math.nan, 100.0, math.nan, math.nan, 0.010000000000000009, math.nan]
_CX = [1.0, math.nan, math.nan, 1.0, math.nan, math.nan, 1.0, math.nan]
_M = 17
_PERM = {1: [11, 2, 16, 10, 3, 12, 17, 8, 4, 7, 14, 1, 15, 13, 9, 5, 6], 5: [11, 8, 12, 1, 9, 17, 2, 4, 14, 16, 15, 6, 10, 3, 5, 13, 7], 7: [17, 11, 5, 13, 16, 14, 7, 15, 3, 12, 8, 2, 9, 6, 10, 1, 4]}

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

def tenn_local_minima_cdk_8D():
    PB = SimpleNamespace()
    PB.name = _NAME
    PB.nvar = 8
    PB.var_types = ['C', 'K', 'D', 'C', 'D', 'K', 'C', 'K']
    PB.cardinality = [math.inf, 17.0, 17.0, math.inf, 17.0, 17.0, math.inf, 17.0]
    PB.f_star = 5.889542894199896e-32
    PB.f_star_justification = 'E'
    PB.n_minimizers_documented = 1
    PB.description = 'DFO Benchmark Suite v2.1.0 core (Python): CDK adaptation of tenn_local_minima (F = f o T).'
    PB.x0_idx = [-7.576920282840728, 10.0, 3.0, -152.00956980128183, 8.0, 10.0, -0.033612351453892825, 14.0]
    PB.x0_alt_idx = [[3.2892088294029254, 11.0, 14.0, -79.50224474357401, 3.0, 8.0, -0.08112575688424083, 3.0], [-5.363449537754058, 11.0, 6.0, -239.01835987053153, 15.0, 5.0, 0.03035994605363458, 8.0], [-0.9365080475807179, 10.0, 17.0, -166.5110348128236, 10.0, 3.0, -0.017153459376713448, 13.0], [-9.68977872133255, 11.0, 9.0, -94.00370975511555, 7.0, 17.0, -0.06466686480706146, 2.0], [6.810639560222626, 12.0, 1.0, -286.66603062273964, 2.0, 13.0, 0.05675628240382792, 8.0], [-1.8420188069343553, 10.0, 12.0, -214.1587055650317, 14.0, 11.0, 0.009242876973479905, 13.0], [2.584922683238984, 11.0, 4.0, -141.65138050732378, 9.0, 9.0, -0.04820797272988208, 1.0], [-6.218654143810272, 11.0, 15.0, -69.14405544961596, 5.0, 6.0, -0.09572137816023009, 6.0], [4.6474749684333805, 10.0, 7.0, -228.66017057657325, 16.0, 4.0, 0.02570176905065928, 12.0]]
    PB.suite_release = '2.1.0'
    PB.base_function = 'tenn_local_minima'
    PB.m = 17
    PB.pattern = 'CDK'
    PB.strategy = 'moderate'
    PB.scale_factors_C = [1.0, 100.0, 0.01]
    PB.lb_orig = [-10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0]
    PB.ub_orig = [7.1711669921875, 6.722013920642178, 7.2627328, 6.580008329862557, 6.216201079161259, 7.20804593676692, 6.952313789346391, 6.55068177807107]
    PB.lb_C_work = [-10.0, -10.0, -10.0, -1000.0, -10.0, -10.0, -0.09999999999999999, -10.0]
    PB.ub_C_work = [7.1711669921875005, 6.722013920642178, 7.2627328, 658.0008329862558, 6.216201079161259, 7.20804593676692, 0.06952313789346391, 6.55068177807107]
    PB.ProblemData = _problem_data()
    PB.func_idx = _objective
    PB.func_F = lambda xr: _objective(_cell2idx(xr))
    PB.idx2cell = _idx2cell
    PB.cell2idx = _cell2idx
    PB.decode = lambda xr: _decode(_cell2idx(xr))
    return PB
