"""ackley_cdk_10D  DFO Benchmark Suite v2.1.0 core, Python version (self-contained; standard library only).
ACKLEY_CDK_10D  DFO Benchmark Suite v2.1.0 core: ackley, dimension 10, pattern CDK, strategy baseline.

 Copyright (C) 2026 J. F. A. Madeira

 Objective F(z) = f(T(z)), where f is the original function (embedded below, translated from the MATLAB original) and T acts on the finite
 coordinates only: uniform grid of 17 levels on the interval of the domain, piecewise-affine map with |T - Id| <=
 half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
 Global optimum of F on this domain: f* = 4.440892098500626e-16 (E: exact by structure); the optimal solution is listed in
 reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
 Variables: x1:C x2:K x3:D x4:C x5:D x6:K x7:C x8:K x9:D x10:C


Python version of the MATLAB wrapper ackley_cdk_10D.m: same instance, same values (level indices 1..m as in MATLAB).
Usage:  from ackley_cdk_10D import ackley_cdk_10D;  PB = ackley_cdk_10D();  f = PB.func_idx(PB.x0_idx)
"""
# Copyright (C) 2026 J. F. A. Madeira. LGPL-3.0-or-later.
import math
from types import SimpleNamespace

_NAME = 'ackley_cdk_10D'
_N = 10
_TYPES = ['C', 'K', 'D', 'C', 'D', 'K', 'C', 'K', 'D', 'C']
_CARD = [math.inf, 17.0, 17.0, math.inf, 17.0, 17.0, math.inf, 17.0, 17.0, math.inf]
_LB = [-19.659594726562503, -30.0, -21.3513984, -21.050229071220325, -19.970767024110373, -30.0, -18.69154344416374, -20.12905364446252, -20.512336648311006, -20.415744802984953]
_UB = [30.0, 18.374272214601433, 30.0, 30.0, 30.0, 19.809936626868804, 30.0, 30.0, 30.0, 30.0]
_LBW = [-19.659594726562503, -30.0, -21.3513984, -21.050229071220325, -19.970767024110373, -30.0, -18.69154344416374, -20.12905364446252, -20.512336648311006, -20.415744802984953]
_UBW = [30.0, 18.374272214601433, 30.0, 30.0, 30.0, 19.809936626868804, 30.0, 30.0, 30.0, 30.0]
_TB = [math.nan, 0.23392013412589563, 1.1148384, math.nan, -1.2317293900689847, 1.1312103917930045, math.nan, -1.3306585277890761, -1.5702104051943806, math.nan]
_TI = [math.nan, 0.0, 0.0, math.nan, 0.0, 0.0, math.nan, 0.0, 0.0, math.nan]
_CW = [0.0, math.nan, math.nan, 0.0, math.nan, math.nan, 0.0, math.nan, math.nan, 0.0]
_CX = [0.0, math.nan, math.nan, 0.0, math.nan, math.nan, 0.0, math.nan, math.nan, 0.0]
_M = 17
_PERM = {1: [4, 11, 9, 7, 8, 3, 1, 12, 13, 5, 2, 16, 14, 6, 15, 10, 17], 5: [2, 13, 12, 10, 17, 7, 16, 1, 5, 8, 14, 3, 6, 9, 11, 4, 15], 7: [4, 8, 13, 17, 3, 14, 15, 7, 5, 10, 1, 11, 16, 12, 2, 6, 9]}

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
    s1 = 0.0
    for v in x: s1 += v ** 2
    s2 = 0.0
    for v in x: s2 += math.cos(2 * math.pi * v)
    return -20 * math.exp(-0.02 * math.sqrt(1 / n * s1)) - math.exp(1 / n * s2) + 20 + math.exp(1)

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

def ackley_cdk_10D():
    PB = SimpleNamespace()
    PB.name = _NAME
    PB.nvar = 10
    PB.var_types = ['C', 'K', 'D', 'C', 'D', 'K', 'C', 'K', 'D', 'C']
    PB.cardinality = [math.inf, 17.0, 17.0, math.inf, 17.0, 17.0, math.inf, 17.0, 17.0, math.inf]
    PB.f_star = 4.440892098500626e-16
    PB.f_star_justification = 'E'
    PB.n_minimizers_documented = 1
    PB.description = 'DFO Benchmark Suite v2.1.0 core (Python): CDK adaptation of ackley (F = f o T).'
    PB.x0_idx = [9.777293944358824, 12.0, 4.0, -16.989173780673436, 1.0, 2.0, 25.48069533675175, 13.0, 10.0, 5.963168618266064]
    PB.x0_alt_idx = [[-15.755691039562228, 12.0, 15.0, -15.777235814384573, 14.0, 15.0, 11.833584544442882, 2.0, 5.0, 14.7754331676582], [24.98069525957108, 13.0, 8.0, -20.68877388829207, 9.0, 13.0, -1.8135262478659904, 6.0, 1.0, 23.58769771705034], [-0.04308489561081075, 12.0, 3.0, -18.45625658197048, 5.0, 10.0, -18.314934722226386, 11.0, 13.0, -17.955835158655457], [12.759779369831087, 13.0, 11.0, -16.22373927564889, 16.0, 8.0, 16.561015080340724, 17.0, 9.0, -9.14357060926332], [-12.554974973201755, 13.0, 6.0, -13.991221969327302, 12.0, 6.0, 2.9139042880318513, 6.0, 5.0, -0.33130605987118145], [18.87023731470108, 12.0, 14.0, -19.923339383267525, 7.0, 2.0, -10.733206504277016, 11.0, 17.0, 8.480958489520958], [-6.153542840480807, 13.0, 10.0, -17.690822076945935, 4.0, 16.0, 24.142743298290092, 16.0, 13.0, 17.293223038913098], [6.649321424961087, 13.0, 2.0, -15.458304770624347, 15.0, 14.0, 7.6413348239297, 5.0, 9.0, 26.105487588305238], [-18.810920011997226, 12.0, 13.0, -20.369842844531842, 11.0, 11.0, -6.005775968379174, 10.0, 4.0, -15.438045287400561]]
    PB.suite_release = '2.1.0'
    PB.base_function = 'ackley'
    PB.m = 17
    PB.pattern = 'CDK'
    PB.strategy = 'baseline'
    PB.scale_factors_C = [1.0, 1.0, 1.0, 1.0]
    PB.lb_orig = [-19.659594726562503, -30.0, -21.3513984, -21.050229071220325, -19.970767024110373, -30.0, -18.69154344416374, -20.12905364446252, -20.512336648311006, -20.415744802984953]
    PB.ub_orig = [30.0, 18.374272214601433, 30.0, 30.0, 30.0, 19.809936626868804, 30.0, 30.0, 30.0, 30.0]
    PB.lb_C_work = [-19.659594726562503, -30.0, -21.3513984, -21.050229071220325, -19.970767024110373, -30.0, -18.69154344416374, -20.12905364446252, -20.512336648311006, -20.415744802984953]
    PB.ub_C_work = [30.0, 18.374272214601433, 30.0, 30.0, 30.0, 19.809936626868804, 30.0, 30.0, 30.0, 30.0]
    PB.ProblemData = _problem_data()
    PB.func_idx = _objective
    PB.func_F = lambda xr: _objective(_cell2idx(xr))
    PB.idx2cell = _idx2cell
    PB.cell2idx = _cell2idx
    PB.decode = lambda xr: _decode(_cell2idx(xr))
    return PB
