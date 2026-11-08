"""schwefel_cd_10D  DFO Benchmark Suite v2.1.0 core, Python version (self-contained; standard library only).
SCHWEFEL_CD_10D  DFO Benchmark Suite v2.1.0 core: schwefel, dimension 10, pattern CD, strategy extreme.

 Copyright (C) 2026 J. F. A. Madeira

 Objective F(z) = f(T(z)), where f is the original function (embedded below, translated from the MATLAB original) and T acts on the finite
 coordinates only: uniform grid of 17 levels on the interval of the domain, piecewise-affine map with |T - Id| <=
 half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
 Global optimum of F on this domain: f* = -4189.828872721624 (A: analytic reduction, 1-D root to machine precision); the optimal solution is listed in
 reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
 Variables: x1:C x2:D x3:C x4:D x5:C x6:D x7:C x8:D x9:C x10:D


Python version of the MATLAB wrapper schwefel_cd_10D.m: same instance, same values (level indices 1..m as in MATLAB).
Usage:  from schwefel_cd_10D import schwefel_cd_10D;  PB = schwefel_cd_10D();  f = PB.func_idx(PB.x0_idx)
"""
# Copyright (C) 2026 J. F. A. Madeira. LGPL-3.0-or-later.
import math
from types import SimpleNamespace

_NAME = 'schwefel_cd_10D'
_N = 10
_TYPES = ['C', 'D', 'C', 'D', 'C', 'D', 'C', 'D', 'C', 'D']
_CARD = [math.inf, 17.0, math.inf, 17.0, math.inf, 17.0, math.inf, 17.0, math.inf, 17.0]
_LB = [-343.245849609375, -300.1884875273078, -344.09664, -322.0616409829238, -330.47947544566625, -366.20111340639335, -335.504962823721, -309.3069420891491, -344.6511054491658, -301.55848948296364]
_UB = [500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0]
_LBW = [-343.245849609375, -300.1884875273078, -344.09664, -322.0616409829238, -33047947544.566624, -366.20111340639335, -33550496282.372093, -309.3069420891491, -34465110544.91658, -301.55848948296364]
_UBW = [500.0, 500.0, 500.0, 500.0, 50000000000.0, 500.0, 49999999999.99999, 500.0, 50000000000.0, 500.0]
_TB = [math.nan, 399.97643905908654, math.nan, 397.24229487713455, math.nan, 445.86243041210037, math.nan, 398.83663223885634, math.nan, 399.80518881462956]
_TI = [math.nan, 420.9687, math.nan, 420.9687, math.nan, 420.9687, math.nan, 420.9687, math.nan, 420.9687]
_CW = [420.9687, math.nan, 420.9687, math.nan, 42096870000.0, math.nan, 42096870000.00001, math.nan, 42096870000.0, math.nan]
_CX = [420.9687, math.nan, 420.9687, math.nan, 420.9687, math.nan, 420.9687, math.nan, 420.9687, math.nan]
_M = 17
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

def _orig(x):
    s = 0.0
    for v in x: s += v * math.sin(math.sqrt(abs(v)))
    return -s

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

def schwefel_cd_10D():
    PB = SimpleNamespace()
    PB.name = _NAME
    PB.nvar = 10
    PB.var_types = ['C', 'D', 'C', 'D', 'C', 'D', 'C', 'D', 'C', 'D']
    PB.cardinality = [math.inf, 17.0, math.inf, 17.0, math.inf, 17.0, math.inf, 17.0, math.inf, 17.0]
    PB.f_star = -4189.828872721624
    PB.f_star_justification = 'A'
    PB.n_minimizers_documented = 1
    PB.description = 'DFO Benchmark Suite v2.1.0 core (Python): CD adaptation of schwefel (F = f o T).'
    PB.x0_idx = [-0.2654820680618286, 13.0, -29.1473016832, 9.0, -20194546196.57209, 9.0, 13624158040.302475, 3.0, 9450915534.408012, 2.0]
    PB.x0_alt_idx = [[217.13383853435516, 12.0, 396.2774048768001, 9.0, 36148276217.20978, 7.0, -9793109083.82257, 7.0, -15769906904.129173, 5.0], [-207.78301537036896, 13.0, 141.02258094080003, 10.0, 17616916021.31474, 5.0, 45153288939.45124, 12.0, 46980165769.126434, 8.0], [325.83349883556366, 13.0, -243.88548689919998, 8.0, -8401862435.547985, 1.0, 21736021815.326195, 1.0, 25424797293.580605, 11.0], [-104.0242487192154, 12.0, 304.7098013696001, 9.0, -26933222631.44301, 15.0, -1681245308.7988663, 7.0, 3869428818.0347977, 13.0], [113.3750718832016, 13.0, -80.19826647039997, 10.0, 28728925394.978146, 12.0, -25098512432.923912, 12.0, -17685939657.511024, 16.0], [-311.5417820215225, 13.0, 474.8796839936001, 8.0, 10197565199.083118, 10.0, 34745614662.062325, 17.0, 45064133015.74458, 3.0], [459.2376273870468, 12.0, 89.97161615360005, 9.0, -15821213257.779598, 8.0, 6430618466.224846, 5.0, 23508764540.19876, 6.0], [26.909433007240295, 12.0, -327.3497626624, 9.0, 48009027416.96997, 5.0, -16986648657.900204, 10.0, 1953396064.6529388, 9.0], [244.3087536096573, 13.0, 260.1414987776001, 10.0, 21990248960.107227, 3.0, 42857478437.08604, 16.0, -23267426373.884247, 11.0]]
    PB.suite_release = '2.1.0'
    PB.base_function = 'schwefel'
    PB.m = 17
    PB.pattern = 'CD'
    PB.strategy = 'extreme'
    PB.scale_factors_C = [1.0, 1.0, 100000000.0, 100000000.0, 100000000.0]
    PB.lb_orig = [-343.245849609375, -300.1884875273078, -344.09664, -322.0616409829238, -330.47947544566625, -366.20111340639335, -335.504962823721, -309.3069420891491, -344.6511054491658, -301.55848948296364]
    PB.ub_orig = [500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0]
    PB.lb_C_work = [-343.245849609375, -300.1884875273078, -344.09664, -322.0616409829238, -33047947544.566624, -366.20111340639335, -33550496282.372093, -309.3069420891491, -34465110544.91658, -301.55848948296364]
    PB.ub_C_work = [500.0, 500.0, 500.0, 500.0, 50000000000.0, 500.0, 49999999999.99999, 500.0, 50000000000.0, 500.0]
    PB.ProblemData = _problem_data()
    PB.func_idx = _objective
    PB.func_F = lambda xr: _objective(_cell2idx(xr))
    PB.idx2cell = _idx2cell
    PB.cell2idx = _cell2idx
    PB.decode = lambda xr: _decode(_cell2idx(xr))
    return PB
