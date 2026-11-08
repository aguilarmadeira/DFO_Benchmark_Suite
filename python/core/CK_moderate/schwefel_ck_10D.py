"""schwefel_ck_10D  DFO Benchmark Suite v2.1.0 core, Python version (self-contained; standard library only).
SCHWEFEL_CK_10D  DFO Benchmark Suite v2.1.0 core: schwefel, dimension 10, pattern CK, strategy moderate.

 Copyright (C) 2026 J. F. A. Madeira

 Objective F(z) = f(T(z)), where f is the original function (embedded below, translated from the MATLAB original) and T acts on the finite
 coordinates only: uniform grid of 17 levels on the interval of the domain, piecewise-affine map with |T - Id| <=
 half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
 Global optimum of F on this domain: f* = -4189.828872721624 (A: analytic reduction, 1-D root to machine precision); the optimal solution is listed in
 reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
 Variables: x1:C x2:K x3:C x4:K x5:C x6:K x7:C x8:K x9:C x10:K


Python version of the MATLAB wrapper schwefel_ck_10D.m: same instance, same values (level indices 1..m as in MATLAB).
Usage:  from schwefel_ck_10D import schwefel_ck_10D;  PB = schwefel_ck_10D();  f = PB.func_idx(PB.x0_idx)
"""
# Copyright (C) 2026 J. F. A. Madeira. LGPL-3.0-or-later.
import math
from types import SimpleNamespace

_NAME = 'schwefel_ck_10D'
_N = 10
_TYPES = ['C', 'K', 'C', 'K', 'C', 'K', 'C', 'K', 'C', 'K']
_CARD = [math.inf, 17.0, math.inf, 17.0, math.inf, 17.0, math.inf, 17.0, math.inf, 17.0]
_LB = [-343.245849609375, -300.1884875273078, -344.09664, -322.0616409829238, -330.47947544566625, -366.20111340639335, -335.504962823721, -309.3069420891491, -344.6511054491658, -301.55848948296364]
_UB = [500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0]
_LBW = [-343.245849609375, -300.1884875273078, -34409.664000000004, -322.0616409829238, -3.3047947544566623, -366.20111340639335, -335.5049628237209, -309.3069420891491, -34465.11054491658, -301.55848948296364]
_UBW = [500.0, 500.0, 50000.0, 500.0, 5.0, 500.0, 500.0, 500.0, 50000.0, 500.0]
_TB = [math.nan, 399.97643905908654, math.nan, 397.24229487713455, math.nan, 445.86243041210037, math.nan, 398.83663223885634, math.nan, 399.80518881462956]
_TI = [math.nan, 420.9687, math.nan, 420.9687, math.nan, 420.9687, math.nan, 420.9687, math.nan, 420.9687]
_CW = [420.9687, math.nan, 42096.869999999995, math.nan, 4.209687000000001, math.nan, 420.9687000000001, math.nan, 42096.869999999995, math.nan]
_CX = [420.9687, math.nan, 420.9687, math.nan, 420.9687, math.nan, 420.9687, math.nan, 420.9687, math.nan]
_M = 17
_PERM = {1: [8, 5, 16, 6, 10, 13, 12, 1, 17, 2, 11, 7, 3, 15, 14, 9, 4], 3: [2, 16, 9, 14, 8, 6, 17, 7, 10, 12, 3, 13, 15, 1, 4, 11, 5], 5: [1, 17, 11, 8, 7, 13, 9, 15, 3, 16, 6, 12, 14, 4, 10, 5, 2], 7: [11, 12, 13, 14, 5, 10, 15, 3, 2, 16, 9, 17, 4, 6, 1, 7, 8], 9: [8, 4, 10, 16, 13, 12, 6, 11, 2, 15, 1, 5, 17, 7, 9, 14, 3]}

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

def schwefel_ck_10D():
    PB = SimpleNamespace()
    PB.name = _NAME
    PB.nvar = 10
    PB.var_types = ['C', 'K', 'C', 'K', 'C', 'K', 'C', 'K', 'C', 'K']
    PB.cardinality = [math.inf, 17.0, math.inf, 17.0, math.inf, 17.0, math.inf, 17.0, math.inf, 17.0]
    PB.f_star = -4189.828872721624
    PB.f_star_justification = 'A'
    PB.n_minimizers_documented = 1
    PB.description = 'DFO Benchmark Suite v2.1.0 core (Python): CK adaptation of schwefel (F = f o T).'
    PB.x0_idx = [-0.2654820680618286, 13.0, -2914.730168320002, 9.0, -2.0194546196572087, 9.0, 136.24158040302484, 3.0, 9450.915534408014, 2.0]
    PB.x0_alt_idx = [[217.13383853435516, 12.0, 39627.740487680014, 9.0, 3.6148276217209774, 7.0, -97.93109083822569, 7.0, -15769.906904129173, 5.0], [-207.78301537036896, 13.0, 14102.258094079996, 10.0, 1.761691602131474, 5.0, 451.5328893945125, 12.0, 46980.165769126426, 8.0], [325.83349883556366, 13.0, -24388.548689920004, 8.0, -0.8401862435547982, 1.0, 217.36021815326194, 1.0, 25424.79729358061, 11.0], [-104.0242487192154, 12.0, 30470.98013696001, 9.0, -2.693322263144301, 15.0, -16.812453087988615, 7.0, 3869.428818034794, 13.0], [113.3750718832016, 13.0, -8019.826647040001, 10.0, 2.872892539497815, 12.0, -250.98512432923908, 12.0, -17685.939657511022, 16.0], [-311.5417820215225, 13.0, 47487.968399360005, 8.0, 1.0197565199083125, 10.0, 347.4561466206234, 17.0, 45064.13301574458, 3.0], [459.2376273870468, 12.0, 8997.161615359997, 9.0, -1.5821213257779598, 8.0, 64.30618466224848, 5.0, 23508.76454019876, 6.0], [26.909433007240295, 12.0, -32734.976266240003, 9.0, 4.800902741696996, 5.0, -169.86648657900201, 10.0, 1953.3960646529376, 9.0], [244.3087536096573, 13.0, 26014.14987776, 10.0, 2.1990248960107226, 3.0, 428.57478437086047, 16.0, -23267.426373884242, 11.0]]
    PB.suite_release = '2.1.0'
    PB.base_function = 'schwefel'
    PB.m = 17
    PB.pattern = 'CK'
    PB.strategy = 'moderate'
    PB.scale_factors_C = [1.0, 100.0, 0.01, 1.0, 100.0]
    PB.lb_orig = [-343.245849609375, -300.1884875273078, -344.09664, -322.0616409829238, -330.47947544566625, -366.20111340639335, -335.504962823721, -309.3069420891491, -344.6511054491658, -301.55848948296364]
    PB.ub_orig = [500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0]
    PB.lb_C_work = [-343.245849609375, -300.1884875273078, -34409.664000000004, -322.0616409829238, -3.3047947544566623, -366.20111340639335, -335.5049628237209, -309.3069420891491, -34465.11054491658, -301.55848948296364]
    PB.ub_C_work = [500.0, 500.0, 50000.0, 500.0, 5.0, 500.0, 500.0, 500.0, 50000.0, 500.0]
    PB.ProblemData = _problem_data()
    PB.func_idx = _objective
    PB.func_F = lambda xr: _objective(_cell2idx(xr))
    PB.idx2cell = _idx2cell
    PB.cell2idx = _cell2idx
    PB.decode = lambda xr: _decode(_cell2idx(xr))
    return PB
