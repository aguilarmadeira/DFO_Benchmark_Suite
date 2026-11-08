"""schwefel_c_10D  DFO Benchmark Suite v2.1.0 core, Python version (self-contained; standard library only).
SCHWEFEL_C_10D  DFO Benchmark Suite v2.1.0 core: schwefel, dimension 10, pattern C, strategy extreme.

 Copyright (C) 2026 J. F. A. Madeira

 Objective F(z) = f(T(z)), where f is the original function (embedded below, translated from the MATLAB original) and T acts on the finite
 coordinates only: uniform grid of 0 levels on the interval of the domain, piecewise-affine map with |T - Id| <=
 half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
 Global optimum of F on this domain: f* = -4189.828872721624 (A: analytic reduction, 1-D root to machine precision); the optimal solution is listed in
 reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
 Variables: x1:C x2:C x3:C x4:C x5:C x6:C x7:C x8:C x9:C x10:C


Python version of the MATLAB wrapper schwefel_c_10D.m: same instance, same values (level indices 1..m as in MATLAB).
Usage:  from schwefel_c_10D import schwefel_c_10D;  PB = schwefel_c_10D();  f = PB.func_idx(PB.x0_idx)
"""
# Copyright (C) 2026 J. F. A. Madeira. LGPL-3.0-or-later.
import math
from types import SimpleNamespace

_NAME = 'schwefel_c_10D'
_N = 10
_TYPES = ['C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C']
_CARD = [math.inf, math.inf, math.inf, math.inf, math.inf, math.inf, math.inf, math.inf, math.inf, math.inf]
_LB = [-343.245849609375, -300.1884875273078, -344.09664, -322.0616409829238, -330.47947544566625, -366.20111340639335, -335.504962823721, -309.3069420891491, -344.6511054491658, -301.55848948296364]
_UB = [500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0]
_LBW = [-343.245849609375, -300.18848752730787, -344.09664, -322.0616409829238, -330.47947544566625, -36620111340.63933, -33550496282.372093, -30930694208.91491, -34465110544.91658, -30155848948.29636]
_UBW = [500.0, 500.0, 500.0, 500.0, 500.0, 50000000000.0, 49999999999.99999, 49999999999.99999, 50000000000.0, 50000000000.0]
_TB = [math.nan, math.nan, math.nan, math.nan, math.nan, math.nan, math.nan, math.nan, math.nan, math.nan]
_TI = [math.nan, math.nan, math.nan, math.nan, math.nan, math.nan, math.nan, math.nan, math.nan, math.nan]
_CW = [420.9687, 420.9686999999999, 420.9687, 420.9687, 420.9687, 42096870000.0, 42096870000.00001, 42096870000.00001, 42096870000.0, 42096870000.0]
_CX = [420.9687, 420.9687, 420.9687, 420.9687, 420.9687, 420.9687, 420.9687, 420.9687, 420.9687, 420.9687]
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

def schwefel_c_10D():
    PB = SimpleNamespace()
    PB.name = _NAME
    PB.nvar = 10
    PB.var_types = ['C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C']
    PB.cardinality = [math.inf, math.inf, math.inf, math.inf, math.inf, math.inf, math.inf, math.inf, math.inf, math.inf]
    PB.f_star = -4189.828872721624
    PB.f_star_justification = 'A'
    PB.n_minimizers_documented = 1
    PB.description = 'DFO Benchmark Suite v2.1.0 core (Python): C adaptation of schwefel (F = f o T).'
    PB.x0_idx = [-0.2654820680618286, 313.76500221700974, -29.1473016832, 64.83117630599907, -201.9454619657209, 4620050908.932243, 13624158040.302475, -22836444868.294464, 9450915534.408012, -23957401651.708267]
    PB.x0_alt_idx = [[217.13383853435516, 257.78473765748606, 396.2774048768001, 100.78139384169555, 361.48276217209775, -6143394878.957855, -9793109083.82257, -854540317.8923035, -15769906904.129173, -9946807364.669552], [-207.78301537036896, 280.83543482905463, 141.02258094080003, 136.73161137739186, 176.1691602131474, -16906840666.847946, 45153288939.45124, 25375075256.742165, 46980165769.126434, 4063786922.369171], [325.83349883556366, 295.10491403050185, -243.88548689919998, 57.641132798859815, -84.01862435547983, -34293945401.131943, 21736021815.326195, -29101818629.037113, 25424797293.580605, 18074381209.407894], [-104.0242487192154, 239.12464947097828, 304.7098013696001, 93.5913503345563, -269.3322263144301, 41050175114.098724, -1681245308.7988663, -2872203054.402645, 3869428818.0347977, 32084975496.446613], [113.3750718832016, 288.51900055291094, -80.19826647039997, 129.5415678702526, 287.2892539497815, 23663070379.814713, -25098512432.923912, 23357412520.231827, -17685939657.511024, 46095569783.48534], [-311.5417820215225, 311.5696977244794, 474.8796839936001, 34.016704132545044, 101.97565199083124, 12899624591.92463, 34745614662.062325, 49587028094.866295, 45064133015.74458, -19954374712.55435], [459.2376273870468, 243.88114253812728, 89.97161615360005, 69.96692166824141, -158.21213257779598, 2136178804.0345306, 6430618466.224846, -9137576815.14529, 23508764540.19876, -5943780425.515629], [26.909433007240295, 266.93183970969596, -327.3497626624, 105.91713920393778, 480.0902741696997, -15250925930.24947, -16986648657.900204, 17092038759.489182, 1953396064.6529388, 8066813861.52309], [244.3087536096573, 316.3261907916285, 260.1414987776001, 141.8673567396341, 219.90248960107226, -26014371718.139565, 42857478437.08604, 43321654334.123634, -23267426373.884247, 22077408148.561817]]
    PB.suite_release = '2.1.0'
    PB.base_function = 'schwefel'
    PB.m = 0
    PB.pattern = 'C'
    PB.strategy = 'extreme'
    PB.scale_factors_C = [1.0, 1.0, 1.0, 1.0, 1.0, 100000000.0, 100000000.0, 100000000.0, 100000000.0, 100000000.0]
    PB.lb_orig = [-343.245849609375, -300.1884875273078, -344.09664, -322.0616409829238, -330.47947544566625, -366.20111340639335, -335.504962823721, -309.3069420891491, -344.6511054491658, -301.55848948296364]
    PB.ub_orig = [500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0]
    PB.lb_C_work = [-343.245849609375, -300.18848752730787, -344.09664, -322.0616409829238, -330.47947544566625, -36620111340.63933, -33550496282.372093, -30930694208.91491, -34465110544.91658, -30155848948.29636]
    PB.ub_C_work = [500.0, 500.0, 500.0, 500.0, 500.0, 50000000000.0, 49999999999.99999, 49999999999.99999, 50000000000.0, 50000000000.0]
    PB.ProblemData = _problem_data()
    PB.func_idx = _objective
    PB.func_F = lambda xr: _objective(_cell2idx(xr))
    PB.idx2cell = _idx2cell
    PB.cell2idx = _cell2idx
    PB.decode = lambda xr: _decode(_cell2idx(xr))
    return PB
