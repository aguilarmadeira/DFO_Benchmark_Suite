"""ackley_c_10D  DFO Benchmark Suite v2.1.0 core, Python version (self-contained; standard library only).
ACKLEY_C_10D  DFO Benchmark Suite v2.1.0 core: ackley, dimension 10, pattern C, strategy moderate.

 Copyright (C) 2026 J. F. A. Madeira

 Objective F(z) = f(T(z)), where f is the original function (embedded below, translated from the MATLAB original) and T acts on the finite
 coordinates only: uniform grid of 0 levels on the interval of the domain, piecewise-affine map with |T - Id| <=
 half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
 Global optimum of F on this domain: f* = 4.440892098500626e-16 (E: exact by structure); the optimal solution is listed in
 reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
 Variables: x1:C x2:C x3:C x4:C x5:C x6:C x7:C x8:C x9:C x10:C


Python version of the MATLAB wrapper ackley_c_10D.m: same instance, same values (level indices 1..m as in MATLAB).
Usage:  from ackley_c_10D import ackley_c_10D;  PB = ackley_c_10D();  f = PB.func_idx(PB.x0_idx)
"""
# Copyright (C) 2026 J. F. A. Madeira. LGPL-3.0-or-later.
import math
from types import SimpleNamespace

_NAME = 'ackley_c_10D'
_N = 10
_TYPES = ['C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C']
_CARD = [math.inf, math.inf, math.inf, math.inf, math.inf, math.inf, math.inf, math.inf, math.inf, math.inf]
_LB = [-19.659594726562503, -30.0, -21.3513984, -21.050229071220325, -19.970767024110373, -30.0, -18.69154344416374, -20.12905364446252, -20.512336648311006, -20.415744802984953]
_UB = [30.0, 18.374272214601433, 30.0, 30.0, 30.0, 19.809936626868804, 30.0, 30.0, 30.0, 30.0]
_LBW = [-19.659594726562503, -3000.0, -0.21351398400000002, -21.05022907122033, -1997.0767024110373, -0.3, -18.69154344416374, -2012.905364446252, -0.20512336648311008, -20.415744802984953]
_UBW = [30.0, 1837.4272214601433, 0.30000000000000004, 30.0, 3000.0, 0.19809936626868804, 30.0, 3000.0, 0.3, 30.0]
_TB = [math.nan, math.nan, math.nan, math.nan, math.nan, math.nan, math.nan, math.nan, math.nan, math.nan]
_TI = [math.nan, math.nan, math.nan, math.nan, math.nan, math.nan, math.nan, math.nan, math.nan, math.nan]
_CW = [0.0, 4.547473508864641e-13, 2.7755575615628914e-17, -3.552713678800501e-15, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
_CX = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
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

def ackley_c_10D():
    PB = SimpleNamespace()
    PB.name = _NAME
    PB.nvar = 10
    PB.var_types = ['C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C']
    PB.cardinality = [math.inf, math.inf, math.inf, math.inf, math.inf, math.inf, math.inf, math.inf, math.inf, math.inf]
    PB.f_star = 4.440892098500626e-16
    PB.f_star_justification = 'E'
    PB.n_minimizers_documented = 1
    PB.description = 'DFO Benchmark Suite v2.1.0 core (Python): C adaptation of ackley (F = f o T).'
    PB.x0_idx = [9.777293944358824, 293.5204082094888, -0.10818199560192002, -16.98917378067344, -1843.147217212801, -0.2634984078428499, 25.48069533675175, 1612.8452570755226, 0.07755806709817509, 5.963168618266064]
    PB.x0_alt_idx = [[-15.755691039562228, 432.87016356019285, 0.24922373726208005, -15.777235814384577, 1997.5811573675833, 0.13167100290194905, 11.833584544442882, -1761.4927029256935, -0.0732690217843053, 14.7754331676582], [24.98069525957108, 731.4767821688438, -0.004657576427519988, -20.688773888292072, 432.00566157088724, 0.06977699880939014, -1.8135262478659904, -399.91773660941317, -0.20217574101912547, 23.58769771705034], [-0.04308489561081075, 393.0559477457059, -0.14022526820352002, -18.456256581970482, -683.0445116943858, -0.030205623186281783, -18.314934722226386, 1224.763503391173, 0.17308604164223979, -17.955835158655457], [12.759779369831087, 479.32008201042845, 0.09492305534976006, -16.223739275648896, 2707.158540354576, -0.09209962727884072, 16.561015080340724, 2849.4447433917585, 0.0441793224074197, -9.14357060926332], [-12.554974973201755, 618.6698373611312, -0.04064463642624, -13.991221969327306, 1592.1083670893022, -0.15399363137139957, 2.9139042880318513, -524.8932166094582, -0.08472739682740044, -0.33130605987118145], [18.87023731470108, 280.24900293799374, 0.19844747452416003, -19.92333938326753, -14.423976210250885, -0.2539762533670716, -10.733206504277016, 1099.7880233911283, 0.2905343858339649, 8.480958489520958], [-6.153542840480807, 578.8556215466447, 0.06287978274816, -17.69082207694594, -1129.4741494755235, 0.17928177528084044, 24.142743298290092, 2724.469263391713, 0.16162766659914468, 17.293223038913098], [6.649321424961087, 647.4245487827056, -0.19100153094144, -15.45830477062435, 2260.728902573437, 0.11738777118828159, 7.6413348239297, -912.9749702938086, 0.0327209473643246, 26.105487588305238], [-18.810920011997226, 309.0037143595673, 0.16640420192256003, -20.369842844531846, 1145.6787293081638, 0.017405149192609604, -6.005775968379174, 711.7062697067768, -0.09618577187049555, -15.438045287400561]]
    PB.suite_release = '2.1.0'
    PB.base_function = 'ackley'
    PB.m = 0
    PB.pattern = 'C'
    PB.strategy = 'moderate'
    PB.scale_factors_C = [1.0, 100.0, 0.01, 1.0, 100.0, 0.01, 1.0, 100.0, 0.01, 1.0]
    PB.lb_orig = [-19.659594726562503, -30.0, -21.3513984, -21.050229071220325, -19.970767024110373, -30.0, -18.69154344416374, -20.12905364446252, -20.512336648311006, -20.415744802984953]
    PB.ub_orig = [30.0, 18.374272214601433, 30.0, 30.0, 30.0, 19.809936626868804, 30.0, 30.0, 30.0, 30.0]
    PB.lb_C_work = [-19.659594726562503, -3000.0, -0.21351398400000002, -21.05022907122033, -1997.0767024110373, -0.3, -18.69154344416374, -2012.905364446252, -0.20512336648311008, -20.415744802984953]
    PB.ub_C_work = [30.0, 1837.4272214601433, 0.30000000000000004, 30.0, 3000.0, 0.19809936626868804, 30.0, 3000.0, 0.3, 30.0]
    PB.ProblemData = _problem_data()
    PB.func_idx = _objective
    PB.func_F = lambda xr: _objective(_cell2idx(xr))
    PB.idx2cell = _idx2cell
    PB.cell2idx = _cell2idx
    PB.decode = lambda xr: _decode(_cell2idx(xr))
    return PB
