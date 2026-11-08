"""salomon_c_10D  DFO Benchmark Suite v2.1.0 core, Python version (self-contained; standard library only).
SALOMON_C_10D  DFO Benchmark Suite v2.1.0 core: salomon, dimension 10, pattern C, strategy baseline.

 Copyright (C) 2026 J. F. A. Madeira

 Objective F(z) = f(T(z)), where f is the original function (embedded below, translated from the MATLAB original) and T acts on the finite
 coordinates only: uniform grid of 0 levels on the interval of the domain, piecewise-affine map with |T - Id| <=
 half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
 Global optimum of F on this domain: f* = 0 (E: exact by structure); the optimal solution is listed in
 reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
 Variables: x1:C x2:C x3:C x4:C x5:C x6:C x7:C x8:C x9:C x10:C


Python version of the MATLAB wrapper salomon_c_10D.m: same instance, same values (level indices 1..m as in MATLAB).
Usage:  from salomon_c_10D import salomon_c_10D;  PB = salomon_c_10D();  f = PB.func_idx(PB.x0_idx)
"""
# Copyright (C) 2026 J. F. A. Madeira. LGPL-3.0-or-later.
import math
from types import SimpleNamespace

_NAME = 'salomon_c_10D'
_N = 10
_TYPES = ['C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C']
_CARD = [math.inf, math.inf, math.inf, math.inf, math.inf, math.inf, math.inf, math.inf, math.inf, math.inf]
_LB = [-100.0, -100.0, -63.779328, -70.41232819658475, -69.91407690731508, -62.55383214873429, -100.0, -100.0, -100.0, -100.0]
_UB = [63.399169921875, 61.59325306101712, 100.0, 100.0, 100.0, 100.0, 69.5715808000383, 64.0719147336193, 70.75630804635489, 61.75997375866169]
_LBW = [-100.0, -100.0, -63.779328, -70.41232819658475, -69.91407690731508, -62.55383214873429, -100.0, -100.0, -100.0, -100.0]
_UBW = [63.399169921875, 61.59325306101712, 100.0, 100.0, 100.0, 100.0, 69.5715808000383, 64.0719147336193, 70.75630804635489, 61.75997375866169]
_TB = [math.nan, math.nan, math.nan, math.nan, math.nan, math.nan, math.nan, math.nan, math.nan, math.nan]
_TI = [math.nan, math.nan, math.nan, math.nan, math.nan, math.nan, math.nan, math.nan, math.nan, math.nan]
_CW = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.4210854715202004e-14, 0.0]
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
    s = 0.0
    for v in x: s += v ** 2
    r = math.sqrt(s)
    return 1 - math.cos(2 * math.pi * r) + 0.1 * r

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

def salomon_c_10D():
    PB = SimpleNamespace()
    PB.name = _NAME
    PB.nvar = 10
    PB.var_types = ['C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C']
    PB.cardinality = [math.inf, math.inf, math.inf, math.inf, math.inf, math.inf, math.inf, math.inf, math.inf, math.inf]
    PB.f_star = 0.0
    PB.f_star_justification = 'E'
    PB.n_minimizers_documented = 1
    PB.description = 'DFO Benchmark Suite v2.1.0 core (Python): C adaptation of salomon (F = f o T).'
    PB.x0_idx = [27.73538625240326, 53.90889169003137, 56.29057294336001, -63.243796273279855, 78.55329457518486, -22.67376135757918, -34.18013340409669, 54.12091363591031, -33.49107061463995, 57.34273063580841]
    PB.x0_alt_idx = [[-52.68764269351959, 58.56384136668618, -18.392800624639996, -55.79136209558664, 25.319507895733025, -42.87275825180059, -81.70711625808664, -64.93228940086225, -84.47797512128967, -75.95056522002102], [-12.476128220558167, 47.25896358052452, 89.30848546816, -48.338927917893436, -12.595203480279416, 98.52022000774933, 29.81044481761532, -11.756627284980098, 42.37879059178681, -47.676229735451145], [-92.89915716648102, 50.14059909464416, 14.62511190016, -64.7342831088185, -65.82899015973125, 65.89107117862243, -17.716538036374658, 41.41903483090206, -1.1979609890409932, -19.401894250881284], [53.58564555644989, 60.115491258904484, -36.15958212608, -57.28184893112528, 64.76612680208943, 45.692074284401016, -65.24352089036464, -69.0227249482378, -44.774712569868775, 8.872441233688605], [-26.83738338947296, 48.8106134727428, 46.38519918592002, -49.829414753432076, 10.139696913234019, 13.062925455274097, 56.214324180943095, -15.847062832355633, -88.35146415069659, 37.14677671825851], [13.374131083488464, 53.46556314939761, -3.1416696012800074, -69.6315969970169, -27.77501446277843, -7.136071438947319, 8.687341326953117, 37.3285992835265, 38.50530156237991, -96.14651913757093], [-67.04889786243439, 55.75609394267221, 79.40311171072003, -62.17916281932368, 87.50102719560287, -27.335068333168728, -48.779925522642614, -81.72460375324609, -5.071450018447905, -67.87218365300106], [33.001179814338684, 44.4512161565105, 29.876242923520003, -54.72672864163047, 49.58631581959041, -59.96421716229564, -96.30690837663259, -28.54894163736394, -48.6482015992757, -39.59784816843118], [-47.42184913158417, 54.426108320770794, -44.80713064448, -47.27429446393726, -3.6474708598614143, 81.42876109725431, 25.150936694675124, 24.62672047851821, -99.63510610592543, -11.323512683861296]]
    PB.suite_release = '2.1.0'
    PB.base_function = 'salomon'
    PB.m = 0
    PB.pattern = 'C'
    PB.strategy = 'baseline'
    PB.scale_factors_C = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
    PB.lb_orig = [-100.0, -100.0, -63.779328, -70.41232819658475, -69.91407690731508, -62.55383214873429, -100.0, -100.0, -100.0, -100.0]
    PB.ub_orig = [63.399169921875, 61.59325306101712, 100.0, 100.0, 100.0, 100.0, 69.5715808000383, 64.0719147336193, 70.75630804635489, 61.75997375866169]
    PB.lb_C_work = [-100.0, -100.0, -63.779328, -70.41232819658475, -69.91407690731508, -62.55383214873429, -100.0, -100.0, -100.0, -100.0]
    PB.ub_C_work = [63.399169921875, 61.59325306101712, 100.0, 100.0, 100.0, 100.0, 69.5715808000383, 64.0719147336193, 70.75630804635489, 61.75997375866169]
    PB.ProblemData = _problem_data()
    PB.func_idx = _objective
    PB.func_F = lambda xr: _objective(_cell2idx(xr))
    PB.idx2cell = _idx2cell
    PB.cell2idx = _cell2idx
    PB.decode = lambda xr: _decode(_cell2idx(xr))
    return PB
