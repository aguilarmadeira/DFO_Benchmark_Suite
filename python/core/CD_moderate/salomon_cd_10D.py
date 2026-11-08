"""salomon_cd_10D  DFO Benchmark Suite v2.1.0 core, Python version (self-contained; standard library only).
SALOMON_CD_10D  DFO Benchmark Suite v2.1.0 core: salomon, dimension 10, pattern CD, strategy moderate.

 Copyright (C) 2026 J. F. A. Madeira

 Objective F(z) = f(T(z)), where f is the original function (embedded below, translated from the MATLAB original) and T acts on the finite
 coordinates only: uniform grid of 17 levels on the interval of the domain, piecewise-affine map with |T - Id| <=
 half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
 Global optimum of F on this domain: f* = 0 (E: exact by structure); the optimal solution is listed in
 reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
 Variables: x1:C x2:D x3:C x4:D x5:C x6:D x7:C x8:D x9:C x10:D


Python version of the MATLAB wrapper salomon_cd_10D.m: same instance, same values (level indices 1..m as in MATLAB).
Usage:  from salomon_cd_10D import salomon_cd_10D;  PB = salomon_cd_10D();  f = PB.func_idx(PB.x0_idx)
"""
# Copyright (C) 2026 J. F. A. Madeira. LGPL-3.0-or-later.
import math
from types import SimpleNamespace

_NAME = 'salomon_cd_10D'
_N = 10
_TYPES = ['C', 'D', 'C', 'D', 'C', 'D', 'C', 'D', 'C', 'D']
_CARD = [math.inf, 17.0, math.inf, 17.0, math.inf, 17.0, math.inf, 17.0, math.inf, 17.0]
_LB = [-100.0, -100.0, -63.779328, -70.41232819658475, -69.91407690731508, -62.55383214873429, -100.0, -100.0, -100.0, -100.0]
_UB = [63.399169921875, 61.59325306101712, 100.0, 100.0, 100.0, 100.0, 69.5715808000383, 64.0719147336193, 70.75630804635489, 61.75997375866169]
_LBW = [-100.0, -100.0, -6377.9328, -70.41232819658475, -0.6991407690731507, -62.55383214873429, -100.0, -100.0, -10000.0, -100.0]
_UBW = [63.399169921875, 61.59325306101712, 10000.0, 100.0, 0.9999999999999999, 100.0, 69.5715808000383, 64.0719147336193, 7075.63080463549, 61.75997375866169]
_TB = [math.nan, 0.9957831631356981, math.nan, 4.143065389421082, math.nan, -1.5961450929589347, math.nan, 2.544946708512043, math.nan, 1.0999835991635507]
_TI = [math.nan, 0.0, math.nan, 0.0, math.nan, 0.0, math.nan, 0.0, math.nan, 0.0]
_CW = [0.0, math.nan, 9.094947017729282e-13, math.nan, 0.0, math.nan, 0.0, math.nan, 0.0, math.nan]
_CX = [0.0, math.nan, 0.0, math.nan, 0.0, math.nan, 0.0, math.nan, 0.0, math.nan]
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

def salomon_cd_10D():
    PB = SimpleNamespace()
    PB.name = _NAME
    PB.nvar = 10
    PB.var_types = ['C', 'D', 'C', 'D', 'C', 'D', 'C', 'D', 'C', 'D']
    PB.cardinality = [math.inf, 17.0, math.inf, 17.0, math.inf, 17.0, math.inf, 17.0, math.inf, 17.0]
    PB.f_star = 0.0
    PB.f_star_justification = 'E'
    PB.n_minimizers_documented = 1
    PB.description = 'DFO Benchmark Suite v2.1.0 core (Python): CD adaptation of salomon (F = f o T).'
    PB.x0_idx = [27.73538625240326, 16.0, 5629.057294336001, 2.0, 0.7855329457518486, 5.0, -34.18013340409669, 16.0, -3349.107061463994, 17.0]
    PB.x0_alt_idx = [[-52.68764269351959, 17.0, -1839.2800624639995, 2.0, 0.2531950789573303, 3.0, -81.70711625808664, 4.0, -8447.797512128967, 3.0], [-12.476128220558167, 16.0, 8930.848546816, 3.0, -0.12595203480279415, 17.0, 29.81044481761532, 10.0, 4237.879059178682, 6.0], [-92.89915716648102, 16.0, 1462.511190015999, 2.0, -0.6582899015973124, 14.0, -17.716538036374658, 15.0, -119.79609890409847, 9.0], [53.58564555644989, 17.0, -3615.9582126079995, 2.0, 0.6476612680208945, 12.0, -65.24352089036464, 4.0, -4477.471256986877, 12.0], [-26.83738338947296, 16.0, 4638.519918592002, 3.0, 0.10139696913234031, 8.0, 56.214324180943095, 9.0, -8835.14641506966, 15.0], [13.374131083488464, 16.0, -314.16696012800094, 1.0, -0.27775014462778425, 6.0, 8.687341326953117, 14.0, 3850.53015623799, 1.0], [-67.04889786243439, 16.0, 7940.311171072001, 2.0, 0.8750102719560288, 4.0, -48.779925522642614, 3.0, -507.1450018447904, 4.0], [33.001179814338684, 15.0, 2987.6242923519994, 2.0, 0.4958631581959042, 1.0, -96.30690837663259, 8.0, -4864.82015992757, 7.0], [-47.42184913158417, 16.0, -4480.713064448, 3.0, -0.036474708598614125, 15.0, 25.150936694675124, 13.0, -9963.510610592544, 10.0]]
    PB.suite_release = '2.1.0'
    PB.base_function = 'salomon'
    PB.m = 17
    PB.pattern = 'CD'
    PB.strategy = 'moderate'
    PB.scale_factors_C = [1.0, 100.0, 0.01, 1.0, 100.0]
    PB.lb_orig = [-100.0, -100.0, -63.779328, -70.41232819658475, -69.91407690731508, -62.55383214873429, -100.0, -100.0, -100.0, -100.0]
    PB.ub_orig = [63.399169921875, 61.59325306101712, 100.0, 100.0, 100.0, 100.0, 69.5715808000383, 64.0719147336193, 70.75630804635489, 61.75997375866169]
    PB.lb_C_work = [-100.0, -100.0, -6377.9328, -70.41232819658475, -0.6991407690731507, -62.55383214873429, -100.0, -100.0, -10000.0, -100.0]
    PB.ub_C_work = [63.399169921875, 61.59325306101712, 10000.0, 100.0, 0.9999999999999999, 100.0, 69.5715808000383, 64.0719147336193, 7075.63080463549, 61.75997375866169]
    PB.ProblemData = _problem_data()
    PB.func_idx = _objective
    PB.func_F = lambda xr: _objective(_cell2idx(xr))
    PB.idx2cell = _idx2cell
    PB.cell2idx = _cell2idx
    PB.decode = lambda xr: _decode(_cell2idx(xr))
    return PB
