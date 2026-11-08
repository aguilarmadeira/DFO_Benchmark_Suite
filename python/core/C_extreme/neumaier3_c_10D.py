"""neumaier3_c_10D  DFO Benchmark Suite v2.1.0 core, Python version (self-contained; standard library only).
NEUMAIER3_C_10D  DFO Benchmark Suite v2.1.0 core: neumaier3, dimension 10, pattern C, strategy extreme.

 Copyright (C) 2026 J. F. A. Madeira

 Objective F(z) = f(T(z)), where f is the original function (embedded below, translated from the MATLAB original) and T acts on the finite
 coordinates only: uniform grid of 0 levels on the interval of the domain, piecewise-affine map with |T - Id| <=
 half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
 Global optimum of F on this domain: f* = -210 (E: exact by structure); the optimal solution is listed in
 reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
 Variables: x1:C x2:C x3:C x4:C x5:C x6:C x7:C x8:C x9:C x10:C


Python version of the MATLAB wrapper neumaier3_c_10D.m: same instance, same values (level indices 1..m as in MATLAB).
Usage:  from neumaier3_c_10D import neumaier3_c_10D;  PB = neumaier3_c_10D();  f = PB.func_idx(PB.x0_idx)
"""
# Copyright (C) 2026 J. F. A. Madeira. LGPL-3.0-or-later.
import math
from types import SimpleNamespace

_NAME = 'neumaier3_c_10D'
_N = 10
_TYPES = ['C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C']
_CARD = [math.inf, math.inf, math.inf, math.inf, math.inf, math.inf, math.inf, math.inf, math.inf, math.inf]
_LB = [-100.0, -100.0, -100.0, -72.69804248229904, -66.21159756847209, -71.1692167641189, -100.0, -100.0, -100.0, -100.0]
_UB = [61.813232421875, 69.88954935731343, 72.739328, 100.0, 100.0, 100.0, 62.20825900073035, 69.96665157572455, 61.65233829210158, 65.62204272417893]
_LBW = [-100.0, -100.0, -100.0, -72.69804248229906, -66.21159756847209, -7116921676.411891, -10000000000.0, -10000000000.0, -10000000000.0, -10000000000.0]
_UBW = [61.813232421875, 69.88954935731343, 72.739328, 100.0, 100.0, 10000000000.0, 6220825900.073034, 6996665157.572454, 6165233829.210157, 6562204272.417892]
_TB = [math.nan, math.nan, math.nan, math.nan, math.nan, math.nan, math.nan, math.nan, math.nan, math.nan]
_TI = [math.nan, math.nan, math.nan, math.nan, math.nan, math.nan, math.nan, math.nan, math.nan, math.nan]
_CW = [10.0, 18.0, 24.0, 27.999999999999986, 30.000000000000014, 3000000000.000001, 2799999999.999998, 2400000000.0, 1800000000.0, 1000000000.0]
_CX = [10.0, 18.0, 24.0, 28.0, 30.0, 30.0, 28.0, 24.0, 18.0, 10.0]
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
    d = len(x)
    s1 = 0.0
    for v in x: s1 += (v - 1) ** 2
    s2 = 0.0
    for i in range(1, d): s2 += x[i] * x[i - 1]
    return s1 - s2

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

def neumaier3_c_10D():
    PB = SimpleNamespace()
    PB.name = _NAME
    PB.nvar = 10
    PB.var_types = ['C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C']
    PB.cardinality = [math.inf, math.inf, math.inf, math.inf, math.inf, math.inf, math.inf, math.inf, math.inf, math.inf]
    PB.f_star = -210.0
    PB.f_star_justification = 'E'
    PB.n_minimizers_documented = 1
    PB.description = 'DFO Benchmark Suite v2.1.0 core (Python): C adaptation of neumaier3 (F = f o T).'
    PB.x0_idx = [31.55226171016693, -89.20226458131387, -83.91451377664, 75.54463371762529, -41.86053706820758, 3564598404.771863, 4827545108.330551, -1190677265.6115952, -9326393231.576431, 1541049719.5351238]
    PB.x0_alt_idx = [[-48.09018862247467, -84.3083269455065, 36.31205851136002, 83.09702624600573, 70.9036046330173, 1437643640.2473068, 281223662.2893219, 3425853888.865368, 1981431632.433401, 4435989110.790092], [-8.26896345615387, -96.1936040610387, -49.090265251839995, 90.64941877438613, 33.81506633261442, -689311124.2772474, -5215962399.9173975, -8015114624.4036255, -2143911594.3026485, -9211582305.126198], [-87.91141378879547, -93.1640236198246, 71.13630703616002, 98.20181130276657, -18.25873996795123, -4125161128.50922, -9762283845.958628, -2506500155.0518885, -6269254821.038702, -6316642913.871227], [41.98163020610809, -82.67701440023737, -15.592654766080003, 81.58654774032965, -55.347278268354074, -6252115893.033775, 1856093182.813406, 3002114314.2998486, 5740077683.459583, -3421703522.6162567], [-37.66082012653351, -94.56229151576957, -94.36178833408, 89.13894026871006, 56.054566341946895, 8636567458.638111, -2690228263.227825, -8438854198.969143, 1614734456.7235298, -526764131.36128426], [2.1604050397872925, -89.6683538799622, 19.23159375872001, 93.23881049840227, 18.96602804154405, 5200717454.406135, -7236549709.269056, -2930239729.617407, -2510608770.012521, 2368175259.8936863], [-77.48204529285431, -87.26022583694586, -59.53753980928, 76.62354693596532, -33.1077782590216, 3073762689.8815813, 3430962703.337494, 1686291424.8595543, -6635951996.748572, 5263114651.148655], [21.59695541858673, -99.14550295247808, 54.05584228352001, 84.17593946434576, 94.64163144236605, -362087314.3503933, -1115358742.7037373, -9754677088.409437, 5373380507.74971, -8384456764.767634], [-58.04549491405487, -88.65849373289083, -24.713291284479993, 91.72833199272617, 42.5678251418004, -2489042078.8749475, -5661680188.744969, -4246062619.0577, 546529640.5252094, -5489517373.512664]]
    PB.suite_release = '2.1.0'
    PB.base_function = 'neumaier3'
    PB.m = 0
    PB.pattern = 'C'
    PB.strategy = 'extreme'
    PB.scale_factors_C = [1.0, 1.0, 1.0, 1.0, 1.0, 100000000.0, 100000000.0, 100000000.0, 100000000.0, 100000000.0]
    PB.lb_orig = [-100.0, -100.0, -100.0, -72.69804248229904, -66.21159756847209, -71.1692167641189, -100.0, -100.0, -100.0, -100.0]
    PB.ub_orig = [61.813232421875, 69.88954935731343, 72.739328, 100.0, 100.0, 100.0, 62.20825900073035, 69.96665157572455, 61.65233829210158, 65.62204272417893]
    PB.lb_C_work = [-100.0, -100.0, -100.0, -72.69804248229906, -66.21159756847209, -7116921676.411891, -10000000000.0, -10000000000.0, -10000000000.0, -10000000000.0]
    PB.ub_C_work = [61.813232421875, 69.88954935731343, 72.739328, 100.0, 100.0, 10000000000.0, 6220825900.073034, 6996665157.572454, 6165233829.210157, 6562204272.417892]
    PB.ProblemData = _problem_data()
    PB.func_idx = _objective
    PB.func_F = lambda xr: _objective(_cell2idx(xr))
    PB.idx2cell = _idx2cell
    PB.cell2idx = _cell2idx
    PB.decode = lambda xr: _decode(_cell2idx(xr))
    return PB
