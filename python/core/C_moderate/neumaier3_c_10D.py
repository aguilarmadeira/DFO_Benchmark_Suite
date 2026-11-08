"""neumaier3_c_10D  DFO Benchmark Suite v2.1.0 core, Python version (self-contained; standard library only).
NEUMAIER3_C_10D  DFO Benchmark Suite v2.1.0 core: neumaier3, dimension 10, pattern C, strategy moderate.

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
_LBW = [-100.0, -10000.0, -1.0, -72.69804248229906, -6621.15975684721, -0.711692167641189, -100.0, -10000.0, -1.0, -100.0]
_UBW = [61.813232421875, 6988.9549357313435, 0.72739328, 100.0, 10000.0, 1.0, 62.20825900073035, 6996.6651575724545, 0.6165233829210157, 65.62204272417893]
_TB = [math.nan, math.nan, math.nan, math.nan, math.nan, math.nan, math.nan, math.nan, math.nan, math.nan]
_TI = [math.nan, math.nan, math.nan, math.nan, math.nan, math.nan, math.nan, math.nan, math.nan, math.nan]
_CW = [10.0, 1800.0, 0.24, 27.999999999999986, 3000.0, 0.29999999999999993, 27.999999999999986, 2399.999999999998, 0.17999999999999994, 10.0]
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
    PB.x0_idx = [31.55226171016693, -8920.226458131387, -0.8391451377664, 75.54463371762529, -4186.0537068207595, 0.35645984047718626, 48.27545108330551, -1190.6772656115954, -0.9326393231576432, 15.41049719535124]
    PB.x0_alt_idx = [[-48.09018862247467, -8430.832694550649, 0.36312058511360035, 83.09702624600573, 7090.360463301729, 0.1437643640247308, 2.8122366228932094, 3425.8538888653675, 0.19814316324334014, 44.35989110790092], [-8.26896345615387, -9619.36040610387, -0.49090265251839993, 90.64941877438613, 3381.5066332614424, -0.06893111242772476, -52.159623999173974, -8015.114624403625, -0.21439115943026488, -92.11582305126197], [-87.91141378879547, -9316.40236198246, 0.7113630703616001, 98.20181130276657, -1825.8739967951233, -0.412516112850922, -97.62283845958629, -2506.500155051889, -0.6269254821038701, -63.16642913871227], [41.98163020610809, -8267.701440023737, -0.15592654766079994, 81.58654774032965, -5534.727826835408, -0.6252115893033775, 18.560931828134073, 3002.114314299848, 0.5740077683459583, -34.217035226162565], [-37.66082012653351, -9456.229151576958, -0.9436178833408, 89.13894026871006, 5605.456634194688, 0.863656745863811, -26.902282632278244, -8438.854198969144, 0.161473445672353, -5.267641313612842], [2.1604050397872925, -8966.83538799622, 0.19231593758720011, 93.23881049840227, 1896.6028041544032, 0.5200717454406133, -72.36549709269056, -2930.2397296174067, -0.2510608770012521, 23.681752598936853], [-77.48204529285431, -8726.022583694587, -0.5953753980928, 76.62354693596532, -3310.7778259021607, 0.3073762689881582, 34.30962703337494, 1686.2914248595534, -0.6635951996748573, 52.63114651148655], [21.59695541858673, -9914.550295247807, 0.5405584228352001, 84.17593946434576, 9464.163144236605, -0.03620873143503933, -11.153587427037365, -9754.677088409437, 0.537338050774971, -83.84456764767634], [-58.04549491405487, -8865.849373289082, -0.2471329128448, 91.72833199272617, 4256.78251418004, -0.24890420788749473, -56.6168018874497, -4246.0626190577, 0.05465296405252107, -54.89517373512664]]
    PB.suite_release = '2.1.0'
    PB.base_function = 'neumaier3'
    PB.m = 0
    PB.pattern = 'C'
    PB.strategy = 'moderate'
    PB.scale_factors_C = [1.0, 100.0, 0.01, 1.0, 100.0, 0.01, 1.0, 100.0, 0.01, 1.0]
    PB.lb_orig = [-100.0, -100.0, -100.0, -72.69804248229904, -66.21159756847209, -71.1692167641189, -100.0, -100.0, -100.0, -100.0]
    PB.ub_orig = [61.813232421875, 69.88954935731343, 72.739328, 100.0, 100.0, 100.0, 62.20825900073035, 69.96665157572455, 61.65233829210158, 65.62204272417893]
    PB.lb_C_work = [-100.0, -10000.0, -1.0, -72.69804248229906, -6621.15975684721, -0.711692167641189, -100.0, -10000.0, -1.0, -100.0]
    PB.ub_C_work = [61.813232421875, 6988.9549357313435, 0.72739328, 100.0, 10000.0, 1.0, 62.20825900073035, 6996.6651575724545, 0.6165233829210157, 65.62204272417893]
    PB.ProblemData = _problem_data()
    PB.func_idx = _objective
    PB.func_F = lambda xr: _objective(_cell2idx(xr))
    PB.idx2cell = _idx2cell
    PB.cell2idx = _cell2idx
    PB.decode = lambda xr: _decode(_cell2idx(xr))
    return PB
