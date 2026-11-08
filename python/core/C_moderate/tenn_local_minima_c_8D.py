"""tenn_local_minima_c_8D  DFO Benchmark Suite v2.1.0 core, Python version (self-contained; standard library only).
TENN_LOCAL_MINIMA_C_8D  DFO Benchmark Suite v2.1.0 core: tenn_local_minima, dimension 8, pattern C, strategy moderate.

 Copyright (C) 2026 J. F. A. Madeira

 Objective F(z) = f(T(z)), where f is the original function (embedded below, translated from the MATLAB original) and T acts on the finite
 coordinates only: uniform grid of 0 levels on the interval of the domain, piecewise-affine map with |T - Id| <=
 half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
 Global optimum of F on this domain: f* = 5.889542894199896e-32 (E: exact by structure); the optimal solution is listed in
 reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
 Variables: x1:C x2:C x3:C x4:C x5:C x6:C x7:C x8:C


Python version of the MATLAB wrapper tenn_local_minima_c_8D.m: same instance, same values (level indices 1..m as in MATLAB).
Usage:  from tenn_local_minima_c_8D import tenn_local_minima_c_8D;  PB = tenn_local_minima_c_8D();  f = PB.func_idx(PB.x0_idx)
"""
# Copyright (C) 2026 J. F. A. Madeira. LGPL-3.0-or-later.
import math
from types import SimpleNamespace

_NAME = 'tenn_local_minima_c_8D'
_N = 8
_TYPES = ['C', 'C', 'C', 'C', 'C', 'C', 'C', 'C']
_CARD = [math.inf, math.inf, math.inf, math.inf, math.inf, math.inf, math.inf, math.inf]
_LB = [-10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0]
_UB = [7.1711669921875, 6.722013920642178, 7.2627328, 6.580008329862557, 6.216201079161259, 7.20804593676692, 6.952313789346391, 6.55068177807107]
_LBW = [-10.0, -999.9999999999999, -0.1, -10.0, -1000.0, -0.1, -10.0, -1000.0]
_UBW = [7.1711669921875005, 672.2013920642177, 0.07262732800000002, 6.580008329862558, 621.6201079161258, 0.07208045936766921, 6.952313789346391, 655.0681778071071]
_TB = [math.nan, math.nan, math.nan, math.nan, math.nan, math.nan, math.nan, math.nan]
_TI = [math.nan, math.nan, math.nan, math.nan, math.nan, math.nan, math.nan, math.nan]
_CW = [1.0000000000000018, 100.00000000000011, 0.010000000000000009, 1.0, 100.0, 0.009999999999999995, 1.0, 100.0]
_CX = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
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

def tenn_local_minima_c_8D():
    PB = SimpleNamespace()
    PB.name = _NAME
    PB.nvar = 8
    PB.var_types = ['C', 'C', 'C', 'C', 'C', 'C', 'C', 'C']
    PB.cardinality = [math.inf, math.inf, math.inf, math.inf, math.inf, math.inf, math.inf, math.inf]
    PB.f_star = 5.889542894199896e-32
    PB.f_star_justification = 'E'
    PB.n_minimizers_documented = 1
    PB.description = 'DFO Benchmark Suite v2.1.0 core (Python): C adaptation of tenn_local_minima (F = f o T).'
    PB.x0_idx = [-7.576920282840728, -61.05930066078679, -0.08254392459264001, -1.5200956980128186, -293.35863066014053, -0.0038166572127911608, -3.3612351453892826, 385.0548681459097]
    PB.x0_alt_idx = [[3.2892088294029254, 42.1630074913254, 0.03760469569536004, -0.7950224474357395, -801.4094082717291, -0.02519943618747196, -8.112575688424084, -815.8890480147511], [-5.363449537754058, 90.33341796231105, -0.04774225526783999, -2.390183598705315, 444.95976558116104, -0.059740848377340935, 3.0359946053634577, -279.48191005510705], [-0.9365080475807179, -45.002497170458014, 0.07240636502016004, -1.6651103481282359, -63.09101203042769, -0.08112362735202172, -1.7153459376713442, 256.9252279045372], [-9.68977872133255, 3.167913300527516, -0.014266363822079992, -0.9400370975511549, -424.94012701997656, 0.06855582547074385, -6.466686480706146, -857.1511355501083], [6.810639560222626, 106.39022145263982, -0.09961331478528, -2.866660306227396, -946.2819648699323, 0.034014413280874845, 5.675628240382791, -320.74399759046423], [-1.8420188069343553, -10.595061119754064, 0.020535305502720022, -2.1415870556503167, 300.0872089829577, 0.012631634306194045, 0.9242876973479905, 215.66314036917993], [2.584922683238984, 13.107839270731006, -0.06481164546047999, -1.4165138050732367, -207.9635686286308, -0.00875114466848674, -4.820797272988209, -985.2807757914808], [-6.218654143810272, 61.27824974171665, 0.05533697482752001, -0.6914405544961593, -569.8126836181796, -0.043292556858355716, -9.57213781602301, -448.8736378318366], [4.6474749684333805, -55.70703283067712, -0.03000997613567999, -2.286601705765732, 530.3548276126708, -0.0646753358330365, 2.570176905065926, 87.53350012780743]]
    PB.suite_release = '2.1.0'
    PB.base_function = 'tenn_local_minima'
    PB.m = 0
    PB.pattern = 'C'
    PB.strategy = 'moderate'
    PB.scale_factors_C = [1.0, 100.0, 0.01, 1.0, 100.0, 0.01, 1.0, 100.0]
    PB.lb_orig = [-10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0]
    PB.ub_orig = [7.1711669921875, 6.722013920642178, 7.2627328, 6.580008329862557, 6.216201079161259, 7.20804593676692, 6.952313789346391, 6.55068177807107]
    PB.lb_C_work = [-10.0, -999.9999999999999, -0.1, -10.0, -1000.0, -0.1, -10.0, -1000.0]
    PB.ub_C_work = [7.1711669921875005, 672.2013920642177, 0.07262732800000002, 6.580008329862558, 621.6201079161258, 0.07208045936766921, 6.952313789346391, 655.0681778071071]
    PB.ProblemData = _problem_data()
    PB.func_idx = _objective
    PB.func_F = lambda xr: _objective(_cell2idx(xr))
    PB.idx2cell = _idx2cell
    PB.cell2idx = _cell2idx
    PB.decode = lambda xr: _decode(_cell2idx(xr))
    return PB
