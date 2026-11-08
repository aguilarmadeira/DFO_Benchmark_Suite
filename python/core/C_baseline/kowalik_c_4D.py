"""kowalik_c_4D  DFO Benchmark Suite v2.1.0 core, Python version (self-contained; standard library only).
KOWALIK_C_4D  DFO Benchmark Suite v2.1.0 core: kowalik, dimension 4, pattern C, strategy baseline.

 Copyright (C) 2026 J. F. A. Madeira

 Objective F(z) = f(T(z)), where f is the original function (embedded below, translated from the MATLAB original) and T acts on the finite
 coordinates only: uniform grid of 0 levels on the interval of the domain, piecewise-affine map with |T - Id| <=
 half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
 Global optimum of F on this domain: f* = 0.0003074859878056062 (L: literature value at a refined point); the optimal solution is listed in
 reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
 Variables: x1:C x2:C x3:C x4:C


Python version of the MATLAB wrapper kowalik_c_4D.m: same instance, same values (level indices 1..m as in MATLAB).
Usage:  from kowalik_c_4D import kowalik_c_4D;  PB = kowalik_c_4D();  f = PB.func_idx(PB.x0_idx)
"""
# Copyright (C) 2026 J. F. A. Madeira. LGPL-3.0-or-later.
import math
from types import SimpleNamespace

_NAME = 'kowalik_c_4D'
_N = 4
_TYPES = ['C', 'C', 'C', 'C']
_CARD = [math.inf, math.inf, math.inf, math.inf]
_LB = [0.059979711914062504, 0.08174305746075294, 0.0, 0.05673411078717201]
_UB = [0.42, 0.42, 0.3474645888, 0.42]
_LBW = [0.059979711914062504, 0.08174305746075294, 0.0, 0.05673411078717201]
_UBW = [0.42, 0.42, 0.3474645888, 0.42]
_TB = [math.nan, math.nan, math.nan, math.nan]
_TI = [math.nan, math.nan, math.nan, math.nan]
_CW = [0.1928334529007549, 0.1908362411132952, 0.1231173000994744, 0.1357659904379078]
_CX = [0.1928334529007549, 0.1908362411132952, 0.1231173000994744, 0.1357659904379078]
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

_a = [0.1957, 0.1947, 0.1735, 0.16, 0.0844, 0.0627, 0.0456, 0.0342, 0.0323, 0.0235, 0.0246]
_b = [0.25, 0.5, 1, 2, 4, 6, 8, 10, 12, 14, 16]
def _orig(x):
    s = 0.0
    for a, b in zip(_a, _b): s += (a - x[0] * (1 + x[1] * b) / (1 + x[2] * b + x[3] * b ** 2)) ** 2
    return s

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

def kowalik_c_4D():
    PB = SimpleNamespace()
    PB.name = _NAME
    PB.nvar = 4
    PB.var_types = ['C', 'C', 'C', 'C']
    PB.cardinality = [math.inf, math.inf, math.inf, math.inf]
    PB.f_star = 0.00030748598780560617
    PB.f_star_justification = 'L'
    PB.n_minimizers_documented = 1
    PB.description = 'DFO Benchmark Suite v2.1.0 core (Python): C adaptation of kowalik (F = f o T).'
    PB.x0_idx = [0.10515803908109665, 0.316218377483386, 0.213037488685056, 0.41924351126777837]
    PB.x0_alt_idx = [[0.33298337763547897, 0.32596240463472237, 0.12130683724185601, 0.37703144000981126], [0.15156690434217454, 0.33570643178605863, 0.283086349787136, 0.39291770338646553], [0.2443846348643303, 0.3083303555037328, 0.19135569834393598, 0.4088039667631198], [0.06085866769552231, 0.3292104136851678, 0.016900677599232, 0.3738541873344804], [0.415253638780117, 0.3389544408365041, 0.258736031404032, 0.3897404507111347], [0.23805615323781967, 0.3152903748975444, 0.086949538701312, 0.39836442225846125], [0.32665489600896835, 0.32008505492439254, 0.328784892506112, 0.4142506856351155], [0.149457410466671, 0.3409651131058275, 0.15699839980339197, 0.3793009062064761], [0.3698995204567909, 0.3173010471668678, 0.065267748360192, 0.3951871695831304]]
    PB.suite_release = '2.1.0'
    PB.base_function = 'kowalik'
    PB.m = 0
    PB.pattern = 'C'
    PB.strategy = 'baseline'
    PB.scale_factors_C = [1.0, 1.0, 1.0, 1.0]
    PB.lb_orig = [0.059979711914062504, 0.08174305746075294, 0.0, 0.05673411078717201]
    PB.ub_orig = [0.42, 0.42, 0.3474645888, 0.42]
    PB.lb_C_work = [0.059979711914062504, 0.08174305746075294, 0.0, 0.05673411078717201]
    PB.ub_C_work = [0.42, 0.42, 0.3474645888, 0.42]
    PB.ProblemData = _problem_data()
    PB.func_idx = _objective
    PB.func_F = lambda xr: _objective(_cell2idx(xr))
    PB.idx2cell = _idx2cell
    PB.cell2idx = _cell2idx
    PB.decode = lambda xr: _decode(_cell2idx(xr))
    return PB
