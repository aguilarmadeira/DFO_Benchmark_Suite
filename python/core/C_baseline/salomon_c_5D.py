"""salomon_c_5D  DFO Benchmark Suite v2.1.0 core, Python version (self-contained; standard library only).
SALOMON_C_5D  DFO Benchmark Suite v2.1.0 core: salomon, dimension 5, pattern C, strategy baseline.

 Copyright (C) 2026 J. F. A. Madeira

 Objective F(z) = f(T(z)), where f is the original function (embedded below, translated from the MATLAB original) and T acts on the finite
 coordinates only: uniform grid of 0 levels on the interval of the domain, piecewise-affine map with |T - Id| <=
 half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
 Global optimum of F on this domain: f* = 0 (E: exact by structure); the optimal solution is listed in
 reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
 Variables: x1:C x2:C x3:C x4:C x5:C


Python version of the MATLAB wrapper salomon_c_5D.m: same instance, same values (level indices 1..m as in MATLAB).
Usage:  from salomon_c_5D import salomon_c_5D;  PB = salomon_c_5D();  f = PB.func_idx(PB.x0_idx)
"""
# Copyright (C) 2026 J. F. A. Madeira. LGPL-3.0-or-later.
import math
from types import SimpleNamespace

_NAME = 'salomon_c_5D'
_N = 5
_TYPES = ['C', 'C', 'C', 'C', 'C']
_CARD = [math.inf, math.inf, math.inf, math.inf, math.inf]
_LB = [-70.399169921875, -100.0, -66.579328, -100.0, -71.18680418004234]
_UB = [100.0, 66.2599197276838, 100.0, 72.41232819658475, 100.0]
_LBW = [-70.399169921875, -100.0, -66.579328, -100.0, -71.18680418004234]
_UBW = [100.0, 66.2599197276838, 100.0, 72.41232819658475, 100.0]
_TB = [math.nan, math.nan, math.nan, math.nan, math.nan]
_TI = [math.nan, math.nan, math.nan, math.nan, math.nan]
_CW = [0.0, 0.0, 0.0, 1.4210854715202004e-14, 0.0]
_CX = [0.0, 0.0, 0.0, 0.0, 0.0]
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

def salomon_c_5D():
    PB = SimpleNamespace()
    PB.name = _NAME
    PB.nvar = 5
    PB.var_types = ['C', 'C', 'C', 'C', 'C']
    PB.cardinality = [math.inf, math.inf, math.inf, math.inf, math.inf]
    PB.f_star = 0.0
    PB.f_star_justification = 'E'
    PB.n_minimizers_documented = 1
    PB.description = 'DFO Benchmark Suite v2.1.0 core (Python): C adaptation of salomon (F = f o T).'
    PB.x0_idx = [-22.39120066165924, 2.9336677234951196, 22.227443343360008, 51.516040189418504, 62.830213066842816]
    PB.x0_alt_idx = [[21.539835333824158, 7.723048127173243, -53.73273022464, 59.05593792396303, 9.197683132148839], [-64.32537138462067, -3.9083042817593423, 55.80983586816001, 66.59583565850755, -29.001025238604413], [91.4301198720932, -0.9434497461490707, -20.15033769984, 50.00806064250958, 87.13848203004943], [4.56648051738739, 9.31950826173265, 88.11289915392001, 57.54795837705413, 48.939773659296165], [48.49751651287079, -2.311844147199963, 12.152725585919995, 65.08785611159865, -6.095830963028774], [-37.367690205574036, 2.4775362564781602, -38.220863201280004, 69.18094345320853, -44.29453933378202], [70.46303451061249, 4.834215502732491, 45.73511811072002, 52.59316843721055, 71.84496793487182], [-16.899821162223816, -6.797136906200123, -4.638470676480011, 60.133066171755075, 33.64625956411855], [27.031214833259583, 3.465821101681584, 79.31751063552002, 67.6729639062996, -19.98627037057541]]
    PB.suite_release = '2.1.0'
    PB.base_function = 'salomon'
    PB.m = 0
    PB.pattern = 'C'
    PB.strategy = 'baseline'
    PB.scale_factors_C = [1.0, 1.0, 1.0, 1.0, 1.0]
    PB.lb_orig = [-70.399169921875, -100.0, -66.579328, -100.0, -71.18680418004234]
    PB.ub_orig = [100.0, 66.2599197276838, 100.0, 72.41232819658475, 100.0]
    PB.lb_C_work = [-70.399169921875, -100.0, -66.579328, -100.0, -71.18680418004234]
    PB.ub_C_work = [100.0, 66.2599197276838, 100.0, 72.41232819658475, 100.0]
    PB.ProblemData = _problem_data()
    PB.func_idx = _objective
    PB.func_F = lambda xr: _objective(_cell2idx(xr))
    PB.idx2cell = _idx2cell
    PB.cell2idx = _cell2idx
    PB.decode = lambda xr: _decode(_cell2idx(xr))
    return PB
