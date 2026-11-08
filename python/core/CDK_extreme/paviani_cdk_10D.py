"""paviani_cdk_10D  DFO Benchmark Suite v2.1.0 core, Python version (self-contained; standard library only).
PAVIANI_CDK_10D  DFO Benchmark Suite v2.1.0 core: paviani, dimension 10, pattern CDK, strategy extreme.

 Copyright (C) 2026 J. F. A. Madeira

 Objective F(z) = f(T(z)), where f is the original function (embedded below, translated from the MATLAB original) and T acts on the finite
 coordinates only: uniform grid of 17 levels on the interval of the domain, piecewise-affine map with |T - Id| <=
 half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
 Global optimum of F on this domain: f* = -45.77845205382886 (A: analytic reduction, 1-D root to machine precision); the optimal solution is listed in
 reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
 Variables: x1:C x2:K x3:D x4:C x5:D x6:K x7:C x8:K x9:D x10:C


Python version of the MATLAB wrapper paviani_cdk_10D.m: same instance, same values (level indices 1..m as in MATLAB).
Usage:  from paviani_cdk_10D import paviani_cdk_10D;  PB = paviani_cdk_10D();  f = PB.func_idx(PB.x0_idx)
"""
# Copyright (C) 2026 J. F. A. Madeira. LGPL-3.0-or-later.
import math
from types import SimpleNamespace

_NAME = 'paviani_cdk_10D'
_N = 10
_TYPES = ['C', 'K', 'D', 'C', 'D', 'K', 'C', 'K', 'D', 'C']
_CARD = [math.inf, 17.0, 17.0, math.inf, 17.0, 17.0, math.inf, 17.0, 17.0, math.inf]
_LB = [3.1781763354492187, 3.3917369212010366, 3.20312627328, 3.172785281132861, 3.4030945768731646, 3.197009175449039, 3.545224663737264, 3.2314999192762484, 3.5588647308292924, 3.395080028701464]
_UB = [9.999, 9.999, 9.999, 9.999, 9.999, 9.999, 9.999, 9.999, 9.999, 9.999]
_LBW = [3.178176335449219, 3.3917369212010366, 3.20312627328, 3.1727852811328616, 3.4030945768731646, 3.197009175449039, 354522466.3737263, 3.2314999192762484, 3.5588647308292924, 339508002.87014645]
_UBW = [9.999, 9.999, 9.999, 9.999, 9.999, 9.999, 999900000.0, 9.999, 9.999, 999900000.0]
_TB = [math.nan, 9.17309211515013, 9.14951578416, math.nan, 9.174511822109146, 9.14875114693113, math.nan, 9.153062489909532, 9.193983091353662, math.nan]
_TI = [math.nan, 9.351, 9.351, math.nan, 9.351, 9.351, math.nan, 9.351, 9.351, math.nan]
_CW = [9.351, math.nan, math.nan, 9.351, math.nan, math.nan, 935100000.0, math.nan, math.nan, 935100000.0]
_CX = [9.351, math.nan, math.nan, 9.351, math.nan, math.nan, 9.351, math.nan, math.nan, 9.351]
_M = 17
_PERM = {1: [3, 9, 7, 15, 4, 13, 8, 17, 11, 2, 10, 12, 16, 6, 14, 1, 5], 5: [6, 3, 14, 7, 4, 2, 17, 9, 13, 11, 12, 5, 1, 15, 10, 16, 8], 7: [14, 16, 5, 17, 13, 8, 6, 10, 12, 7, 9, 1, 15, 11, 2, 4, 3]}

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
    for v in x: s += (math.log(v - 2)) ** 2 + (math.log(10 - v)) ** 2
    p = 1.0
    for v in x: p *= v
    return s - p ** 0.2

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

def paviani_cdk_10D():
    PB = SimpleNamespace()
    PB.name = _NAME
    PB.nvar = 10
    PB.var_types = ['C', 'K', 'D', 'C', 'D', 'K', 'C', 'K', 'D', 'C']
    PB.cardinality = [math.inf, 17.0, 17.0, math.inf, 17.0, 17.0, math.inf, 17.0, 17.0, math.inf]
    PB.f_star = -45.778452053828865
    PB.f_star_justification = 'A'
    PB.n_minimizers_documented = 1
    PB.description = 'DFO Benchmark Suite v2.1.0 core (Python): CDK adaptation of paviani (F = f o T).'
    PB.x0_idx = [4.4604112528574476, 7.0, 6.0, 3.320625000158722, 5.0, 12.0, 982428961.5362723, 10.0, 2.0, 822461451.5891619]
    PB.x0_alt_idx = [[8.776713728080988, 8.0, 17.0, 3.6191475097301713, 16.0, 10.0, 801544600.9005346, 14.0, 14.0, 937892633.0137618], [5.339658053365946, 7.0, 9.0, 3.9176700193016205, 12.0, 8.0, 582828217.1252964, 4.0, 9.0, 393717063.4406482], [7.098151654382945, 7.0, 2.0, 3.2609204982444324, 7.0, 5.0, 401943856.48955876, 9.0, 5.0, 509148244.86524796], [3.6211302160084253, 8.0, 12.0, 3.5594430078158816, 4.0, 3.0, 864203889.2253326, 14.0, 1.0, 624579426.2898479], [9.216337128335237, 7.0, 5.0, 3.8579655173873304, 14.0, 17.0, 683319528.5895951, 3.0, 13.0, 740010607.7144477], [5.779281453620196, 7.0, 15.0, 4.020020594011831, 11.0, 14.0, 502435167.9538574, 8.0, 9.0, 855441789.1390475], [7.537775054637194, 8.0, 8.0, 3.3632710729546433, 6.0, 12.0, 926863177.5501306, 13.0, 5.0, 970872970.5636473], [4.040770734432936, 6.0, 3.0, 3.6617935825260925, 2.0, 9.0, 745978816.9143931, 2.0, 17.0, 426697400.9905338], [8.357073209656477, 7.0, 11.0, 3.9603160920975418, 13.0, 7.0, 565094456.2786554, 7.0, 12.0, 542128582.4151337]]
    PB.suite_release = '2.1.0'
    PB.base_function = 'paviani'
    PB.m = 17
    PB.pattern = 'CDK'
    PB.strategy = 'extreme'
    PB.scale_factors_C = [1.0, 1.0, 100000000.0, 100000000.0]
    PB.lb_orig = [3.1781763354492187, 3.3917369212010366, 3.20312627328, 3.172785281132861, 3.4030945768731646, 3.197009175449039, 3.545224663737264, 3.2314999192762484, 3.5588647308292924, 3.395080028701464]
    PB.ub_orig = [9.999, 9.999, 9.999, 9.999, 9.999, 9.999, 9.999, 9.999, 9.999, 9.999]
    PB.lb_C_work = [3.178176335449219, 3.3917369212010366, 3.20312627328, 3.1727852811328616, 3.4030945768731646, 3.197009175449039, 354522466.3737263, 3.2314999192762484, 3.5588647308292924, 339508002.87014645]
    PB.ub_C_work = [9.999, 9.999, 9.999, 9.999, 9.999, 9.999, 999900000.0, 9.999, 9.999, 999900000.0]
    PB.ProblemData = _problem_data()
    PB.func_idx = _objective
    PB.func_F = lambda xr: _objective(_cell2idx(xr))
    PB.idx2cell = _idx2cell
    PB.cell2idx = _cell2idx
    PB.decode = lambda xr: _decode(_cell2idx(xr))
    return PB
