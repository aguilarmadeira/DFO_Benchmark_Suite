"""paviani_c_10D  DFO Benchmark Suite v2.1.0 core, Python version (self-contained; standard library only).
PAVIANI_C_10D  DFO Benchmark Suite v2.1.0 core: paviani, dimension 10, pattern C, strategy extreme.

 Copyright (C) 2026 J. F. A. Madeira

 Objective F(z) = f(T(z)), where f is the original function (embedded below, translated from the MATLAB original) and T acts on the finite
 coordinates only: uniform grid of 0 levels on the interval of the domain, piecewise-affine map with |T - Id| <=
 half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
 Global optimum of F on this domain: f* = -45.77845205382886 (A: analytic reduction, 1-D root to machine precision); the optimal solution is listed in
 reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
 Variables: x1:C x2:C x3:C x4:C x5:C x6:C x7:C x8:C x9:C x10:C


Python version of the MATLAB wrapper paviani_c_10D.m: same instance, same values (level indices 1..m as in MATLAB).
Usage:  from paviani_c_10D import paviani_c_10D;  PB = paviani_c_10D();  f = PB.func_idx(PB.x0_idx)
"""
# Copyright (C) 2026 J. F. A. Madeira. LGPL-3.0-or-later.
import math
from types import SimpleNamespace

_NAME = 'paviani_c_10D'
_N = 10
_TYPES = ['C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C']
_CARD = [math.inf, math.inf, math.inf, math.inf, math.inf, math.inf, math.inf, math.inf, math.inf, math.inf]
_LB = [3.1781763354492187, 3.3917369212010366, 3.20312627328, 3.172785281132861, 3.4030945768731646, 3.197009175449039, 3.545224663737264, 3.2314999192762484, 3.5588647308292924, 3.395080028701464]
_UB = [9.999, 9.999, 9.999, 9.999, 9.999, 9.999, 9.999, 9.999, 9.999, 9.999]
_LBW = [3.178176335449219, 3.391736921201037, 3.2031262732800005, 3.1727852811328616, 3.403094576873164, 319700917.5449038, 354522466.3737263, 323149991.92762476, 355886473.08292925, 339508002.87014645]
_UBW = [9.999, 9.999, 9.999, 9.999, 9.999, 999900000.0, 999900000.0, 999900000.0, 999900000.0000001, 999900000.0]
_TB = [math.nan, math.nan, math.nan, math.nan, math.nan, math.nan, math.nan, math.nan, math.nan, math.nan]
_TI = [math.nan, math.nan, math.nan, math.nan, math.nan, math.nan, math.nan, math.nan, math.nan, math.nan]
_CW = [9.351, 9.351, 9.351, 9.351, 9.351, 935100000.0, 935100000.0, 935100000.0, 935100000.0000001, 935100000.0]
_CX = [9.351, 9.351, 9.351, 9.351, 9.351, 9.351, 9.351, 9.351, 9.351, 9.351]
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

def paviani_c_10D():
    PB = SimpleNamespace()
    PB.name = _NAME
    PB.nvar = 10
    PB.var_types = ['C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C']
    PB.cardinality = [math.inf, math.inf, math.inf, math.inf, math.inf, math.inf, math.inf, math.inf, math.inf, math.inf]
    PB.f_star = -45.778452053828865
    PB.f_star_justification = 'A'
    PB.n_minimizers_documented = 1
    PB.description = 'DFO Benchmark Suite v2.1.0 core (Python): C adaptation of paviani (F = f o T).'
    PB.x0_idx = [4.4604112528574476, 6.0140983077568215, 5.195132780056167, 3.320625000158722, 4.969064609711692, 796490442.7979071, 982428961.5362723, 709525998.8690771, 410723163.6299898, 822461451.5891619]
    PB.x0_alt_idx = [[8.776713728080988, 6.204430989039096, 9.925060893853288, 3.6191475097301713, 9.443972646588877, 711968663.3212383, 801544600.9005346, 893340733.5299366, 861220338.5778972, 937892633.0137618], [5.339658053365946, 5.742194477353572, 6.565180923362919, 3.9176700193016205, 7.97215903977545, 627446883.8445696, 582828217.1252964, 437799869.3704153, 696869060.4421608, 393717063.4406482], [7.098151654382945, 5.8600194705283135, 3.466262503978599, 3.2609204982444324, 5.905673268592963, 490911701.61302775, 401943856.48955876, 657134359.5212958, 532517782.3064244, 509148244.86524796], [3.6211302160084253, 6.267875216133188, 7.883036756448462, 3.5594430078158816, 4.433859661779537, 406389922.136359, 864203889.2253326, 876468849.6721766, 368166504.170688, 624579426.2898479], [9.216337128335237, 5.805638704447663, 4.784118337064141, 3.8579655173873304, 8.85470659280499, 998042378.4730403, 683319528.5895951, 420927985.5126552, 846611336.0769429, 740010607.7144477], [5.779281453620196, 5.995971385729939, 9.253084899755214, 4.020020594011831, 7.382892985991562, 861507196.2414985, 502435167.9538574, 640262475.6635358, 682260057.9412065, 855441789.1390475], [7.537775054637194, 6.089627149535501, 6.154166480370893, 3.3632710729546433, 5.316407214809075, 776985416.7648298, 926863177.5501306, 824077210.3243953, 517908779.8054701, 970872970.5636473], [4.040770734432936, 5.627390637849977, 4.099094265410765, 3.6617935825260925, 3.844593607995649, 640450234.5332878, 745978816.9143931, 368536346.164874, 996353611.711725, 426697400.9905338], [8.357073209656477, 6.035246383454852, 7.524214623677646, 3.9603160920975418, 8.319501644872833, 555928455.056619, 565094456.2786554, 587870836.3157547, 804054676.6176411, 542128582.4151337]]
    PB.suite_release = '2.1.0'
    PB.base_function = 'paviani'
    PB.m = 0
    PB.pattern = 'C'
    PB.strategy = 'extreme'
    PB.scale_factors_C = [1.0, 1.0, 1.0, 1.0, 1.0, 100000000.0, 100000000.0, 100000000.0, 100000000.0, 100000000.0]
    PB.lb_orig = [3.1781763354492187, 3.3917369212010366, 3.20312627328, 3.172785281132861, 3.4030945768731646, 3.197009175449039, 3.545224663737264, 3.2314999192762484, 3.5588647308292924, 3.395080028701464]
    PB.ub_orig = [9.999, 9.999, 9.999, 9.999, 9.999, 9.999, 9.999, 9.999, 9.999, 9.999]
    PB.lb_C_work = [3.178176335449219, 3.391736921201037, 3.2031262732800005, 3.1727852811328616, 3.403094576873164, 319700917.5449038, 354522466.3737263, 323149991.92762476, 355886473.08292925, 339508002.87014645]
    PB.ub_C_work = [9.999, 9.999, 9.999, 9.999, 9.999, 999900000.0, 999900000.0, 999900000.0, 999900000.0000001, 999900000.0]
    PB.ProblemData = _problem_data()
    PB.func_idx = _objective
    PB.func_F = lambda xr: _objective(_cell2idx(xr))
    PB.idx2cell = _idx2cell
    PB.cell2idx = _cell2idx
    PB.decode = lambda xr: _decode(_cell2idx(xr))
    return PB
