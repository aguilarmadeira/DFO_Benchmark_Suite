"""paviani_c_10D  DFO Benchmark Suite v2.1.0 core, Python version (self-contained; standard library only).
PAVIANI_C_10D  DFO Benchmark Suite v2.1.0 core: paviani, dimension 10, pattern C, strategy moderate.

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
_LBW = [3.178176335449219, 339.1736921201037, 0.032031262732800006, 3.1727852811328616, 340.3094576873164, 0.031970091754490386, 3.5452246637372635, 323.14999192762485, 0.03558864730829291, 3.3950800287014644]
_UBW = [9.999, 999.9000000000001, 0.09999000000000001, 9.999, 999.9000000000001, 0.09999, 9.999, 999.9000000000001, 0.09999, 9.999]
_TB = [math.nan, math.nan, math.nan, math.nan, math.nan, math.nan, math.nan, math.nan, math.nan, math.nan]
_TI = [math.nan, math.nan, math.nan, math.nan, math.nan, math.nan, math.nan, math.nan, math.nan, math.nan]
_CW = [9.351, 935.1000000000001, 0.09351000000000001, 9.351, 935.1000000000001, 0.09351, 9.351, 935.1000000000001, 0.09350999999999998, 9.351]
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
    PB.x0_idx = [4.4604112528574476, 601.4098307756822, 0.05195132780056168, 3.320625000158722, 496.90646097116917, 0.07964904427979069, 9.824289615362723, 709.5259988690772, 0.04107231636299896, 8.22461451589162]
    PB.x0_alt_idx = [[8.776713728080988, 620.4430989039097, 0.09925060893853288, 3.6191475097301713, 944.3972646588877, 0.07119686633212383, 8.015446009005348, 893.3407335299366, 0.08612203385778969, 9.378926330137618], [5.339658053365946, 574.2194477353572, 0.0656518092336292, 3.9176700193016205, 797.215903977545, 0.06274468838445696, 5.8282821712529636, 437.79986937041537, 0.06968690604421607, 3.937170634406482], [7.098151654382945, 586.0019470528314, 0.03466262503978599, 3.2609204982444324, 590.5673268592964, 0.04909117016130277, 4.019438564895588, 657.1343595212959, 0.05325177823064242, 5.09148244865248], [3.6211302160084253, 626.7875216133189, 0.07883036756448461, 3.5594430078158816, 443.3859661779537, 0.0406389922136359, 8.642038892253328, 876.4688496721767, 0.03681665041706879, 6.245794262898479], [9.216337128335237, 580.5638704447664, 0.047841183370641416, 3.8579655173873304, 885.4706592804989, 0.09980423784730402, 6.833195285895951, 420.9279855126553, 0.08466113360769426, 7.4001060771444775], [5.779281453620196, 599.5971385729938, 0.09253084899755214, 4.020020594011831, 738.2892985991564, 0.08615071962414983, 5.024351679538575, 640.262475663536, 0.06822600579412064, 8.554417891390475], [7.537775054637194, 608.9627149535502, 0.061541664803708936, 3.3632710729546433, 531.6407214809075, 0.07769854167648296, 9.268631775501309, 824.0772103243953, 0.05179087798054699, 9.708729705636475], [4.040770734432936, 562.7390637849977, 0.04099094265410766, 3.6617935825260925, 384.4593607995649, 0.06404502345332877, 7.459788169143931, 368.5363461648741, 0.09963536117117247, 4.266974009905338], [8.357073209656477, 603.5246383454852, 0.07524214623677647, 3.9603160920975418, 831.9501644872835, 0.0555928455056619, 5.650944562786554, 587.8708363157548, 0.0804054676617641, 5.421285824151337]]
    PB.suite_release = '2.1.0'
    PB.base_function = 'paviani'
    PB.m = 0
    PB.pattern = 'C'
    PB.strategy = 'moderate'
    PB.scale_factors_C = [1.0, 100.0, 0.01, 1.0, 100.0, 0.01, 1.0, 100.0, 0.01, 1.0]
    PB.lb_orig = [3.1781763354492187, 3.3917369212010366, 3.20312627328, 3.172785281132861, 3.4030945768731646, 3.197009175449039, 3.545224663737264, 3.2314999192762484, 3.5588647308292924, 3.395080028701464]
    PB.ub_orig = [9.999, 9.999, 9.999, 9.999, 9.999, 9.999, 9.999, 9.999, 9.999, 9.999]
    PB.lb_C_work = [3.178176335449219, 339.1736921201037, 0.032031262732800006, 3.1727852811328616, 340.3094576873164, 0.031970091754490386, 3.5452246637372635, 323.14999192762485, 0.03558864730829291, 3.3950800287014644]
    PB.ub_C_work = [9.999, 999.9000000000001, 0.09999000000000001, 9.999, 999.9000000000001, 0.09999, 9.999, 999.9000000000001, 0.09999, 9.999]
    PB.ProblemData = _problem_data()
    PB.func_idx = _objective
    PB.func_F = lambda xr: _objective(_cell2idx(xr))
    PB.idx2cell = _idx2cell
    PB.cell2idx = _cell2idx
    PB.decode = lambda xr: _decode(_cell2idx(xr))
    return PB
