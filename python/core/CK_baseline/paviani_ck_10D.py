"""paviani_ck_10D  DFO Benchmark Suite v2.1.0 core, Python version (self-contained; standard library only).
PAVIANI_CK_10D  DFO Benchmark Suite v2.1.0 core: paviani, dimension 10, pattern CK, strategy baseline.

 Copyright (C) 2026 J. F. A. Madeira

 Objective F(z) = f(T(z)), where f is the original function (embedded below, translated from the MATLAB original) and T acts on the finite
 coordinates only: uniform grid of 17 levels on the interval of the domain, piecewise-affine map with |T - Id| <=
 half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
 Global optimum of F on this domain: f* = -45.77845205382886 (A: analytic reduction, 1-D root to machine precision); the optimal solution is listed in
 reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
 Variables: x1:C x2:K x3:C x4:K x5:C x6:K x7:C x8:K x9:C x10:K


Python version of the MATLAB wrapper paviani_ck_10D.m: same instance, same values (level indices 1..m as in MATLAB).
Usage:  from paviani_ck_10D import paviani_ck_10D;  PB = paviani_ck_10D();  f = PB.func_idx(PB.x0_idx)
"""
# Copyright (C) 2026 J. F. A. Madeira. LGPL-3.0-or-later.
import math
from types import SimpleNamespace

_NAME = 'paviani_ck_10D'
_N = 10
_TYPES = ['C', 'K', 'C', 'K', 'C', 'K', 'C', 'K', 'C', 'K']
_CARD = [math.inf, 17.0, math.inf, 17.0, math.inf, 17.0, math.inf, 17.0, math.inf, 17.0]
_LB = [3.1781763354492187, 3.3917369212010366, 3.20312627328, 3.172785281132861, 3.4030945768731646, 3.197009175449039, 3.545224663737264, 3.2314999192762484, 3.5588647308292924, 3.395080028701464]
_UB = [9.999, 9.999, 9.999, 9.999, 9.999, 9.999, 9.999, 9.999, 9.999, 9.999]
_LBW = [3.1781763354492187, 3.3917369212010366, 3.20312627328, 3.172785281132861, 3.4030945768731646, 3.197009175449039, 3.545224663737264, 3.2314999192762484, 3.5588647308292924, 3.395080028701464]
_UBW = [9.999, 9.999, 9.999, 9.999, 9.999, 9.999, 9.999, 9.999, 9.999, 9.999]
_TB = [math.nan, 9.17309211515013, math.nan, 9.145723160141607, math.nan, 9.14875114693113, math.nan, 9.153062489909532, math.nan, 9.173510003587683]
_TI = [math.nan, 9.351, math.nan, 9.351, math.nan, 9.351, math.nan, 9.351, math.nan, 9.351]
_CW = [9.351, math.nan, 9.351, math.nan, 9.351, math.nan, 9.351, math.nan, 9.351, math.nan]
_CX = [9.351, math.nan, 9.351, math.nan, 9.351, math.nan, 9.351, math.nan, 9.351, math.nan]
_M = 17
_PERM = {1: [3, 9, 7, 15, 4, 13, 8, 17, 11, 2, 10, 12, 16, 6, 14, 1, 5], 3: [4, 6, 9, 14, 13, 15, 12, 10, 11, 3, 17, 16, 1, 7, 5, 8, 2], 5: [6, 3, 14, 7, 4, 2, 17, 9, 13, 11, 12, 5, 1, 15, 10, 16, 8], 7: [14, 16, 5, 17, 13, 8, 6, 10, 12, 7, 9, 1, 15, 11, 2, 4, 3], 9: [2, 12, 15, 9, 5, 3, 4, 17, 11, 8, 14, 16, 1, 13, 7, 10, 6]}

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

def paviani_ck_10D():
    PB = SimpleNamespace()
    PB.name = _NAME
    PB.nvar = 10
    PB.var_types = ['C', 'K', 'C', 'K', 'C', 'K', 'C', 'K', 'C', 'K']
    PB.cardinality = [math.inf, 17.0, math.inf, 17.0, math.inf, 17.0, math.inf, 17.0, math.inf, 17.0]
    PB.f_star = -45.778452053828865
    PB.f_star_justification = 'A'
    PB.n_minimizers_documented = 1
    PB.description = 'DFO Benchmark Suite v2.1.0 core (Python): CK adaptation of paviani (F = f o T).'
    PB.x0_idx = [4.460411252857447, 7.0, 5.195132780056166, 1.0, 4.969064609711692, 12.0, 9.824289615362725, 10.0, 4.107231636299898, 13.0]
    PB.x0_alt_idx = [[8.776713728080988, 8.0, 9.925060893853287, 2.0, 9.443972646588877, 10.0, 8.015446009005348, 14.0, 8.612203385778972, 15.0], [5.339658053365946, 7.0, 6.565180923362918, 3.0, 7.97215903977545, 8.0, 5.828282171252964, 4.0, 6.9686906044216075, 2.0], [7.098151654382944, 7.0, 3.4662625039785984, 1.0, 5.905673268592963, 5.0, 4.019438564895588, 9.0, 5.325177823064244, 5.0], [3.621130216008425, 8.0, 7.883036756448462, 2.0, 4.433859661779537, 3.0, 8.642038892253328, 14.0, 3.68166504170688, 8.0], [9.216337128335237, 7.0, 4.784118337064141, 3.0, 8.85470659280499, 17.0, 6.833195285895952, 3.0, 8.466113360769429, 11.0], [5.779281453620195, 7.0, 9.253084899755214, 3.0, 7.382892985991562, 14.0, 5.024351679538575, 8.0, 6.8226005794120645, 14.0], [7.537775054637194, 8.0, 6.154166480370893, 1.0, 5.316407214809075, 12.0, 9.268631775501309, 13.0, 5.1790877980547005, 16.0], [4.040770734432936, 6.0, 4.0990942654107645, 2.0, 3.8445936079956495, 9.0, 7.459788169143932, 2.0, 9.963536117117249, 3.0], [8.357073209656477, 7.0, 7.524214623677645, 3.0, 8.319501644872833, 7.0, 5.650944562786555, 7.0, 8.04054676617641, 6.0]]
    PB.suite_release = '2.1.0'
    PB.base_function = 'paviani'
    PB.m = 17
    PB.pattern = 'CK'
    PB.strategy = 'baseline'
    PB.scale_factors_C = [1.0, 1.0, 1.0, 1.0, 1.0]
    PB.lb_orig = [3.1781763354492187, 3.3917369212010366, 3.20312627328, 3.172785281132861, 3.4030945768731646, 3.197009175449039, 3.545224663737264, 3.2314999192762484, 3.5588647308292924, 3.395080028701464]
    PB.ub_orig = [9.999, 9.999, 9.999, 9.999, 9.999, 9.999, 9.999, 9.999, 9.999, 9.999]
    PB.lb_C_work = [3.1781763354492187, 3.3917369212010366, 3.20312627328, 3.172785281132861, 3.4030945768731646, 3.197009175449039, 3.545224663737264, 3.2314999192762484, 3.5588647308292924, 3.395080028701464]
    PB.ub_C_work = [9.999, 9.999, 9.999, 9.999, 9.999, 9.999, 9.999, 9.999, 9.999, 9.999]
    PB.ProblemData = _problem_data()
    PB.func_idx = _objective
    PB.func_F = lambda xr: _objective(_cell2idx(xr))
    PB.idx2cell = _idx2cell
    PB.cell2idx = _cell2idx
    PB.decode = lambda xr: _decode(_cell2idx(xr))
    return PB
