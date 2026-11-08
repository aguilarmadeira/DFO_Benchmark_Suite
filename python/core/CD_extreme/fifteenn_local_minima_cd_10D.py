"""fifteenn_local_minima_cd_10D  DFO Benchmark Suite v2.1.0 core, Python version (self-contained; standard library only).
FIFTEENN_LOCAL_MINIMA_CD_10D  DFO Benchmark Suite v2.1.0 core: fifteenn_local_minima, dimension 10, pattern CD, strategy extreme.

 Copyright (C) 2026 J. F. A. Madeira

 Objective F(z) = f(T(z)), where f is the original function (embedded below, translated from the MATLAB original) and T acts on the finite
 coordinates only: uniform grid of 17 levels on the interval of the domain, piecewise-affine map with |T - Id| <=
 half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
 Global optimum of F on this domain: f* = 1.349783804395672e-32 (E: exact by structure); the optimal solution is listed in
 reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
 Variables: x1:C x2:D x3:C x4:D x5:C x6:D x7:C x8:D x9:C x10:D


Python version of the MATLAB wrapper fifteenn_local_minima_cd_10D.m: same instance, same values (level indices 1..m as in MATLAB).
Usage:  from fifteenn_local_minima_cd_10D import fifteenn_local_minima_cd_10D;  PB = fifteenn_local_minima_cd_10D();  f = PB.func_idx(PB.x0_idx)
"""
# Copyright (C) 2026 J. F. A. Madeira. LGPL-3.0-or-later.
import math
from types import SimpleNamespace

_NAME = 'fifteenn_local_minima_cd_10D'
_N = 10
_TYPES = ['C', 'D', 'C', 'D', 'C', 'D', 'C', 'D', 'C', 'D']
_CARD = [math.inf, 17.0, math.inf, 17.0, math.inf, 17.0, math.inf, 17.0, math.inf, 17.0]
_LB = [-10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0]
_UB = [6.2250732421875, 6.02105370116344, 6.1091328, 6.35551853394419, 7.26909364114473, 6.0565666468260915, 6.143317249554005, 6.779490642336999, 7.1391468726884195, 7.334618065521342]
_LBW = [-10.0, -10.0, -10.0, -10.0, -1000000000.0, -10.0, -1000000000.0, -10.0, -1000000000.0000001, -10.0]
_UBW = [6.2250732421874995, 6.02105370116344, 6.109132800000001, 6.35551853394419, 726909364.1144729, 6.0565666468260915, 614331724.9554005, 6.779490642336999, 713914687.2688421, 7.334618065521342]
_TB = [math.nan, 1.0144744195498632, math.nan, 1.2444189920866329, math.nan, 1.0388895696929374, math.nan, 0.48718165146062375, math.nan, 0.8341362909508376]
_TI = [math.nan, 1.0, math.nan, 1.0, math.nan, 1.0, math.nan, 1.0, math.nan, 1.0]
_CW = [1.0, math.nan, 1.0, math.nan, 100000000.0, math.nan, 100000000.0, math.nan, 100000000.00000012, math.nan]
_CX = [1.0, math.nan, 1.0, math.nan, 1.0, math.nan, 1.0, math.nan, 1.0, math.nan]
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
    d = len(x)
    s = 0.0
    for i in range(d - 1): s += ((x[i] - 1) ** 2) * (1 + 10 * math.sin(3 * math.pi * x[i + 1]) ** 2)
    f = 0.1 * s
    return f + 0.1 * (math.sin(3 * math.pi * x[0]) ** 2 + (x[d - 1] - 1) ** 2 * (1 + math.sin(2 * math.pi * x[d - 1]) ** 2))

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

def fifteenn_local_minima_cd_10D():
    PB = SimpleNamespace()
    PB.name = _NAME
    PB.nvar = 10
    PB.var_types = ['C', 'D', 'C', 'D', 'C', 'D', 'C', 'D', 'C', 'D']
    PB.cardinality = [math.inf, 17.0, math.inf, 17.0, math.inf, 17.0, math.inf, 17.0, math.inf, 17.0]
    PB.f_star = 1.3497838043956716e-32
    PB.f_star_justification = 'E'
    PB.n_minimizers_documented = 1
    PB.description = 'DFO Benchmark Suite v2.1.0 core (Python): CD adaptation of fifteenn_local_minima (F = f o T).'
    PB.x0_idx = [4.181094288825987, 13.0, 4.902880935936002, 12.0, 166409856.00218725, 7.0, 565044169.7460966, 12.0, -410194723.79430914, 3.0]
    PB.x0_alt_idx = [[-4.161191904544831, 12.0, -2.4428836208639986, 10.0, -218933555.82500923, 5.0, 112584412.92468691, 17.0, -921960324.0941122, 6.0], [0.02183479070663452, 12.0, -7.314285379584, 11.0, -759971275.461174, 3.0, -339875343.89672315, 5.0, 351326012.5725323, 9.0], [-8.15408102273941, 13.0, 0.8047175516160028, 12.0, 567322698.6102805, 16.0, -886967206.7199963, 11.0, -86062423.68700206, 12.0], [5.155549371242522, 12.0, -4.190402347008, 12.0, 26284978.97411573, 14.0, 269318838.490273, 16.0, -523450859.94653636, 15.0], [-3.115435230731965, 13.0, 3.9286005841920026, 11.0, -528906768.6446825, 12.0, -183140918.3311367, 5.0, -960839296.2060708, 1.0], [1.0675914645195004, 13.0, -0.9428011745279985, 11.0, -914250180.4718789, 9.0, -635600675.1525465, 10.0, 312447040.4605738, 4.0], [-7.108324348926544, 12.0, -8.288565731328, 12.0, 257349485.7906072, 7.0, 426053264.05585957, 14.0, -124941395.79896069, 7.0], [3.159104812145232, 13.0, 2.304799997952003, 10.0, -127993926.03658903, 4.0, -26406492.765550137, 4.0, -562329832.058495, 10.0], [-5.1594141840934755, 13.0, -5.040964558848, 11.0, -669031645.6727538, 2.0, -478866249.5869601, 9.0, 710956504.6081494, 13.0]]
    PB.suite_release = '2.1.0'
    PB.base_function = 'fifteenn_local_minima'
    PB.m = 17
    PB.pattern = 'CD'
    PB.strategy = 'extreme'
    PB.scale_factors_C = [1.0, 1.0, 100000000.0, 100000000.0, 100000000.0]
    PB.lb_orig = [-10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0]
    PB.ub_orig = [6.2250732421875, 6.02105370116344, 6.1091328, 6.35551853394419, 7.26909364114473, 6.0565666468260915, 6.143317249554005, 6.779490642336999, 7.1391468726884195, 7.334618065521342]
    PB.lb_C_work = [-10.0, -10.0, -10.0, -10.0, -1000000000.0, -10.0, -1000000000.0, -10.0, -1000000000.0000001, -10.0]
    PB.ub_C_work = [6.2250732421874995, 6.02105370116344, 6.109132800000001, 6.35551853394419, 726909364.1144729, 6.0565666468260915, 614331724.9554005, 6.779490642336999, 713914687.2688421, 7.334618065521342]
    PB.ProblemData = _problem_data()
    PB.func_idx = _objective
    PB.func_F = lambda xr: _objective(_cell2idx(xr))
    PB.idx2cell = _idx2cell
    PB.cell2idx = _cell2idx
    PB.decode = lambda xr: _decode(_cell2idx(xr))
    return PB
