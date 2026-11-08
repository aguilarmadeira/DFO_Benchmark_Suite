"""fifteenn_local_minima_c_10D  DFO Benchmark Suite v2.1.0 core, Python version (self-contained; standard library only).
FIFTEENN_LOCAL_MINIMA_C_10D  DFO Benchmark Suite v2.1.0 core: fifteenn_local_minima, dimension 10, pattern C, strategy moderate.

 Copyright (C) 2026 J. F. A. Madeira

 Objective F(z) = f(T(z)), where f is the original function (embedded below, translated from the MATLAB original) and T acts on the finite
 coordinates only: uniform grid of 0 levels on the interval of the domain, piecewise-affine map with |T - Id| <=
 half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
 Global optimum of F on this domain: f* = 1.349783804395672e-32 (E: exact by structure); the optimal solution is listed in
 reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
 Variables: x1:C x2:C x3:C x4:C x5:C x6:C x7:C x8:C x9:C x10:C


Python version of the MATLAB wrapper fifteenn_local_minima_c_10D.m: same instance, same values (level indices 1..m as in MATLAB).
Usage:  from fifteenn_local_minima_c_10D import fifteenn_local_minima_c_10D;  PB = fifteenn_local_minima_c_10D();  f = PB.func_idx(PB.x0_idx)
"""
# Copyright (C) 2026 J. F. A. Madeira. LGPL-3.0-or-later.
import math
from types import SimpleNamespace

_NAME = 'fifteenn_local_minima_c_10D'
_N = 10
_TYPES = ['C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C']
_CARD = [math.inf, math.inf, math.inf, math.inf, math.inf, math.inf, math.inf, math.inf, math.inf, math.inf]
_LB = [-10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0]
_UB = [6.2250732421875, 6.02105370116344, 6.1091328, 6.35551853394419, 7.26909364114473, 6.0565666468260915, 6.143317249554005, 6.779490642336999, 7.1391468726884195, 7.334618065521342]
_LBW = [-10.0, -1000.0, -0.1, -10.0, -999.9999999999999, -0.1, -10.0, -1000.0, -0.1, -10.0]
_UBW = [6.2250732421874995, 602.1053701163439, 0.061091328, 6.355518533944191, 726.9093641144729, 0.06056566646826092, 6.143317249554006, 677.9490642336998, 0.0713914687268842, 7.334618065521341]
_TB = [math.nan, math.nan, math.nan, math.nan, math.nan, math.nan, math.nan, math.nan, math.nan, math.nan]
_TI = [math.nan, math.nan, math.nan, math.nan, math.nan, math.nan, math.nan, math.nan, math.nan, math.nan]
_CW = [1.0, 100.0, 0.009999999999999981, 1.0, 100.00000000000011, 0.010000000000000009, 1.0, 100.0, 0.009999999999999995, 1.0]
_CX = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
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

def fifteenn_local_minima_c_10D():
    PB = SimpleNamespace()
    PB.name = _NAME
    PB.nvar = 10
    PB.var_types = ['C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C']
    PB.cardinality = [math.inf, math.inf, math.inf, math.inf, math.inf, math.inf, math.inf, math.inf, math.inf, math.inf]
    PB.f_star = 1.3497838043956716e-32
    PB.f_star_justification = 'E'
    PB.n_minimizers_documented = 1
    PB.description = 'DFO Benchmark Suite v2.1.0 core (Python): C adaptation of fifteenn_local_minima (F = f o T).'
    PB.x0_idx = [4.181094288825987, 209.4540311212088, 0.04902880935936002, 1.3146673406419431, 166.4098560021872, -0.03780547011174782, 5.650441697460966, 129.96744827168095, -0.04101947237943091, -7.484635969745376]
    PB.x0_alt_idx = [[-4.161191904544831, 97.37258547520923, -0.024428836208640003, -0.5858697984544463, -218.9335558250092, -0.057757416832655985, 1.1258441292468682, 673.7902751839879, -0.09219603240941121, -4.4546848928944405], [0.02183479070663452, 143.52376897650333, -0.07314285379584001, 0.12938611410871026, -759.9712754611738, -0.07770936355356414, -3.3987534389672307, -543.7563777816227, 0.03513260125725323, -1.4247338160435046], [-8.15408102273941, 242.41916219356153, 0.008047175516160007, 0.8446420266718668, 567.3226986102806, 0.04967614243377261, -8.869672067199964, 0.06644913068453207, -0.00860624236870021, 1.6052172608074322], [5.155549371242522, 112.75631330897409, -0.04190402347008001, 1.5598979392350234, 26.28497897411569, 0.029724195712864415, 2.693188384902731, 543.8892760429915, -0.05234508599465364, 4.635168337658367], [-3.115435230731965, 158.90749681026819, 0.03928600584192002, -0.013665068403920344, -528.9067686446824, 0.009772248991956264, -1.8314091833113668, -585.5889029287232, -0.09608392962060708, -9.648886738924617], [1.0675914645195004, 205.05868031156183, -0.00942801174527999, 0.37461671270179053, -914.2501804718787, -0.022457818787972317, -6.356006751525465, -41.76607601641604, 0.031244704046057364, -6.61893566207368], [-7.108324348926544, 69.53536368077835, -0.08288565731328, 1.0898726252649489, 257.34948579060745, -0.042409765508880475, 4.260532640558596, 413.98827690199505, -0.012494139579896063, -3.588984585222744], [3.159104812145232, 168.43075689783655, 0.023047999979520006, -0.4836903823739984, -127.99392603658896, -0.07463983328880905, -0.26406492765550027, -715.4899020697196, -0.0562329832058495, -0.5590335083718063], [-5.1594141840934755, 214.58194039913042, -0.050409645588480004, 0.23156553018915993, -669.0316456727537, -0.0945917800097172, -4.788662495869601, -171.66707515741257, 0.07109565046081495, 2.4709175684791287]]
    PB.suite_release = '2.1.0'
    PB.base_function = 'fifteenn_local_minima'
    PB.m = 0
    PB.pattern = 'C'
    PB.strategy = 'moderate'
    PB.scale_factors_C = [1.0, 100.0, 0.01, 1.0, 100.0, 0.01, 1.0, 100.0, 0.01, 1.0]
    PB.lb_orig = [-10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0]
    PB.ub_orig = [6.2250732421875, 6.02105370116344, 6.1091328, 6.35551853394419, 7.26909364114473, 6.0565666468260915, 6.143317249554005, 6.779490642336999, 7.1391468726884195, 7.334618065521342]
    PB.lb_C_work = [-10.0, -1000.0, -0.1, -10.0, -999.9999999999999, -0.1, -10.0, -1000.0, -0.1, -10.0]
    PB.ub_C_work = [6.2250732421874995, 602.1053701163439, 0.061091328, 6.355518533944191, 726.9093641144729, 0.06056566646826092, 6.143317249554006, 677.9490642336998, 0.0713914687268842, 7.334618065521341]
    PB.ProblemData = _problem_data()
    PB.func_idx = _objective
    PB.func_F = lambda xr: _objective(_cell2idx(xr))
    PB.idx2cell = _idx2cell
    PB.cell2idx = _cell2idx
    PB.decode = lambda xr: _decode(_cell2idx(xr))
    return PB
