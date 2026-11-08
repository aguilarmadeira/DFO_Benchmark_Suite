"""fifteenn_local_minima_ck_10D  DFO Benchmark Suite v2.1.0 core, Python version (self-contained; standard library only).
FIFTEENN_LOCAL_MINIMA_CK_10D  DFO Benchmark Suite v2.1.0 core: fifteenn_local_minima, dimension 10, pattern CK, strategy moderate.

 Copyright (C) 2026 J. F. A. Madeira

 Objective F(z) = f(T(z)), where f is the original function (embedded below, translated from the MATLAB original) and T acts on the finite
 coordinates only: uniform grid of 17 levels on the interval of the domain, piecewise-affine map with |T - Id| <=
 half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
 Global optimum of F on this domain: f* = 1.349783804395672e-32 (E: exact by structure); the optimal solution is listed in
 reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
 Variables: x1:C x2:K x3:C x4:K x5:C x6:K x7:C x8:K x9:C x10:K


Python version of the MATLAB wrapper fifteenn_local_minima_ck_10D.m: same instance, same values (level indices 1..m as in MATLAB).
Usage:  from fifteenn_local_minima_ck_10D import fifteenn_local_minima_ck_10D;  PB = fifteenn_local_minima_ck_10D();  f = PB.func_idx(PB.x0_idx)
"""
# Copyright (C) 2026 J. F. A. Madeira. LGPL-3.0-or-later.
import math
from types import SimpleNamespace

_NAME = 'fifteenn_local_minima_ck_10D'
_N = 10
_TYPES = ['C', 'K', 'C', 'K', 'C', 'K', 'C', 'K', 'C', 'K']
_CARD = [math.inf, 17.0, math.inf, 17.0, math.inf, 17.0, math.inf, 17.0, math.inf, 17.0]
_LB = [-10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0]
_UB = [6.2250732421875, 6.02105370116344, 6.1091328, 6.35551853394419, 7.26909364114473, 6.0565666468260915, 6.143317249554005, 6.779490642336999, 7.1391468726884195, 7.334618065521342]
_LBW = [-10.0, -10.0, -1000.0000000000001, -10.0, -0.09999999999999999, -10.0, -10.0, -10.0, -1000.0, -10.0]
_UBW = [6.2250732421874995, 6.02105370116344, 610.9132800000001, 6.35551853394419, 0.07269093641144729, 6.0565666468260915, 6.143317249554006, 6.779490642336999, 713.914687268842, 7.334618065521342]
_TB = [math.nan, 1.0144744195498632, math.nan, 1.2444189920866329, math.nan, 1.0388895696929374, math.nan, 0.48718165146062375, math.nan, 0.8341362909508376]
_TI = [math.nan, 1.0, math.nan, 1.0, math.nan, 1.0, math.nan, 1.0, math.nan, 1.0]
_CW = [1.0, math.nan, 99.99999999999989, math.nan, 0.010000000000000009, math.nan, 1.0, math.nan, 100.0, math.nan]
_CX = [1.0, math.nan, 1.0, math.nan, 1.0, math.nan, 1.0, math.nan, 1.0, math.nan]
_M = 17
_PERM = {1: [8, 12, 1, 9, 13, 7, 14, 6, 5, 4, 16, 2, 15, 11, 3, 17, 10], 3: [15, 14, 7, 6, 12, 1, 17, 9, 4, 5, 2, 10, 3, 13, 8, 11, 16], 5: [8, 9, 14, 11, 3, 16, 10, 15, 4, 12, 6, 7, 2, 1, 17, 5, 13], 7: [11, 17, 6, 16, 10, 12, 5, 1, 15, 8, 9, 13, 2, 3, 4, 7, 14], 9: [13, 12, 4, 5, 8, 16, 9, 7, 3, 15, 1, 11, 2, 6, 17, 14, 10]}

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

def fifteenn_local_minima_ck_10D():
    PB = SimpleNamespace()
    PB.name = _NAME
    PB.nvar = 10
    PB.var_types = ['C', 'K', 'C', 'K', 'C', 'K', 'C', 'K', 'C', 'K']
    PB.cardinality = [math.inf, 17.0, math.inf, 17.0, math.inf, 17.0, math.inf, 17.0, math.inf, 17.0]
    PB.f_star = 1.3497838043956716e-32
    PB.f_star_justification = 'E'
    PB.n_minimizers_documented = 1
    PB.description = 'DFO Benchmark Suite v2.1.0 core (Python): CK adaptation of fifteenn_local_minima (F = f o T).'
    PB.x0_idx = [4.181094288825987, 13.0, 490.28809359360014, 12.0, 0.016640985600218733, 7.0, 5.650441697460966, 12.0, -410.1947237943091, 3.0]
    PB.x0_alt_idx = [[-4.161191904544831, 12.0, -244.28836208639996, 10.0, -0.021893355582500926, 5.0, 1.1258441292468682, 17.0, -921.960324094112, 6.0], [0.02183479070663452, 12.0, -731.4285379584, 11.0, -0.07599712754611739, 3.0, -3.3987534389672307, 5.0, 351.32601257253236, 9.0], [-8.15408102273941, 13.0, 80.47175516160007, 12.0, 0.05673226986102807, 16.0, -8.869672067199964, 11.0, -86.06242368700202, 12.0], [5.155549371242522, 12.0, -419.04023470080006, 12.0, 0.002628497897411569, 14.0, 2.693188384902731, 16.0, -523.4508599465363, 15.0], [-3.115435230731965, 13.0, 392.8600584192003, 11.0, -0.052890676864468245, 12.0, -1.8314091833113668, 5.0, -960.8392962060707, 1.0], [1.0675914645195004, 13.0, -94.28011745279991, 11.0, -0.09142501804718788, 9.0, -6.356006751525465, 10.0, 312.4470404605736, 4.0], [-7.108324348926544, 12.0, -828.8565731328001, 12.0, 0.025734948579060737, 7.0, 4.260532640558596, 14.0, -124.94139579896068, 7.0], [3.159104812145232, 13.0, 230.47999979520011, 10.0, -0.012799392603658893, 4.0, -0.26406492765550027, 4.0, -562.329832058495, 10.0], [-5.1594141840934755, 13.0, -504.09645588480004, 11.0, -0.06690316456727538, 2.0, -4.788662495869601, 9.0, 710.9565046081495, 13.0]]
    PB.suite_release = '2.1.0'
    PB.base_function = 'fifteenn_local_minima'
    PB.m = 17
    PB.pattern = 'CK'
    PB.strategy = 'moderate'
    PB.scale_factors_C = [1.0, 100.0, 0.01, 1.0, 100.0]
    PB.lb_orig = [-10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0]
    PB.ub_orig = [6.2250732421875, 6.02105370116344, 6.1091328, 6.35551853394419, 7.26909364114473, 6.0565666468260915, 6.143317249554005, 6.779490642336999, 7.1391468726884195, 7.334618065521342]
    PB.lb_C_work = [-10.0, -10.0, -1000.0000000000001, -10.0, -0.09999999999999999, -10.0, -10.0, -10.0, -1000.0, -10.0]
    PB.ub_C_work = [6.2250732421874995, 6.02105370116344, 610.9132800000001, 6.35551853394419, 0.07269093641144729, 6.0565666468260915, 6.143317249554006, 6.779490642336999, 713.914687268842, 7.334618065521342]
    PB.ProblemData = _problem_data()
    PB.func_idx = _objective
    PB.func_F = lambda xr: _objective(_cell2idx(xr))
    PB.idx2cell = _idx2cell
    PB.cell2idx = _cell2idx
    PB.decode = lambda xr: _decode(_cell2idx(xr))
    return PB
