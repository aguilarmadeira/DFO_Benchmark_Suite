"""rastrigin_c_10D  DFO Benchmark Suite v2.1.0 core, Python version (self-contained; standard library only).
RASTRIGIN_C_10D  DFO Benchmark Suite v2.1.0 core: rastrigin, dimension 10, pattern C, strategy extreme.

 Copyright (C) 2026 J. F. A. Madeira

 Objective F(z) = f(T(z)), where f is the original function (embedded below, translated from the MATLAB original) and T acts on the finite
 coordinates only: uniform grid of 0 levels on the interval of the domain, piecewise-affine map with |T - Id| <=
 half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
 Global optimum of F on this domain: f* = 0 (E: exact by structure); the optimal solution is listed in
 reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
 Variables: x1:C x2:C x3:C x4:C x5:C x6:C x7:C x8:C x9:C x10:C


Python version of the MATLAB wrapper rastrigin_c_10D.m: same instance, same values (level indices 1..m as in MATLAB).
Usage:  from rastrigin_c_10D import rastrigin_c_10D;  PB = rastrigin_c_10D();  f = PB.func_idx(PB.x0_idx)
"""
# Copyright (C) 2026 J. F. A. Madeira. LGPL-3.0-or-later.
import math
from types import SimpleNamespace

_NAME = 'rastrigin_c_10D'
_N = 10
_TYPES = ['C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C']
_CARD = [math.inf, math.inf, math.inf, math.inf, math.inf, math.inf, math.inf, math.inf, math.inf, math.inf]
_LB = [-5.12, -5.12, -3.6955815936, -3.210139775093711, -3.775091646745441, -3.3681715906305802, -5.12, -5.12, -5.12, -5.12]
_UB = [3.7836375, 3.2332190011685213, 5.12, 5.12, 5.12, 5.12, 3.6885590546090206, 3.393660981729729, 3.7162186241472837, 3.23626238058141]
_LBW = [-5.12, -5.120000000000001, -3.6955815936000005, -3.210139775093711, -3.7750916467454414, -336817159.063058, -512000000.0, -511999999.99999994, -512000000.0, -511999999.99999994]
_UBW = [3.7836375, 3.233219001168522, 5.120000000000001, 5.12, 5.120000000000001, 512000000.0, 368855905.460902, 339366098.17297286, 371621862.4147284, 323626238.05814093]
_TB = [math.nan, math.nan, math.nan, math.nan, math.nan, math.nan, math.nan, math.nan, math.nan, math.nan]
_TI = [math.nan, math.nan, math.nan, math.nan, math.nan, math.nan, math.nan, math.nan, math.nan, math.nan]
_CW = [0.0, -8.881784197001252e-16, -4.440892098500626e-16, 0.0, -4.440892098500626e-16, 0.0, 0.0, 0.0, 5.960464477539063e-08, 0.0]
_CX = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
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
    for v in x: s += v ** 2 - 10 * math.cos(2 * math.pi * v)
    return 10 * len(x) + s

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

def rastrigin_c_10D():
    PB = SimpleNamespace()
    PB.name = _NAME
    PB.nvar = 10
    PB.var_types = ['C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C']
    PB.cardinality = [math.inf, math.inf, math.inf, math.inf, math.inf, math.inf, math.inf, math.inf, math.inf, math.inf]
    PB.f_star = 0.0
    PB.f_star_justification = 'E'
    PB.n_minimizers_documented = 1
    PB.description = 'DFO Benchmark Suite v2.1.0 core (Python): C adaptation of rastrigin (F = f o T).'
    PB.x0_idx = [-4.837413848876953, 1.9078568642661526, -2.5220513718599684, 1.7303512868527466, 1.5713195609152302, -324453868.6260904, -325537524.59203374, 153304313.48697108, -283087687.1594293, 214363370.65203935]
    PB.x0_alt_idx = [[3.301067303466797, 2.1484845721187433, 1.9210017513144337, 2.0946431137518884, -1.2154972810778957, 413866382.1565696, 305386947.4854981, -464460531.33106154, 335019132.18895185, -474208547.2722075], [-1.1855312805175782, 1.5641029959053085, -0.7448301225902081, 2.45893494065103, -3.200352441921919, 308392060.6161895, 6867695.991013765, -188533014.7486853, 109520169.1908642, -328147837.4094885], [1.109937762451172, 1.7130630055283413, 3.6982230005841923, 1.657492921472918, 2.8344092087250634, 138010464.28172964, -240015793.08280313, 87394501.83369094, -115978793.80722332, -182087127.5467695], [-3.4288305725097654, 2.2286938080696084, 0.9646874600407043, 2.0217847483720606, 0.8495540478810404, 32536142.741349578, -486899282.1566199, -485685724.91432124, -341477756.80531085, -36026417.68405044], [2.2055025329589846, 1.6443122318561727, -3.0552177466408965, 2.3860765752712023, -2.0101684014892767, -137845453.59311038, 144025189.9209119, -209758208.331945, 314974779.9224551, 110034292.17866868], [-2.281096051025391, 1.8849399397087634, 2.7419087093104655, 2.5838349955878783, 4.826554930306807, -243319775.13349038, -102858299.15290493, 66169308.25043124, 89475816.92436767, 256095002.04138774], [0.014372991943359281, 2.003344049921944, -1.2779964973711362, 1.7823929764097661, 2.039738088313683, 495000475.64916956, -401377550.6473893, 297412206.7627869, -136023146.07372004, -432476915.88285923], [-4.550480218505859, 1.4189624737085076, 4.5191299585802245, 2.146684803308907, 0.054882927469657616, 324618879.31470966, 229546921.43014264, -275668019.9852252, -361522109.07180756, -286416206.02014023], [2.7532849182128905, 1.9345932762497746, 0.4992247518986237, 2.5109766302080496, -2.7319339145234673, 219144557.77432954, -17336567.643674135, 259496.59715110064, 256584710.27657354, -140355496.15742117]]
    PB.suite_release = '2.1.0'
    PB.base_function = 'rastrigin'
    PB.m = 0
    PB.pattern = 'C'
    PB.strategy = 'extreme'
    PB.scale_factors_C = [1.0, 1.0, 1.0, 1.0, 1.0, 100000000.0, 100000000.0, 100000000.0, 100000000.0, 100000000.0]
    PB.lb_orig = [-5.12, -5.12, -3.6955815936, -3.210139775093711, -3.775091646745441, -3.3681715906305802, -5.12, -5.12, -5.12, -5.12]
    PB.ub_orig = [3.7836375, 3.2332190011685213, 5.12, 5.12, 5.12, 5.12, 3.6885590546090206, 3.393660981729729, 3.7162186241472837, 3.23626238058141]
    PB.lb_C_work = [-5.12, -5.120000000000001, -3.6955815936000005, -3.210139775093711, -3.7750916467454414, -336817159.063058, -512000000.0, -511999999.99999994, -512000000.0, -511999999.99999994]
    PB.ub_C_work = [3.7836375, 3.233219001168522, 5.120000000000001, 5.12, 5.120000000000001, 512000000.0, 368855905.460902, 339366098.17297286, 371621862.4147284, 323626238.05814093]
    PB.ProblemData = _problem_data()
    PB.func_idx = _objective
    PB.func_F = lambda xr: _objective(_cell2idx(xr))
    PB.idx2cell = _idx2cell
    PB.cell2idx = _cell2idx
    PB.decode = lambda xr: _decode(_cell2idx(xr))
    return PB
