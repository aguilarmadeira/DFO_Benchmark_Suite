"""rastrigin_c_10D  DFO Benchmark Suite v2.1.0 core, Python version (self-contained; standard library only).
RASTRIGIN_C_10D  DFO Benchmark Suite v2.1.0 core: rastrigin, dimension 10, pattern C, strategy moderate.

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
_LBW = [-5.12, -512.0, -0.036955815936, -3.210139775093711, -377.50916467454414, -0.0336817159063058, -5.12, -512.0, -0.0512, -5.119999999999999]
_UBW = [3.7836375, 323.3219001168522, 0.05120000000000001, 5.12, 512.0, 0.0512, 3.68855905460902, 339.3660981729729, 0.037162186241472835, 3.2362623805814095]
_TB = [math.nan, math.nan, math.nan, math.nan, math.nan, math.nan, math.nan, math.nan, math.nan, math.nan]
_TI = [math.nan, math.nan, math.nan, math.nan, math.nan, math.nan, math.nan, math.nan, math.nan, math.nan]
_CW = [0.0, 0.0, 0.0, 0.0, -1.1368683772161603e-13, 0.0, 0.0, 0.0, 0.0, 8.881784197001252e-16]
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
    PB.x0_idx = [-4.837413848876953, 190.78568642661537, -0.02522051371859968, 1.7303512868527466, 157.13195609152302, -0.03244538686260904, -3.255375245920338, 153.30431348697107, -0.028308768715942935, 2.1436337065203945]
    PB.x0_alt_idx = [[3.301067303466797, 214.84845721187446, 0.019210017513144335, 2.0946431137518884, -121.54972810778958, 0.04138663821565696, 3.0538694748549817, -464.4605313310616, 0.03350191321889518, -4.742085472722075], [-1.1855312805175782, 156.41029959053094, -0.007448301225902076, 2.45893494065103, -320.03524419219195, 0.030839206061618947, 0.068676959910138, -188.53301474868533, 0.010952016919086413, -3.2814783740948847], [1.109937762451172, 171.30630055283427, 0.03698223000584193, 1.657492921472918, 283.4409208725063, 0.013801046428172967, -2.400157930828031, 87.39450183369092, -0.011597879380722334, -1.8208712754676948], [-3.4288305725097654, 222.86938080696098, 0.009646874600407043, 2.0217847483720606, 84.95540478810398, 0.0032536142741349597, -4.8689928215661995, -485.68572491432127, -0.034147775680531095, -0.36026417684050394], [2.2055025329589846, 164.43122318561734, -0.03055217746640896, 2.3860765752712023, -201.01684014892768, -0.013784545359311038, 1.4402518992091196, -209.75820833194507, 0.03149747799224551, 1.1003429217866874], [-2.281096051025391, 188.49399397087643, 0.027419087093104658, 2.5838349955878783, 482.6554930306807, -0.02433197751334904, -1.0285829915290492, 66.16930825043119, 0.008947581692436757, 2.5609500204138778], [0.014372991943359281, 200.33440499219444, -0.01277996497371136, 1.7823929764097661, 203.97380883136822, 0.04950004756491696, -4.013775506473893, 297.4122067627869, -0.013602314607372004, -4.324769158828592], [-4.550480218505859, 141.89624737085092, 0.045191299585802244, 2.146684803308907, 5.488292746965726, 0.03246188793147097, 2.295469214301427, -275.6680199852252, -0.03615221090718076, -2.864162060201402], [2.7532849182128905, 193.4593276249775, 0.004992247518986234, 2.5109766302080496, -273.19339145234676, 0.021914455777432958, -0.17336567643674172, 0.2594965971510419, 0.025658471027657344, -1.4035549615742111]]
    PB.suite_release = '2.1.0'
    PB.base_function = 'rastrigin'
    PB.m = 0
    PB.pattern = 'C'
    PB.strategy = 'moderate'
    PB.scale_factors_C = [1.0, 100.0, 0.01, 1.0, 100.0, 0.01, 1.0, 100.0, 0.01, 1.0]
    PB.lb_orig = [-5.12, -5.12, -3.6955815936, -3.210139775093711, -3.775091646745441, -3.3681715906305802, -5.12, -5.12, -5.12, -5.12]
    PB.ub_orig = [3.7836375, 3.2332190011685213, 5.12, 5.12, 5.12, 5.12, 3.6885590546090206, 3.393660981729729, 3.7162186241472837, 3.23626238058141]
    PB.lb_C_work = [-5.12, -512.0, -0.036955815936, -3.210139775093711, -377.50916467454414, -0.0336817159063058, -5.12, -512.0, -0.0512, -5.119999999999999]
    PB.ub_C_work = [3.7836375, 323.3219001168522, 0.05120000000000001, 5.12, 512.0, 0.0512, 3.68855905460902, 339.3660981729729, 0.037162186241472835, 3.2362623805814095]
    PB.ProblemData = _problem_data()
    PB.func_idx = _objective
    PB.func_F = lambda xr: _objective(_cell2idx(xr))
    PB.idx2cell = _idx2cell
    PB.cell2idx = _cell2idx
    PB.decode = lambda xr: _decode(_cell2idx(xr))
    return PB
