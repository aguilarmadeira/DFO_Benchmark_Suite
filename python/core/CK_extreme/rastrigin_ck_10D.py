"""rastrigin_ck_10D  DFO Benchmark Suite v2.1.0 core, Python version (self-contained; standard library only).
RASTRIGIN_CK_10D  DFO Benchmark Suite v2.1.0 core: rastrigin, dimension 10, pattern CK, strategy extreme.

 Copyright (C) 2026 J. F. A. Madeira

 Objective F(z) = f(T(z)), where f is the original function (embedded below, translated from the MATLAB original) and T acts on the finite
 coordinates only: uniform grid of 17 levels on the interval of the domain, piecewise-affine map with |T - Id| <=
 half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
 Global optimum of F on this domain: f* = 0 (E: exact by structure); the optimal solution is listed in
 reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
 Variables: x1:C x2:K x3:C x4:K x5:C x6:K x7:C x8:K x9:C x10:K


Python version of the MATLAB wrapper rastrigin_ck_10D.m: same instance, same values (level indices 1..m as in MATLAB).
Usage:  from rastrigin_ck_10D import rastrigin_ck_10D;  PB = rastrigin_ck_10D();  f = PB.func_idx(PB.x0_idx)
"""
# Copyright (C) 2026 J. F. A. Madeira. LGPL-3.0-or-later.
import math
from types import SimpleNamespace

_NAME = 'rastrigin_ck_10D'
_N = 10
_TYPES = ['C', 'K', 'C', 'K', 'C', 'K', 'C', 'K', 'C', 'K']
_CARD = [math.inf, 17.0, math.inf, 17.0, math.inf, 17.0, math.inf, 17.0, math.inf, 17.0]
_LB = [-5.12, -5.12, -3.6955815936, -3.210139775093711, -3.775091646745441, -3.3681715906305802, -5.12, -5.12, -5.12, -5.12]
_UB = [3.7836375, 3.2332190011685213, 5.12, 5.12, 5.12, 5.12, 3.6885590546090206, 3.393660981729729, 3.7162186241472837, 3.23626238058141]
_LBW = [-5.12, -5.12, -3.6955815936000005, -3.210139775093711, -377509164.67454416, -3.3681715906305802, -512000000.0, -5.12, -512000000.0, -5.12]
_UBW = [3.7836375, 3.2332190011685213, 5.120000000000001, 5.12, 512000000.00000006, 5.12, 368855905.460902, 3.393660981729729, 371621862.4147284, 3.23626238058141]
_TB = [math.nan, 0.1007618757303268, math.nan, -0.08633735943356946, math.nan, -0.1851072441441124, math.nan, 0.20103811358108015, math.nan, 0.10266398786338105]
_TI = [math.nan, 0.0, math.nan, 0.0, math.nan, 0.0, math.nan, 0.0, math.nan, 0.0]
_CW = [0.0, math.nan, -4.440892098500626e-16, math.nan, -5.960464477539063e-08, math.nan, 0.0, math.nan, 5.960464477539063e-08, math.nan]
_CX = [0.0, math.nan, 0.0, math.nan, 0.0, math.nan, 0.0, math.nan, 0.0, math.nan]
_M = 17
_PERM = {1: [14, 3, 15, 13, 9, 4, 2, 11, 1, 6, 8, 10, 5, 12, 16, 17, 7], 3: [10, 9, 16, 3, 17, 2, 12, 6, 1, 8, 11, 7, 5, 13, 15, 14, 4], 5: [11, 12, 15, 16, 14, 7, 8, 13, 10, 9, 4, 3, 5, 2, 17, 1, 6], 7: [17, 10, 5, 15, 13, 7, 16, 14, 9, 1, 6, 11, 4, 3, 12, 8, 2], 9: [11, 14, 13, 12, 2, 6, 8, 1, 17, 15, 10, 7, 4, 5, 3, 16, 9]}

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

def rastrigin_ck_10D():
    PB = SimpleNamespace()
    PB.name = _NAME
    PB.nvar = 10
    PB.var_types = ['C', 'K', 'C', 'K', 'C', 'K', 'C', 'K', 'C', 'K']
    PB.cardinality = [math.inf, 17.0, math.inf, 17.0, math.inf, 17.0, math.inf, 17.0, math.inf, 17.0]
    PB.f_star = 0.0
    PB.f_star_justification = 'E'
    PB.n_minimizers_documented = 1
    PB.description = 'DFO Benchmark Suite v2.1.0 core (Python): CK adaptation of rastrigin (F = f o T).'
    PB.x0_idx = [-4.837413848876953, 14.0, -2.5220513718599684, 10.0, 157131956.09152305, 1.0, -325537524.59203374, 14.0, -283087687.1594293, 15.0]
    PB.x0_alt_idx = [[3.301067303466797, 15.0, 1.9210017513144337, 11.0, -121549728.10778958, 15.0, 305386947.4854981, 2.0, 335019132.18895185, 2.0], [-1.1855312805175782, 14.0, -0.7448301225902081, 12.0, -320035244.19219196, 13.0, 6867695.991013765, 7.0, 109520169.1908642, 5.0], [1.109937762451172, 14.0, 3.6982230005841923, 10.0, 283440920.8725063, 10.0, -240015793.08280313, 12.0, -115978793.80722332, 7.0], [-3.4288305725097654, 15.0, 0.9646874600407043, 11.0, 84955404.788104, 8.0, -486899282.1566199, 1.0, -341477756.80531085, 10.0], [2.2055025329589846, 14.0, -3.0552177466408965, 12.0, -201016840.14892766, 5.0, 144025189.9209119, 7.0, 314974779.9224551, 13.0], [-2.281096051025391, 14.0, 2.7419087093104655, 12.0, 482655493.0306807, 3.0, -102858299.15290493, 12.0, 89475816.92436767, 16.0], [0.014372991943359281, 15.0, -1.2779964973711362, 11.0, 203973808.83136827, 17.0, -401377550.6473893, 16.0, -136023146.07372004, 3.0], [-4.550480218505859, 14.0, 4.5191299585802245, 11.0, 5488292.746965766, 13.0, 229546921.43014264, 5.0, -361522109.07180756, 5.0], [2.7532849182128905, 15.0, 0.4992247518986237, 12.0, -273193391.45234674, 11.0, -17336567.643674135, 11.0, 256584710.27657354, 8.0]]
    PB.suite_release = '2.1.0'
    PB.base_function = 'rastrigin'
    PB.m = 17
    PB.pattern = 'CK'
    PB.strategy = 'extreme'
    PB.scale_factors_C = [1.0, 1.0, 100000000.0, 100000000.0, 100000000.0]
    PB.lb_orig = [-5.12, -5.12, -3.6955815936, -3.210139775093711, -3.775091646745441, -3.3681715906305802, -5.12, -5.12, -5.12, -5.12]
    PB.ub_orig = [3.7836375, 3.2332190011685213, 5.12, 5.12, 5.12, 5.12, 3.6885590546090206, 3.393660981729729, 3.7162186241472837, 3.23626238058141]
    PB.lb_C_work = [-5.12, -5.12, -3.6955815936000005, -3.210139775093711, -377509164.67454416, -3.3681715906305802, -512000000.0, -5.12, -512000000.0, -5.12]
    PB.ub_C_work = [3.7836375, 3.2332190011685213, 5.120000000000001, 5.12, 512000000.00000006, 5.12, 368855905.460902, 3.393660981729729, 371621862.4147284, 3.23626238058141]
    PB.ProblemData = _problem_data()
    PB.func_idx = _objective
    PB.func_F = lambda xr: _objective(_cell2idx(xr))
    PB.idx2cell = _idx2cell
    PB.cell2idx = _cell2idx
    PB.decode = lambda xr: _decode(_cell2idx(xr))
    return PB
