"""rastrigin_cd_10D  DFO Benchmark Suite v2.1.0 core, Python version (self-contained; standard library only).
RASTRIGIN_CD_10D  DFO Benchmark Suite v2.1.0 core: rastrigin, dimension 10, pattern CD, strategy moderate.

 Copyright (C) 2026 J. F. A. Madeira

 Objective F(z) = f(T(z)), where f is the original function (embedded below, translated from the MATLAB original) and T acts on the finite
 coordinates only: uniform grid of 17 levels on the interval of the domain, piecewise-affine map with |T - Id| <=
 half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
 Global optimum of F on this domain: f* = 0 (E: exact by structure); the optimal solution is listed in
 reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
 Variables: x1:C x2:D x3:C x4:D x5:C x6:D x7:C x8:D x9:C x10:D


Python version of the MATLAB wrapper rastrigin_cd_10D.m: same instance, same values (level indices 1..m as in MATLAB).
Usage:  from rastrigin_cd_10D import rastrigin_cd_10D;  PB = rastrigin_cd_10D();  f = PB.func_idx(PB.x0_idx)
"""
# Copyright (C) 2026 J. F. A. Madeira. LGPL-3.0-or-later.
import math
from types import SimpleNamespace

_NAME = 'rastrigin_cd_10D'
_N = 10
_TYPES = ['C', 'D', 'C', 'D', 'C', 'D', 'C', 'D', 'C', 'D']
_CARD = [math.inf, 17.0, math.inf, 17.0, math.inf, 17.0, math.inf, 17.0, math.inf, 17.0]
_LB = [-5.12, -5.12, -3.6955815936, -3.210139775093711, -3.775091646745441, -3.3681715906305802, -5.12, -5.12, -5.12, -5.12]
_UB = [3.7836375, 3.2332190011685213, 5.12, 5.12, 5.12, 5.12, 3.6885590546090206, 3.393660981729729, 3.7162186241472837, 3.23626238058141]
_LBW = [-5.12, -5.12, -369.55815936000005, -3.210139775093711, -0.03775091646745442, -3.3681715906305802, -5.12, -5.12, -512.0, -5.12]
_UBW = [3.7836375, 3.2332190011685213, 512.0, 5.12, 0.05120000000000001, 5.12, 3.68855905460902, 3.393660981729729, 371.6218624147284, 3.23626238058141]
_TB = [math.nan, 0.1007618757303268, math.nan, -0.08633735943356946, math.nan, -0.1851072441441124, math.nan, 0.20103811358108015, math.nan, 0.10266398786338105]
_TI = [math.nan, 0.0, math.nan, 0.0, math.nan, 0.0, math.nan, 0.0, math.nan, 0.0]
_CW = [0.0, math.nan, -5.684341886080802e-14, math.nan, -6.938893903907228e-18, math.nan, 0.0, math.nan, 0.0, math.nan]
_CX = [0.0, math.nan, 0.0, math.nan, 0.0, math.nan, 0.0, math.nan, 0.0, math.nan]
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

def rastrigin_cd_10D():
    PB = SimpleNamespace()
    PB.name = _NAME
    PB.nvar = 10
    PB.var_types = ['C', 'D', 'C', 'D', 'C', 'D', 'C', 'D', 'C', 'D']
    PB.cardinality = [math.inf, 17.0, math.inf, 17.0, math.inf, 17.0, math.inf, 17.0, math.inf, 17.0]
    PB.f_star = 0.0
    PB.f_star_justification = 'E'
    PB.n_minimizers_documented = 1
    PB.description = 'DFO Benchmark Suite v2.1.0 core (Python): CD adaptation of rastrigin (F = f o T).'
    PB.x0_idx = [-4.837413848876953, 14.0, -252.20513718599688, 10.0, 0.015713195609152304, 1.0, -3.255375245920338, 14.0, -283.0876871594293, 15.0]
    PB.x0_alt_idx = [[3.301067303466797, 15.0, 192.10017513144322, 11.0, -0.012154972810778958, 15.0, 3.0538694748549817, 2.0, 335.01913218895186, 2.0], [-1.1855312805175782, 14.0, -74.48301225902082, 12.0, -0.03200352441921919, 13.0, 0.068676959910138, 7.0, 109.52016919086418, 5.0], [1.109937762451172, 14.0, 369.8223000584192, 10.0, 0.02834409208725064, 10.0, -2.400157930828031, 12.0, -115.97879380722327, 7.0], [-3.4288305725097654, 15.0, 96.46874600407034, 11.0, 0.008495540478810404, 8.0, -4.8689928215661995, 1.0, -341.4777568053109, 10.0], [2.2055025329589846, 14.0, -305.52177466408966, 12.0, -0.02010168401489277, 5.0, 1.4402518992091196, 7.0, 314.9747799224551, 13.0], [-2.281096051025391, 14.0, 274.1908709310464, 12.0, 0.048265549303068075, 3.0, -1.0285829915290492, 12.0, 89.47581692436768, 16.0], [0.014372991943359281, 15.0, -127.79964973711367, 11.0, 0.020397380883136827, 17.0, -4.013775506473893, 16.0, -136.02314607372, 3.0], [-4.550480218505859, 14.0, 451.9129958580224, 11.0, 0.0005488292746965717, 13.0, 2.295469214301427, 5.0, -361.52210907180756, 5.0], [2.7532849182128905, 15.0, 49.92247518986227, 12.0, -0.027319339145234676, 11.0, -0.17336567643674172, 11.0, 256.58471027657356, 8.0]]
    PB.suite_release = '2.1.0'
    PB.base_function = 'rastrigin'
    PB.m = 17
    PB.pattern = 'CD'
    PB.strategy = 'moderate'
    PB.scale_factors_C = [1.0, 100.0, 0.01, 1.0, 100.0]
    PB.lb_orig = [-5.12, -5.12, -3.6955815936, -3.210139775093711, -3.775091646745441, -3.3681715906305802, -5.12, -5.12, -5.12, -5.12]
    PB.ub_orig = [3.7836375, 3.2332190011685213, 5.12, 5.12, 5.12, 5.12, 3.6885590546090206, 3.393660981729729, 3.7162186241472837, 3.23626238058141]
    PB.lb_C_work = [-5.12, -5.12, -369.55815936000005, -3.210139775093711, -0.03775091646745442, -3.3681715906305802, -5.12, -5.12, -512.0, -5.12]
    PB.ub_C_work = [3.7836375, 3.2332190011685213, 512.0, 5.12, 0.05120000000000001, 5.12, 3.68855905460902, 3.393660981729729, 371.6218624147284, 3.23626238058141]
    PB.ProblemData = _problem_data()
    PB.func_idx = _objective
    PB.func_F = lambda xr: _objective(_cell2idx(xr))
    PB.idx2cell = _idx2cell
    PB.cell2idx = _cell2idx
    PB.decode = lambda xr: _decode(_cell2idx(xr))
    return PB
