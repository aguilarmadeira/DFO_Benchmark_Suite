"""griewank_cdk_10D  DFO Benchmark Suite v2.1.0 core, Python version (self-contained; standard library only).
GRIEWANK_CDK_10D  DFO Benchmark Suite v2.1.0 core: griewank, dimension 10, pattern CDK, strategy extreme.

 Copyright (C) 2026 J. F. A. Madeira

 Objective F(z) = f(T(z)), where f is the original function (embedded below, translated from the MATLAB original) and T acts on the finite
 coordinates only: uniform grid of 17 levels on the interval of the domain, piecewise-affine map with |T - Id| <=
 half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
 Global optimum of F on this domain: f* = 0 (E: exact by structure); the optimal solution is listed in
 reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
 Variables: x1:C x2:K x3:D x4:C x5:D x6:K x7:C x8:K x9:D x10:C


Python version of the MATLAB wrapper griewank_cdk_10D.m: same instance, same values (level indices 1..m as in MATLAB).
Usage:  from griewank_cdk_10D import griewank_cdk_10D;  PB = griewank_cdk_10D();  f = PB.func_idx(PB.x0_idx)
"""
# Copyright (C) 2026 J. F. A. Madeira. LGPL-3.0-or-later.
import math
from types import SimpleNamespace

_NAME = 'griewank_cdk_10D'
_N = 10
_TYPES = ['C', 'K', 'D', 'C', 'D', 'K', 'C', 'K', 'D', 'C']
_CARD = [math.inf, 17.0, 17.0, math.inf, 17.0, 17.0, math.inf, 17.0, 17.0, math.inf]
_LB = [-270.0029296875, -256.7433826144388, -400.0, -400.0, -400.0, -400.0, -400.0, -400.0, -278.2615270814498, -287.5916191725778]
_UB = [400.0, 400.0, 264.525312, 285.07788421491045, 275.49101837306193, 285.0082280032212, 291.6565654146861, 262.337520430322, 400.0, 400.0]
_LBW = [-270.0029296875, -256.7433826144388, -400.0, -400.0, -400.0, -400.0, -40000000000.0, -400.0, -278.2615270814498, -28759161917.257782]
_UBW = [400.0, 400.0, 264.525312, 285.07788421491045, 275.49101837306193, 285.0082280032212, 29165656541.468613, 262.337520430322, 400.0, 40000000000.0]
_TB = [math.nan, -10.46461413402426, 15.328319999999962, math.nan, -20.036302165152676, -14.68287174818812, math.nan, 13.960950268951251, 18.47789101668451, math.nan]
_TI = [math.nan, 0.0, 0.0, math.nan, 0.0, 0.0, math.nan, 0.0, 0.0, math.nan]
_CW = [0.0, math.nan, math.nan, -5.684341886080802e-14, math.nan, math.nan, 0.0, math.nan, math.nan, 0.0]
_CX = [0.0, math.nan, math.nan, 0.0, math.nan, math.nan, 0.0, math.nan, math.nan, 0.0]
_M = 17
_PERM = {1: [15, 14, 1, 9, 13, 3, 12, 10, 2, 7, 11, 16, 8, 5, 17, 6, 4], 5: [13, 5, 9, 12, 6, 4, 16, 3, 15, 14, 7, 11, 8, 10, 17, 2, 1], 7: [3, 5, 6, 4, 12, 14, 7, 15, 1, 2, 10, 13, 16, 11, 9, 17, 8]}

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
    n = len(x)
    s = 0.0
    for v in x: s += v ** 2
    p = 1.0
    for i in range(n): p *= math.cos(x[i] / math.sqrt(i + 1))
    return 1 + 1 / 4000 * s - p

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

def griewank_cdk_10D():
    PB = SimpleNamespace()
    PB.name = _NAME
    PB.nvar = 10
    PB.var_types = ['C', 'K', 'D', 'C', 'D', 'K', 'C', 'K', 'D', 'C']
    PB.cardinality = [math.inf, 17.0, 17.0, math.inf, 17.0, 17.0, math.inf, 17.0, 17.0, math.inf]
    PB.f_star = 0.0
    PB.f_star_justification = 'E'
    PB.n_minimizers_documented = 1
    PB.description = 'DFO Benchmark Suite v2.1.0 core (Python): CDK adaptation of griewank (F = f o T).'
    PB.x0_idx = [64.34423542022705, 15.0, 10.0, -303.55838198057484, 16.0, 11.0, -29666681884.502754, 14.0, 9.0, -11668750685.045858]
    PB.x0_alt_idx = [[222.3576216697693, 14.0, 3.0, -383.1655163812246, 13.0, 9.0, 19874117091.566452, 3.0, 4.0, 349794858.16091156], [-107.40944528579712, 14.0, 13.0, -353.20584214442096, 8.0, 7.0, 488587057.45241547, 8.0, 16.0, 12368340401.367683], [308.7251868247986, 15.0, 6.0, -323.2461679076173, 3.0, 4.0, -22951432918.437107, 13.0, 12.0, 24386885944.574463], [-21.041880130767822, 14.0, 16.0, -389.15745122858533, 15.0, 2.0, 26589366057.63211, 2.0, 8.0, 36405431487.78123], [139.91585493087769, 14.0, 9.0, -359.19777699178167, 10.0, 16.0, 7203836023.518066, 7.0, 3.0, -20253426073.05069], [-189.85121202468872, 15.0, 1.0, -342.9339538346597, 6.0, 12.0, -12181694010.595978, 13.0, 15.0, -8234880529.843922], [351.90896940231323, 13.0, 12.0, -312.974279597856, 1.0, 10.0, -35621713986.4855, 1.0, 11.0, 3783665013.36285], [22.141902446746826, 14.0, 4.0, -378.88556291882406, 13.0, 7.0, 13919084989.58371, 6.0, 7.0, 15802210556.569626], [181.1367383003235, 15.0, 16.0, -348.9258886820204, 8.0, 5.0, -5466445044.530327, 11.0, 2.0, 27820756099.776398]]
    PB.suite_release = '2.1.0'
    PB.base_function = 'griewank'
    PB.m = 17
    PB.pattern = 'CDK'
    PB.strategy = 'extreme'
    PB.scale_factors_C = [1.0, 1.0, 100000000.0, 100000000.0]
    PB.lb_orig = [-270.0029296875, -256.7433826144388, -400.0, -400.0, -400.0, -400.0, -400.0, -400.0, -278.2615270814498, -287.5916191725778]
    PB.ub_orig = [400.0, 400.0, 264.525312, 285.07788421491045, 275.49101837306193, 285.0082280032212, 291.6565654146861, 262.337520430322, 400.0, 400.0]
    PB.lb_C_work = [-270.0029296875, -256.7433826144388, -400.0, -400.0, -400.0, -400.0, -40000000000.0, -400.0, -278.2615270814498, -28759161917.257782]
    PB.ub_C_work = [400.0, 400.0, 264.525312, 285.07788421491045, 275.49101837306193, 285.0082280032212, 29165656541.468613, 262.337520430322, 400.0, 40000000000.0]
    PB.ProblemData = _problem_data()
    PB.func_idx = _objective
    PB.func_F = lambda xr: _objective(_cell2idx(xr))
    PB.idx2cell = _idx2cell
    PB.cell2idx = _cell2idx
    PB.decode = lambda xr: _decode(_cell2idx(xr))
    return PB
