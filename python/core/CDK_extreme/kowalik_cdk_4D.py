"""kowalik_cdk_4D  DFO Benchmark Suite v2.1.0 core, Python version (self-contained; standard library only).
KOWALIK_CDK_4D  DFO Benchmark Suite v2.1.0 core: kowalik, dimension 4, pattern CDK, strategy extreme.

 Copyright (C) 2026 J. F. A. Madeira

 Objective F(z) = f(T(z)), where f is the original function (embedded below, translated from the MATLAB original) and T acts on the finite
 coordinates only: uniform grid of 151 levels on the interval of the domain, piecewise-affine map with |T - Id| <=
 half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
 Global optimum of F on this domain: f* = 0.0003074859878056062 (L: literature value at a refined point); the optimal solution is listed in
 reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
 Variables: x1:C x2:K x3:D x4:C


Python version of the MATLAB wrapper kowalik_cdk_4D.m: same instance, same values (level indices 1..m as in MATLAB).
Usage:  from kowalik_cdk_4D import kowalik_cdk_4D;  PB = kowalik_cdk_4D();  f = PB.func_idx(PB.x0_idx)
"""
# Copyright (C) 2026 J. F. A. Madeira. LGPL-3.0-or-later.
import math
from types import SimpleNamespace

_NAME = 'kowalik_cdk_4D'
_N = 4
_TYPES = ['C', 'K', 'D', 'C']
_CARD = [math.inf, 151.0, 151.0, math.inf]
_LB = [0.059979711914062504, 0.08174305746075294, 0.0, 0.05673411078717201]
_UB = [0.42, 0.42, 0.3474645888, 0.42]
_LBW = [0.059979711914062483, 0.08174305746075294, 0.0, 5673411.078717202]
_UBW = [0.42, 0.42, 0.3474645888, 42000000.0]
_TB = [math.nan, 0.189985279073312, 0.122770821376, math.nan]
_TI = [math.nan, 0.1908362411132952, 0.1231173000994744, math.nan]
_CW = [0.19283345290075488, math.nan, math.nan, 13576599.04379078]
_CX = [0.1928334529007549, math.nan, math.nan, 0.1357659904379078]
_M = 151
_PERM = {1: [14, 97, 127, 96, 88, 54, 51, 125, 18, 5, 91, 126, 21, 114, 52, 131, 86, 4, 70, 150, 28, 105, 149, 46, 23, 61, 6, 72, 13, 98, 121, 24, 48, 77, 136, 44, 64, 20, 56, 117, 69, 49, 40, 3, 120, 137, 111, 109, 74, 142, 8, 65, 146, 47, 12, 95, 116, 62, 106, 103, 35, 118, 92, 68, 76, 31, 144, 60, 66, 27, 134, 84, 101, 90, 138, 16, 85, 122, 7, 45, 75, 15, 25, 102, 145, 36, 63, 78, 19, 37, 32, 50, 141, 53, 139, 130, 99, 112, 108, 38, 135, 81, 10, 73, 71, 94, 129, 11, 57, 148, 26, 82, 39, 151, 34, 29, 110, 147, 58, 9, 113, 59, 128, 41, 17, 79, 100, 123, 30, 124, 1, 132, 143, 115, 140, 2, 89, 133, 104, 22, 87, 119, 67, 42, 33, 43, 80, 83, 107, 93, 55]}

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

_a = [0.1957, 0.1947, 0.1735, 0.16, 0.0844, 0.0627, 0.0456, 0.0342, 0.0323, 0.0235, 0.0246]
_b = [0.25, 0.5, 1, 2, 4, 6, 8, 10, 12, 14, 16]
def _orig(x):
    s = 0.0
    for a, b in zip(_a, _b): s += (a - x[0] * (1 + x[1] * b) / (1 + x[2] * b + x[3] * b ** 2)) ** 2
    return s

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

def kowalik_cdk_4D():
    PB = SimpleNamespace()
    PB.name = _NAME
    PB.nvar = 4
    PB.var_types = ['C', 'K', 'D', 'C']
    PB.cardinality = [math.inf, 151.0, 151.0, math.inf]
    PB.f_star = 0.00030748598780560617
    PB.f_star_justification = 'L'
    PB.n_minimizers_documented = 1
    PB.description = 'DFO Benchmark Suite v2.1.0 core (Python): CDK adaptation of kowalik (F = f o T).'
    PB.x0_idx = [0.10515803908109664, 105.0, 93.0, 41924351.126777835]
    PB.x0_alt_idx = [[0.33298337763547897, 109.0, 53.0, 37703144.00098112], [0.1515669043421745, 114.0, 123.0, 39291770.33864655], [0.24438463486433026, 101.0, 84.0, 40880396.676311985], [0.06085866769552229, 111.0, 8.0, 37385418.73344804], [0.415253638780117, 115.0, 113.0, 38974045.07111347], [0.23805615323781965, 105.0, 39.0, 39836442.22584613], [0.32665489600896835, 107.0, 143.0, 41425068.56351156], [0.149457410466671, 116.0, 69.0, 37930090.62064761], [0.3698995204567909, 105.0, 29.0, 39518716.95831304]]
    PB.suite_release = '2.1.0'
    PB.base_function = 'kowalik'
    PB.m = 151
    PB.pattern = 'CDK'
    PB.strategy = 'extreme'
    PB.scale_factors_C = [1.0, 100000000.0]
    PB.lb_orig = [0.059979711914062504, 0.08174305746075294, 0.0, 0.05673411078717201]
    PB.ub_orig = [0.42, 0.42, 0.3474645888, 0.42]
    PB.lb_C_work = [0.059979711914062483, 0.08174305746075294, 0.0, 5673411.078717202]
    PB.ub_C_work = [0.42, 0.42, 0.3474645888, 42000000.0]
    PB.ProblemData = _problem_data()
    PB.func_idx = _objective
    PB.func_F = lambda xr: _objective(_cell2idx(xr))
    PB.idx2cell = _idx2cell
    PB.cell2idx = _cell2idx
    PB.decode = lambda xr: _decode(_cell2idx(xr))
    return PB
