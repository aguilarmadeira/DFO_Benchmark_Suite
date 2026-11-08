"""hartman_4_c_3D  DFO Benchmark Suite v2.1.0 core, Python version (self-contained; standard library only).
HARTMAN_4_C_3D  DFO Benchmark Suite v2.1.0 core: hartman_4, dimension 3, pattern C, strategy moderate.

 Copyright (C) 2026 J. F. A. Madeira

 Objective F(z) = f(T(z)), where f is the original function (embedded below, translated from the MATLAB original) and T acts on the finite
 coordinates only: uniform grid of 0 levels on the interval of the domain, piecewise-affine map with |T - Id| <=
 half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
 Global optimum of F on this domain: f* = -3.862782147819745 (L: literature value at a refined point); the optimal solution is listed in
 reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
 Variables: x1:C x2:C x3:C


Python version of the MATLAB wrapper hartman_4_c_3D.m: same instance, same values (level indices 1..m as in MATLAB).
Usage:  from hartman_4_c_3D import hartman_4_c_3D;  PB = hartman_4_c_3D();  f = PB.func_idx(PB.x0_idx)
"""
# Copyright (C) 2026 J. F. A. Madeira. LGPL-3.0-or-later.
import math
from types import SimpleNamespace

_NAME = 'hartman_4_c_3D'
_N = 3
_TYPES = ['C', 'C', 'C']
_CARD = [math.inf, math.inf, math.inf]
_LB = [0.0, 0.0, 0.0]
_UB = [0.830941162109375, 0.8131514504902708, 0.86929664]
_LBW = [0.0, 0.0, 0.0]
_UBW = [0.830941162109375, 81.31514504902708, 0.0086929664]
_TB = [math.nan, math.nan, math.nan]
_TI = [math.nan, math.nan, math.nan]
_CW = [0.11461399999999998, 55.56489999999999, 0.00852547]
_CX = [0.114614, 0.555649, 0.852547]
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

_A3 = [[3, 10, 30, 1, 0.3689, 0.1170, 0.2673],
       [0.1, 10, 35, 1.2, 0.4699, 0.4387, 0.7470],
       [3, 10, 30, 3, 0.1091, 0.8732, 0.5547],
       [0.1, 10, 35, 3.2, 0.03815, 0.5743, 0.8828]]
_A6 = [[10, 3, 17, 3.5, 1.7, 8, 1, 0.1312, 0.1696, 0.5569, 0.0124, 0.8283, 0.5886],
       [0.05, 10, 17, 0.1, 8, 14, 1.2, 0.2329, 0.4135, 0.8307, 0.3736, 0.1004, 0.9991],
       [3, 3.5, 1.7, 10, 17, 8, 3, 0.2348, 0.1451, 0.3522, 0.2883, 0.3047, 0.6650],
       [17, 8, 0.05, 10, 0.1, 14, 3.2, 0.4047, 0.8828, 0.8732, 0.5743, 0.1091, 0.0381]]
def _orig(x):
    n = len(x)
    A = _A3 if n == 3 else _A6
    f = 0.0
    for r in A:
        e = 0.0
        for j in range(n): e += r[j] * (x[j] - r[n + 1 + j]) ** 2
        f += r[n] * math.exp(-e)
    return -f

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

def hartman_4_c_3D():
    PB = SimpleNamespace()
    PB.name = _NAME
    PB.nvar = 3
    PB.var_types = ['C', 'C', 'C']
    PB.cardinality = [math.inf, math.inf, math.inf]
    PB.f_star = -3.862782147819745
    PB.f_star_justification = 'L'
    PB.n_minimizers_documented = 1
    PB.description = 'DFO Benchmark Suite v2.1.0 core (Python): C adaptation of hartman_4 (F = f o T).'
    PB.x0_idx = [0.4158763140439987, 79.45608823492037, 0.00011405171916800001]
    PB.x0_alt_idx = [[0.006897460818290711, 73.76737438375388, 0.006164356333568001], [0.8248551672697066, 76.1097859695283, 0.001866553745408], [0.40613872230052944, 81.1292393676164, 0.007916858359808], [0.6203657406568527, 74.54817824567868, 0.003552293789696], [0.19678049981594084, 76.89058983145313, 0.001257350660096], [0.7226104539632797, 79.23300141722757, 0.005304795815936001], [0.30389400899410246, 72.35449120503279, 0.003009852686336], [0.5181210273504256, 77.37394460312088, 0.007057297842176], [0.09210138857364654, 79.71635618889532, 0.004762354712576]]
    PB.suite_release = '2.1.0'
    PB.base_function = 'hartman_4'
    PB.m = 0
    PB.pattern = 'C'
    PB.strategy = 'moderate'
    PB.scale_factors_C = [1.0, 100.0, 0.01]
    PB.lb_orig = [0.0, 0.0, 0.0]
    PB.ub_orig = [0.830941162109375, 0.8131514504902708, 0.86929664]
    PB.lb_C_work = [0.0, 0.0, 0.0]
    PB.ub_C_work = [0.830941162109375, 81.31514504902708, 0.0086929664]
    PB.ProblemData = _problem_data()
    PB.func_idx = _objective
    PB.func_F = lambda xr: _objective(_cell2idx(xr))
    PB.idx2cell = _idx2cell
    PB.cell2idx = _cell2idx
    PB.decode = lambda xr: _decode(_cell2idx(xr))
    return PB
