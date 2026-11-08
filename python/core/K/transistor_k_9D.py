"""transistor_k_9D  DFO Benchmark Suite v2.1.0 core, Python version (self-contained; standard library only).
TRANSISTOR_K_9D  DFO Benchmark Suite v2.1.0 core: transistor, dimension 9, pattern K, strategy baseline.

 Copyright (C) 2026 J. F. A. Madeira

 Objective F(z) = f(T(z)), where f is the original function (embedded below, translated from the MATLAB original) and T acts on the finite
 coordinates only: uniform grid of 9 levels on the interval of the domain, piecewise-affine map with |T - Id| <=
 half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
 Global optimum of F on this domain: f* = 1.60049326654736e-24 (E: exact by structure); the optimal solution is listed in
 reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
 Variables: x1:K x2:K x3:K x4:K x5:K x6:K x7:K x8:K x9:K


Python version of the MATLAB wrapper transistor_k_9D.m: same instance, same values (level indices 1..m as in MATLAB).
Usage:  from transistor_k_9D import transistor_k_9D;  PB = transistor_k_9D();  f = PB.func_idx(PB.x0_idx)
"""
# Copyright (C) 2026 J. F. A. Madeira. LGPL-3.0-or-later.
import math
from types import SimpleNamespace

_NAME = 'transistor_k_9D'
_N = 9
_TYPES = ['K', 'K', 'K', 'K', 'K', 'K', 'K', 'K', 'K']
_CARD = [9.0, 9.0, 9.0, 9.0, 9.0, 9.0, 9.0, 9.0, 9.0]
_LB = [-10.0, -10.0, -10.0, -10.0, -7.350085376681921, -6.992661321382304, -10.0, -10.0, -10.0]
_UB = [6.8211669921875, 7.03312503175329, 6.7027328, 6.180008329862558, 10.0, 10.0, 6.7876079069934505, 6.403313357018439, 7.316462562669516]
_LBW = [-10.0, -10.0, -10.0, -10.0, -7.350085376681921, -6.992661321382304, -10.0, -10.0, -10.0]
_UBW = [6.8211669921875, 7.03312503175329, 6.7027328, 6.180008329862558, 10.0, 10.0, 6.7876079069934505, 6.403313357018439, 7.316462562669516]
_TB = [0.5132293701171875, 0.6457031448458057, 0.4392080000000007, 2.1350062473969196, 7.83123932791476, 7.8759173348272125, 4.689156918619268, 0.25207084813652614, 2.9873469220021356]
_TI = [0.8999999526168565, 0.4499874719815851, 1.000006482465387, 2.000068541623683, 7.999971440507853, 7.999692684218543, 5.000031275929672, 0.9999877234567898, 2.00005248348619]
_CW = [math.nan, math.nan, math.nan, math.nan, math.nan, math.nan, math.nan, math.nan, math.nan]
_CX = [math.nan, math.nan, math.nan, math.nan, math.nan, math.nan, math.nan, math.nan, math.nan]
_M = 9
_PERM = {0: [3, 4, 9, 6, 2, 7, 8, 1, 5], 1: [7, 5, 1, 8, 9, 4, 3, 6, 2], 2: [7, 1, 3, 8, 6, 2, 5, 4, 9], 3: [4, 1, 6, 8, 2, 7, 5, 9, 3], 4: [6, 4, 1, 5, 2, 3, 9, 8, 7], 5: [8, 7, 3, 5, 2, 6, 9, 4, 1], 6: [4, 6, 5, 1, 9, 7, 3, 2, 8], 7: [9, 3, 8, 6, 2, 1, 4, 7, 5], 8: [1, 7, 3, 4, 2, 8, 9, 5, 6]}

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

_G = [[0.485, 0.369, 5.2095, 23.3037, 28.5132], [0.752, 1.254, 10.0677, 101.779, 111.8467],
      [0.869, 0.703, 22.9274, 111.461, 134.3884], [0.982, 1.455, 20.2153, 191.267, 211.4823]]
def _orig(x):
    s = 0.0
    for g in _G:
        alfa = (1 - x[0] * x[1]) * x[2] * (math.exp(x[4] * (g[0] - g[2] * x[6] * 0.001 - g[4] * x[7] * 0.001)) - 1) - g[4] + g[3] * x[1]
        beta = (1 - x[0] * x[1]) * x[3] * (math.exp(x[5] * (g[0] - g[1] - g[2] * x[6] * 0.001 + g[3] * x[8] * 0.001)) - 1) - g[4] * x[0] + g[3]
        s += alfa ** 2 + beta ** 2
    return (x[0] * x[2] - x[1] * x[3]) ** 2 + s

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

def transistor_k_9D():
    PB = SimpleNamespace()
    PB.name = _NAME
    PB.nvar = 9
    PB.var_types = ['K', 'K', 'K', 'K', 'K', 'K', 'K', 'K', 'K']
    PB.cardinality = [9.0, 9.0, 9.0, 9.0, 9.0, 9.0, 9.0, 9.0, 9.0]
    PB.f_star = 1.6004932665473596e-24
    PB.f_star_justification = 'E'
    PB.n_minimizers_documented = 1
    PB.description = 'DFO Benchmark Suite v2.1.0 core (Python): K adaptation of transistor (F = f o T).'
    PB.x0_idx = [4.0, 1.0, 5.0, 7.0, 6.0, 7.0, 5.0, 9.0, 2.0]
    PB.x0_alt_idx = [[6.0, 2.0, 1.0, 8.0, 3.0, 6.0, 3.0, 3.0, 8.0], [2.0, 1.0, 7.0, 7.0, 2.0, 4.0, 8.0, 5.0, 6.0], [7.0, 1.0, 3.0, 7.0, 7.0, 3.0, 6.0, 8.0, 4.0], [3.0, 2.0, 8.0, 8.0, 5.0, 2.0, 4.0, 3.0, 2.0], [5.0, 1.0, 4.0, 7.0, 3.0, 8.0, 1.0, 5.0, 8.0], [1.0, 2.0, 2.0, 7.0, 9.0, 7.0, 7.0, 8.0, 6.0], [9.0, 2.0, 6.0, 7.0, 6.0, 6.0, 4.0, 2.0, 4.0], [5.0, 1.0, 4.0, 8.0, 5.0, 5.0, 2.0, 5.0, 2.0], [7.0, 1.0, 7.0, 7.0, 2.0, 4.0, 8.0, 7.0, 8.0]]
    PB.suite_release = '2.1.0'
    PB.base_function = 'transistor'
    PB.m = 9
    PB.pattern = 'K'
    PB.strategy = 'baseline'
    PB.scale_factors_C = []
    PB.lb_orig = [-10.0, -10.0, -10.0, -10.0, -7.350085376681921, -6.992661321382304, -10.0, -10.0, -10.0]
    PB.ub_orig = [6.8211669921875, 7.03312503175329, 6.7027328, 6.180008329862558, 10.0, 10.0, 6.7876079069934505, 6.403313357018439, 7.316462562669516]
    PB.lb_C_work = [-10.0, -10.0, -10.0, -10.0, -7.350085376681921, -6.992661321382304, -10.0, -10.0, -10.0]
    PB.ub_C_work = [6.8211669921875, 7.03312503175329, 6.7027328, 6.180008329862558, 10.0, 10.0, 6.7876079069934505, 6.403313357018439, 7.316462562669516]
    PB.ProblemData = _problem_data()
    PB.func_idx = _objective
    PB.func_F = lambda xr: _objective(_cell2idx(xr))
    PB.idx2cell = _idx2cell
    PB.cell2idx = _cell2idx
    PB.decode = lambda xr: _decode(_cell2idx(xr))
    return PB
