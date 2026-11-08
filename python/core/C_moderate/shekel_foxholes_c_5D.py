"""shekel_foxholes_c_5D  DFO Benchmark Suite v2.1.0 core, Python version (self-contained; standard library only).
SHEKEL_FOXHOLES_C_5D  DFO Benchmark Suite v2.1.0 core: shekel_foxholes, dimension 5, pattern C, strategy moderate.

 Copyright (C) 2026 J. F. A. Madeira

 Objective F(z) = f(T(z)), where f is the original function (embedded below, translated from the MATLAB original) and T acts on the finite
 coordinates only: uniform grid of 0 levels on the interval of the domain, piecewise-affine map with |T - Id| <=
 half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
 Global optimum of F on this domain: f* = -10.40561723899244 (L: literature value at a refined point); the optimal solution is listed in
 reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
 Variables: x1:C x2:C x3:C x4:C x5:C


Python version of the MATLAB wrapper shekel_foxholes_c_5D.m: same instance, same values (level indices 1..m as in MATLAB).
Usage:  from shekel_foxholes_c_5D import shekel_foxholes_c_5D;  PB = shekel_foxholes_c_5D();  f = PB.func_idx(PB.x0_idx)
"""
# Copyright (C) 2026 J. F. A. Madeira. LGPL-3.0-or-later.
import math
from types import SimpleNamespace

_NAME = 'shekel_foxholes_c_5D'
_N = 5
_TYPES = ['C', 'C', 'C', 'C', 'C']
_CARD = [math.inf, math.inf, math.inf, math.inf, math.inf]
_LB = [0.0, 1.3797337804196517, 0.0, 1.4936693044564764, 1.9497507000887921]
_UB = [8.47620849609375, 10.0, 8.6873664, 10.0, 10.0]
_LBW = [0.0, 137.9733780419652, 0.0, 1.4936693044564766, 194.9750700088792]
_UBW = [8.47620849609375, 1000.0, 0.086873664, 10.0, 1000.0]
_TB = [math.nan, math.nan, math.nan, math.nan, math.nan]
_TI = [math.nan, math.nan, math.nan, math.nan, math.nan]
_CW = [8.024914889028711, 915.1725753063795, 0.05113927807781459, 7.620860957619284, 456.4088392837862]
_CX = [8.024914889028711, 9.151725753063793, 5.113927807781459, 7.620860957619284, 4.564088392837863]
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

_DATA = [[0.806,9.681,0.667,4.783,9.095,3.517,9.325,6.544,0.211,5.122,2.020],[0.517,9.400,2.041,3.788,7.931,2.882,2.672,3.568,1.284,7.033,7.374],
[0.100,8.025,9.152,5.114,7.621,4.564,4.711,2.996,6.126,0.734,4.982],[0.908,2.196,0.415,5.649,6.979,9.510,9.166,6.304,6.054,9.377,1.426],
[0.965,8.074,8.777,3.467,1.863,6.708,6.349,4.534,0.276,7.633,1.567],[0.669,7.650,5.658,0.720,2.764,3.278,5.283,7.474,6.274,1.409,8.208],
[0.524,1.256,3.605,8.623,6.905,4.584,8.133,6.071,6.888,4.187,5.448],[0.902,8.314,2.261,4.224,1.781,4.124,0.932,8.129,8.658,1.208,5.762],
[0.531,0.226,8.858,1.420,0.945,1.622,4.698,6.228,9.096,0.972,7.637],[0.876,7.305,2.228,1.242,5.928,9.133,1.826,4.060,5.204,8.713,8.247],
[0.462,0.652,7.027,0.508,4.876,8.807,4.632,5.808,6.937,3.291,7.016],[0.491,2.699,3.516,5.874,4.119,4.461,7.496,8.817,0.690,6.593,9.789],
[0.463,8.327,3.897,2.017,9.570,9.825,1.150,1.395,3.885,6.354,0.109],[0.714,2.132,7.006,7.136,2.641,1.882,5.943,7.273,7.691,2.880,0.564],
[0.352,4.707,5.579,4.080,0.581,9.698,8.542,8.077,8.515,9.231,4.670],[0.869,8.304,7.559,8.567,0.322,7.128,8.392,1.472,8.524,2.277,7.826],
[0.813,8.632,4.409,4.832,5.768,7.050,6.715,1.711,4.323,4.405,4.591],[0.811,4.887,9.112,0.170,8.967,9.693,9.867,7.508,7.770,8.382,6.740],
[0.828,2.440,6.686,4.299,1.007,7.008,1.427,9.398,8.480,9.950,1.675],[0.964,6.306,8.583,6.084,1.138,4.350,3.134,7.853,6.061,7.457,2.258],
[0.789,0.652,2.343,1.370,0.821,1.310,1.063,0.689,8.819,8.833,9.070],[0.360,5.558,1.272,5.756,9.857,2.279,2.764,1.284,1.677,1.244,1.234],
[0.369,3.352,7.549,9.817,9.437,8.687,4.167,2.570,6.540,0.228,0.027],[0.992,8.798,0.880,2.370,0.168,1.701,3.680,1.231,2.390,2.499,0.064],
[0.332,1.460,8.057,1.336,7.217,7.914,3.615,9.981,9.198,5.292,1.224],[0.817,0.432,8.645,8.774,0.249,8.081,7.461,4.416,0.652,4.002,4.644],
[0.632,0.679,2.800,5.523,3.049,2.968,7.225,6.730,4.199,9.614,9.229],[0.883,4.263,1.074,7.286,5.599,8.291,5.200,9.214,8.272,4.398,4.506],
[0.608,9.496,4.830,3.150,8.270,5.079,1.231,5.731,9.494,1.883,9.732],[0.326,4.138,2.562,2.532,9.661,5.611,5.500,6.886,2.341,9.699,6.500]]
def _orig(x):
    d = len(x)
    s = 0.0
    for r in _DATA:
        e = 0.0
        for j in range(d): e += (x[j] - r[1 + j]) ** 2
        s += 1 / (r[0] + e)
    return -s

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

def shekel_foxholes_c_5D():
    PB = SimpleNamespace()
    PB.name = _NAME
    PB.nvar = 5
    PB.var_types = ['C', 'C', 'C', 'C', 'C']
    PB.cardinality = [math.inf, math.inf, math.inf, math.inf, math.inf]
    PB.f_star = -10.405617238992441
    PB.f_star_justification = 'L'
    PB.n_minimizers_documented = 1
    PB.description = 'DFO Benchmark Suite v2.1.0 core (Python): C adaptation of shekel_foxholes (F = f o T).'
    PB.x0_idx = [2.9178354442119603, 416.2499190124145, 0.00183477178368, 2.02509354616057, 612.3058000568683]
    PB.x0_alt_idx = [[5.1031079471111305, 355.9435298219347, 0.06229884192768001, 2.3970905153534354, 360.0928805930836], [0.8318935096263886, 380.77557242977934, 0.01934850244608, 1.5786971831291317, 978.8310499251022], [7.785033291578293, 396.14778928225456, 0.07981257259008, 1.950694152321997, 726.6181304613174], [3.4641535699367525, 335.84140009177474, 0.036195043368960005, 2.3226911215148625, 546.984468397183], [5.649426072835923, 389.05291996572754, 0.01326039607296, 2.6946880907077277, 288.17343462433917], [1.378211635351181, 413.8849625735722, 0.05370877403136001, 1.706239001138114, 906.9116039563578], [6.742062324285508, 340.96547237593313, 0.030774126735359997, 2.0782359703309794, 654.698684492573], [2.3963499605655674, 365.7975149837778, 0.07122250469376, 2.450232939523845, 475.06502242843857], [4.581622463464737, 419.0090348577305, 0.04828785739776, 1.631839607299541, 222.85210296465394]]
    PB.suite_release = '2.1.0'
    PB.base_function = 'shekel_foxholes'
    PB.m = 0
    PB.pattern = 'C'
    PB.strategy = 'moderate'
    PB.scale_factors_C = [1.0, 100.0, 0.01, 1.0, 100.0]
    PB.lb_orig = [0.0, 1.3797337804196517, 0.0, 1.4936693044564764, 1.9497507000887921]
    PB.ub_orig = [8.47620849609375, 10.0, 8.6873664, 10.0, 10.0]
    PB.lb_C_work = [0.0, 137.9733780419652, 0.0, 1.4936693044564766, 194.9750700088792]
    PB.ub_C_work = [8.47620849609375, 1000.0, 0.086873664, 10.0, 1000.0]
    PB.ProblemData = _problem_data()
    PB.func_idx = _objective
    PB.func_F = lambda xr: _objective(_cell2idx(xr))
    PB.idx2cell = _idx2cell
    PB.cell2idx = _cell2idx
    PB.decode = lambda xr: _decode(_cell2idx(xr))
    return PB
