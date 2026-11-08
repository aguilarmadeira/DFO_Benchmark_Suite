"""shekel_foxholes_cdk_10D  DFO Benchmark Suite v2.1.0 core, Python version (self-contained; standard library only).
SHEKEL_FOXHOLES_CDK_10D  DFO Benchmark Suite v2.1.0 core: shekel_foxholes, dimension 10, pattern CDK, strategy extreme.

 Copyright (C) 2026 J. F. A. Madeira

 Objective F(z) = f(T(z)), where f is the original function (embedded below, translated from the MATLAB original) and T acts on the finite
 coordinates only: uniform grid of 17 levels on the interval of the domain, piecewise-affine map with |T - Id| <=
 half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
 Global optimum of F on this domain: f* = -10.20879279215383 (L: literature value at a refined point); the optimal solution is listed in
 reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
 Variables: x1:C x2:K x3:D x4:C x5:D x6:K x7:C x8:K x9:D x10:C


Python version of the MATLAB wrapper shekel_foxholes_cdk_10D.m: same instance, same values (level indices 1..m as in MATLAB).
Usage:  from shekel_foxholes_cdk_10D import shekel_foxholes_cdk_10D;  PB = shekel_foxholes_cdk_10D();  f = PB.func_idx(PB.x0_idx)
"""
# Copyright (C) 2026 J. F. A. Madeira. LGPL-3.0-or-later.
import math
from types import SimpleNamespace

_NAME = 'shekel_foxholes_cdk_10D'
_N = 10
_TYPES = ['C', 'K', 'D', 'C', 'D', 'K', 'C', 'K', 'D', 'C']
_CARD = [math.inf, 17.0, 17.0, math.inf, 17.0, 17.0, math.inf, 17.0, 17.0, math.inf]
_LB = [0.0, 1.6130671137529848, 1.4526336, 1.5936693044564765, 1.3191721876920974, 1.6072196351668362, 1.8508327247039666, 0.0, 0.0, 0.0]
_UB = [8.12620849609375, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 8.606919836403957, 8.29433714144818, 8.594062897207758]
_LBW = [0.0, 1.6130671137529848, 1.4526336, 1.5936693044564763, 1.3191721876920974, 1.6072196351668362, 185083272.4703967, 0.0, 0.0, 0.0]
_UBW = [8.12620849609375, 10.0, 10.0, 10.0, 10.0, 10.0, 1000000000.0, 8.606919836403957, 8.29433714144818, 859406289.7207758]
_TB = [math.nan, 8.951633389219122, 5.1921064, math.nan, 4.574482617307561, 4.754512271979273, math.nan, 5.91725738752772, 0.5183960713405112, math.nan]
_TI = [math.nan, 9.151928383187155, 5.113991012282246, math.nan, 4.564022051654048, 4.711005066398696, math.nan, 6.125993024510581, 0.7340577389717022, math.nan]
_CW = [8.024966814917121, math.nan, math.nan, 7.620959243851917, math.nan, math.nan, 299603060.20294213, math.nan, math.nan, 498199937.96057355]
_CX = [8.024966814917121, math.nan, math.nan, 7.620959243851917, math.nan, math.nan, 2.996030602029421, math.nan, math.nan, 4.981999379605735]
_M = 17
_PERM = {1: [2, 6, 17, 16, 1, 14, 9, 7, 11, 8, 13, 12, 5, 3, 15, 4, 10], 5: [8, 17, 16, 6, 2, 9, 10, 3, 12, 15, 1, 4, 5, 7, 11, 14, 13], 7: [17, 8, 13, 7, 9, 15, 6, 1, 11, 3, 4, 12, 5, 10, 16, 2, 14]}

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

def shekel_foxholes_cdk_10D():
    PB = SimpleNamespace()
    PB.name = _NAME
    PB.nvar = 10
    PB.var_types = ['C', 'K', 'D', 'C', 'D', 'K', 'C', 'K', 'D', 'C']
    PB.cardinality = [math.inf, 17.0, 17.0, math.inf, 17.0, 17.0, math.inf, 17.0, 17.0, math.inf]
    PB.f_star = -10.208792792153835
    PB.f_star_justification = 'L'
    PB.n_minimizers_documented = 1
    PB.description = 'DFO Benchmark Suite v2.1.0 core (Python): CDK adaptation of shekel_foxholes (F = f o T).'
    PB.x0_idx = [6.8604562938213345, 11.0, 5.0, 3.3197505343202653, 11.0, 15.0, 884886584.794312, 7.0, 13.0, 214631338.3365142]
    PB.x0_alt_idx = [[2.860838049650192, 10.0, 16.0, 3.687374325670565, 6.0, 13.0, 656484318.6008939, 11.0, 8.0, 364848608.95358205], [4.8606471717357635, 11.0, 8.0, 2.8786019846999054, 2.0, 11.0, 380311643.79185873, 16.0, 4.0, 515065879.5706499], [0.8610289275646209, 11.0, 4.0, 3.246225776050205, 13.0, 8.0, 964006324.0639274, 6.0, 16.0, 665283150.1877177], [7.384215825796127, 10.0, 11.0, 3.6138495674005044, 9.0, 6.0, 735604057.8705094, 11.0, 12.0, 815500420.8047856], [3.3845975816249845, 11.0, 7.0, 3.9814733587508044, 4.0, 2.0, 507201791.677091, 16.0, 8.0, 107333287.89575149], [5.384406703710556, 11.0, 14.0, 3.0046444274485795, 17.0, 16.0, 278799525.48367286, 5.0, 4.0, 257550558.5128193], [1.3847884595394133, 10.0, 10.0, 3.3722682187988795, 12.0, 14.0, 814723797.1401248, 10.0, 15.0, 407767829.12988716], [6.360504013299941, 11.0, 2.0, 3.7398920101491795, 8.0, 11.0, 586321530.9467065, 15.0, 11.0, 557985099.746955], [2.3608857691287994, 12.0, 13.0, 2.9311196691785195, 3.0, 9.0, 357919264.75328827, 4.0, 7.0, 708202370.3640229]]
    PB.suite_release = '2.1.0'
    PB.base_function = 'shekel_foxholes'
    PB.m = 17
    PB.pattern = 'CDK'
    PB.strategy = 'extreme'
    PB.scale_factors_C = [1.0, 1.0, 100000000.0, 100000000.0]
    PB.lb_orig = [0.0, 1.6130671137529848, 1.4526336, 1.5936693044564765, 1.3191721876920974, 1.6072196351668362, 1.8508327247039666, 0.0, 0.0, 0.0]
    PB.ub_orig = [8.12620849609375, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 8.606919836403957, 8.29433714144818, 8.594062897207758]
    PB.lb_C_work = [0.0, 1.6130671137529848, 1.4526336, 1.5936693044564763, 1.3191721876920974, 1.6072196351668362, 185083272.4703967, 0.0, 0.0, 0.0]
    PB.ub_C_work = [8.12620849609375, 10.0, 10.0, 10.0, 10.0, 10.0, 1000000000.0, 8.606919836403957, 8.29433714144818, 859406289.7207758]
    PB.ProblemData = _problem_data()
    PB.func_idx = _objective
    PB.func_F = lambda xr: _objective(_cell2idx(xr))
    PB.idx2cell = _idx2cell
    PB.cell2idx = _cell2idx
    PB.decode = lambda xr: _decode(_cell2idx(xr))
    return PB
