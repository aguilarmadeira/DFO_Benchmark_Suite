"""schwefel_cdk_10D  DFO Benchmark Suite v2.1.0 core, Python version (self-contained; standard library only).
SCHWEFEL_CDK_10D  DFO Benchmark Suite v2.1.0 core: schwefel, dimension 10, pattern CDK, strategy moderate.

 Copyright (C) 2026 J. F. A. Madeira

 Objective F(z) = f(T(z)), where f is the original function (embedded below, translated from the MATLAB original) and T acts on the finite
 coordinates only: uniform grid of 17 levels on the interval of the domain, piecewise-affine map with |T - Id| <=
 half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
 Global optimum of F on this domain: f* = -4189.828872721624 (A: analytic reduction, 1-D root to machine precision); the optimal solution is listed in
 reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
 Variables: x1:C x2:K x3:D x4:C x5:D x6:K x7:C x8:K x9:D x10:C


Python version of the MATLAB wrapper schwefel_cdk_10D.m: same instance, same values (level indices 1..m as in MATLAB).
Usage:  from schwefel_cdk_10D import schwefel_cdk_10D;  PB = schwefel_cdk_10D();  f = PB.func_idx(PB.x0_idx)
"""
# Copyright (C) 2026 J. F. A. Madeira. LGPL-3.0-or-later.
import math
from types import SimpleNamespace

_NAME = 'schwefel_cdk_10D'
_N = 10
_TYPES = ['C', 'K', 'D', 'C', 'D', 'K', 'C', 'K', 'D', 'C']
_CARD = [math.inf, 17.0, 17.0, math.inf, 17.0, 17.0, math.inf, 17.0, 17.0, math.inf]
_LB = [-343.245849609375, -300.1884875273078, -344.09664, -322.0616409829238, -330.47947544566625, -366.20111340639335, -335.504962823721, -309.3069420891491, -344.6511054491658, -301.55848948296364]
_UB = [500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0]
_LBW = [-343.245849609375, -300.1884875273078, -344.09664, -32206.164098292382, -330.47947544566625, -366.20111340639335, -3.3550496282372095, -309.3069420891491, -344.6511054491658, -301.55848948296364]
_UBW = [500.0, 500.0, 500.0, 50000.0, 500.0, 500.0, 5.0, 500.0, 500.0, 500.0]
_TB = [math.nan, 399.97643905908654, 447.24396, math.nan, 396.19006556929173, 445.86243041210037, math.nan, 398.83663223885634, 447.20930590942714, math.nan]
_TI = [math.nan, 420.9687, 420.9687, math.nan, 420.9687, 420.9687, math.nan, 420.9687, 420.9687, math.nan]
_CW = [420.9687, math.nan, math.nan, 42096.869999999995, math.nan, math.nan, 4.209687000000001, math.nan, math.nan, 420.9687]
_CX = [420.9687, math.nan, math.nan, 420.9687, math.nan, math.nan, 420.9687, math.nan, math.nan, 420.9687]
_M = 17
_PERM = {1: [8, 5, 16, 6, 10, 13, 12, 1, 17, 2, 11, 7, 3, 15, 14, 9, 4], 5: [1, 17, 11, 8, 7, 13, 9, 15, 3, 16, 6, 12, 14, 4, 10, 5, 2], 7: [11, 12, 13, 14, 5, 10, 15, 3, 2, 16, 9, 17, 4, 6, 1, 7, 8]}

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
    for v in x: s += v * math.sin(math.sqrt(abs(v)))
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

def schwefel_cdk_10D():
    PB = SimpleNamespace()
    PB.name = _NAME
    PB.nvar = 10
    PB.var_types = ['C', 'K', 'D', 'C', 'D', 'K', 'C', 'K', 'D', 'C']
    PB.cardinality = [math.inf, 17.0, 17.0, math.inf, 17.0, 17.0, math.inf, 17.0, 17.0, math.inf]
    PB.f_star = -4189.828872721624
    PB.f_star_justification = 'A'
    PB.n_minimizers_documented = 1
    PB.description = 'DFO Benchmark Suite v2.1.0 core (Python): CDK adaptation of schwefel (F = f o T).'
    PB.x0_idx = [-0.2654820680618286, 13.0, 7.0, 6483.117630599907, 3.0, 9.0, 1.3624158040302476, 3.0, 9.0, -239.57401651708273]
    PB.x0_alt_idx = [[217.13383853435516, 12.0, 15.0, 10078.13938416955, 14.0, 7.0, -0.9793109083822573, 7.0, 5.0, -99.46807364669556], [-207.78301537036896, 13.0, 10.0, 13673.161137739186, 11.0, 5.0, 4.515328893945124, 12.0, 16.0, 40.63786922369172], [325.83349883556366, 13.0, 3.0, 5764.11327988598, 6.0, 1.0, 2.1736021815326194, 1.0, 12.0, 180.74381209407892], [-104.0242487192154, 12.0, 13.0, 9359.13503345563, 2.0, 15.0, -0.16812453087988688, 7.0, 8.0, 320.8497549644661], [113.3750718832016, 13.0, 6.0, 12954.156787025258, 13.0, 12.0, -2.509851243292391, 12.0, 4.0, 460.95569783485337], [-311.5417820215225, 13.0, 17.0, 3401.670413254502, 9.0, 10.0, 3.4745614662062327, 17.0, 16.0, -199.54374712554352], [459.2376273870468, 12.0, 9.0, 6996.6921668241375, 4.0, 8.0, 0.6430618466224849, 5.0, 12.0, -59.43780425515632], [26.909433007240295, 12.0, 1.0, 10591.713920393773, 17.0, 5.0, -1.6986648657900205, 10.0, 8.0, 80.66813861523087], [244.3087536096573, 13.0, 12.0, 14186.735673963409, 12.0, 3.0, 4.285747843708604, 16.0, 3.0, 220.77408148561813]]
    PB.suite_release = '2.1.0'
    PB.base_function = 'schwefel'
    PB.m = 17
    PB.pattern = 'CDK'
    PB.strategy = 'moderate'
    PB.scale_factors_C = [1.0, 100.0, 0.01, 1.0]
    PB.lb_orig = [-343.245849609375, -300.1884875273078, -344.09664, -322.0616409829238, -330.47947544566625, -366.20111340639335, -335.504962823721, -309.3069420891491, -344.6511054491658, -301.55848948296364]
    PB.ub_orig = [500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0]
    PB.lb_C_work = [-343.245849609375, -300.1884875273078, -344.09664, -32206.164098292382, -330.47947544566625, -366.20111340639335, -3.3550496282372095, -309.3069420891491, -344.6511054491658, -301.55848948296364]
    PB.ub_C_work = [500.0, 500.0, 500.0, 50000.0, 500.0, 500.0, 5.0, 500.0, 500.0, 500.0]
    PB.ProblemData = _problem_data()
    PB.func_idx = _objective
    PB.func_F = lambda xr: _objective(_cell2idx(xr))
    PB.idx2cell = _idx2cell
    PB.cell2idx = _cell2idx
    PB.decode = lambda xr: _decode(_cell2idx(xr))
    return PB
