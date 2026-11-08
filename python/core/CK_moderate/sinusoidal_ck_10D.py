"""sinusoidal_ck_10D  DFO Benchmark Suite v2.1.0 core, Python version (self-contained; standard library only).
SINUSOIDAL_CK_10D  DFO Benchmark Suite v2.1.0 core: sinusoidal, dimension 10, pattern CK, strategy moderate.

 Copyright (C) 2026 J. F. A. Madeira

 Objective F(z) = f(T(z)), where f is the original function (embedded below, translated from the MATLAB original) and T acts on the finite
 coordinates only: uniform grid of 17 levels on the interval of the domain, piecewise-affine map with |T - Id| <=
 half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
 Global optimum of F on this domain: f* = -3.5 (E: exact by structure); the optimal solution is listed in
 reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
 Variables: x1:C x2:K x3:C x4:K x5:C x6:K x7:C x8:K x9:C x10:K


Python version of the MATLAB wrapper sinusoidal_ck_10D.m: same instance, same values (level indices 1..m as in MATLAB).
Usage:  from sinusoidal_ck_10D import sinusoidal_ck_10D;  PB = sinusoidal_ck_10D();  f = PB.func_idx(PB.x0_idx)
"""
# Copyright (C) 2026 J. F. A. Madeira. LGPL-3.0-or-later.
import math
from types import SimpleNamespace

_NAME = 'sinusoidal_ck_10D'
_N = 10
_TYPES = ['C', 'K', 'C', 'K', 'C', 'K', 'C', 'K', 'C', 'K']
_CARD = [math.inf, 17.0, math.inf, 17.0, math.inf, 17.0, math.inf, 17.0, math.inf, 17.0]
_LB = [32.1532470703125, 26.235208047553726, 0.0, 0.0, 0.0, 0.0, 34.79734198584787, 26.401758734202474, 0.0, 0.0]
_UB = [180.0, 180.0, 148.8125952, 147.71395251978342, 153.96399153063317, 149.13158502853543, 180.0, 180.0, 148.20241637215418, 153.82416663249825]
_LBW = [32.15324707031249, 26.235208047553726, 0.0, 0.0, 0.0, 0.0, 34.79734198584788, 26.401758734202474, 0.0, 0.0]
_UBW = [180.0, 180.0, 14881.25952, 147.71395251978342, 1.5396399153063316, 149.13158502853543, 180.0, 180.0, 14820.241637215418, 153.82416663249825]
_TB = [math.nan, 122.33820301783265, math.nan, 120.01758642232403, math.nan, 121.16941283568504, math.nan, 122.40065952532592, math.nan, 115.36812497437369]
_TI = [math.nan, 120.0, math.nan, 120.0, math.nan, 120.0, math.nan, 120.0, math.nan, 120.0]
_CW = [119.99999999999999, math.nan, 11999.999999999998, math.nan, 1.2, math.nan, 120.0, math.nan, 12000.0, math.nan]
_CX = [120.0, math.nan, 120.0, math.nan, 120.0, math.nan, 120.0, math.nan, 120.0, math.nan]
_M = 17
_PERM = {1: [3, 9, 1, 2, 7, 5, 6, 4, 8, 15, 16, 10, 17, 13, 11, 12, 14], 3: [12, 17, 16, 7, 3, 9, 4, 11, 8, 1, 14, 6, 13, 15, 2, 5, 10], 5: [10, 4, 1, 15, 11, 16, 2, 9, 8, 12, 5, 14, 13, 3, 7, 6, 17], 7: [5, 1, 2, 10, 17, 3, 9, 16, 7, 4, 8, 11, 15, 6, 13, 14, 12], 9: [12, 3, 9, 16, 10, 11, 5, 15, 17, 6, 1, 8, 13, 7, 4, 2, 14]}

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

def _sind(v):
    # sine of an angle in degrees, exact at multiples of 90 (as MATLAB sind)
    r = math.fmod(v, 360.0)
    n = math.floor(r / 90 + 0.5) if r >= 0 else -math.floor(-r / 90 + 0.5)
    t = (r - 90 * n) * math.pi / 180
    q = int(n) % 4
    return (math.sin(t), math.cos(t), -math.sin(t), -math.cos(t))[q]
def _orig(x):
    A = 2.5; B = 5; z = 30
    p1 = 1.0
    for v in x: p1 *= _sind(v - z)
    p2 = 1.0
    for v in x: p2 *= _sind(B * (v - z))
    return -(A * p1 + p2)

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

def sinusoidal_ck_10D():
    PB = SimpleNamespace()
    PB.name = _NAME
    PB.nvar = 10
    PB.var_types = ['C', 'K', 'C', 'K', 'C', 'K', 'C', 'K', 'C', 'K']
    PB.cardinality = [math.inf, 17.0, math.inf, 17.0, math.inf, 17.0, math.inf, 17.0, math.inf, 17.0]
    PB.f_star = -3.5
    PB.f_star_justification = 'E'
    PB.n_minimizers_documented = 1
    PB.description = 'DFO Benchmark Suite v2.1.0 core (Python): CK adaptation of sinusoidal (F = f o T).'
    PB.x0_idx = [138.49029153585434, 1.0, 9243.0479130624, 9.0, 1.2180622320192092, 1.0, 176.5716449563115, 8.0, 12216.010798030198, 6.0]
    PB.x0_alt_idx = [[65.72196782827376, 2.0, 5314.3953997824, 10.0, 0.7356957070885252, 15.0, 135.87470620494014, 13.0, 7790.767281304333, 9.0], [157.54866203069685, 3.0, 12243.1098322944, 8.0, 0.392139692929261, 13.0, 86.66598941406627, 2.0, 4008.6640279506814, 12.0], [84.7803383231163, 1.0, 8314.4573190144, 9.0, 1.4366887864841953, 10.0, 45.96905066269492, 7.0, 226.5607745970303, 14.0], [120.2982106089592, 2.0, 842.8745392128001, 9.0, 1.0931327723249311, 8.0, 149.97233858286614, 12.0, 11236.683578804324, 1.0], [47.52988690137862, 3.0, 11200.231165132802, 10.0, 0.5981471100513909, 5.0, 109.27539983149481, 2.0, 7454.5803254506745, 4.0], [167.07784727811813, 2.0, 3842.9364584447994, 8.0, 0.2545910958921268, 3.0, 68.57846108012345, 7.0, 3672.477072097024, 7.0], [94.30952357053755, 2.0, 14200.293084364801, 9.0, 1.2991401894470611, 17.0, 164.06997096079218, 11.0, 14682.599876304319, 10.0], [129.39425107240675, 2.0, 6842.9983776767995, 10.0, 0.9555841752877969, 14.0, 123.37303220942084, 16.0, 10900.496622950668, 12.0], [56.6259273648262, 1.0, 2914.3458643968, 8.0, 0.47321765035711294, 12.0, 82.67609345804948, 6.0, 6475.253106224802, 15.0]]
    PB.suite_release = '2.1.0'
    PB.base_function = 'sinusoidal'
    PB.m = 17
    PB.pattern = 'CK'
    PB.strategy = 'moderate'
    PB.scale_factors_C = [1.0, 100.0, 0.01, 1.0, 100.0]
    PB.lb_orig = [32.1532470703125, 26.235208047553726, 0.0, 0.0, 0.0, 0.0, 34.79734198584787, 26.401758734202474, 0.0, 0.0]
    PB.ub_orig = [180.0, 180.0, 148.8125952, 147.71395251978342, 153.96399153063317, 149.13158502853543, 180.0, 180.0, 148.20241637215418, 153.82416663249825]
    PB.lb_C_work = [32.15324707031249, 26.235208047553726, 0.0, 0.0, 0.0, 0.0, 34.79734198584788, 26.401758734202474, 0.0, 0.0]
    PB.ub_C_work = [180.0, 180.0, 14881.25952, 147.71395251978342, 1.5396399153063316, 149.13158502853543, 180.0, 180.0, 14820.241637215418, 153.82416663249825]
    PB.ProblemData = _problem_data()
    PB.func_idx = _objective
    PB.func_F = lambda xr: _objective(_cell2idx(xr))
    PB.idx2cell = _idx2cell
    PB.cell2idx = _cell2idx
    PB.decode = lambda xr: _decode(_cell2idx(xr))
    return PB
