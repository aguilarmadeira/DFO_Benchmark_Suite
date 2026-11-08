"""transistor_cd_9D  DFO Benchmark Suite v2.1.0 core, Python version (self-contained; standard library only).
TRANSISTOR_CD_9D  DFO Benchmark Suite v2.1.0 core: transistor, dimension 9, pattern CD, strategy moderate.

 Copyright (C) 2026 J. F. A. Madeira

 Objective F(z) = f(T(z)), where f is the original function (embedded below, translated from the MATLAB original) and T acts on the finite
 coordinates only: uniform grid of 25 levels on the interval of the domain, piecewise-affine map with |T - Id| <=
 half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
 Global optimum of F on this domain: f* = 1.60049326654736e-24 (E: exact by structure); the optimal solution is listed in
 reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
 Variables: x1:C x2:D x3:C x4:D x5:C x6:D x7:C x8:D x9:C


Python version of the MATLAB wrapper transistor_cd_9D.m: same instance, same values (level indices 1..m as in MATLAB).
Usage:  from transistor_cd_9D import transistor_cd_9D;  PB = transistor_cd_9D();  f = PB.func_idx(PB.x0_idx)
"""
# Copyright (C) 2026 J. F. A. Madeira. LGPL-3.0-or-later.
import math
from types import SimpleNamespace

_NAME = 'transistor_cd_9D'
_N = 9
_TYPES = ['C', 'D', 'C', 'D', 'C', 'D', 'C', 'D', 'C']
_CARD = [math.inf, 25.0, math.inf, 25.0, math.inf, 25.0, math.inf, 25.0, math.inf]
_LB = [-10.0, -10.0, -10.0, -10.0, -7.350085376681921, -6.992661321382304, -10.0, -10.0, -10.0]
_UB = [6.8211669921875, 7.03312503175329, 6.7027328, 6.180008329862558, 10.0, 10.0, 6.7876079069934505, 6.403313357018439, 7.316462562669516]
_LBW = [-10.0, -10.0, -1000.0, -10.0, -0.07350085376681921, -6.992661321382304, -10.0, -10.0, -1000.0]
_UBW = [6.8211669921875, 7.03312503175329, 670.2732799999999, 6.180008329862558, 0.09999999999999999, 10.0, 6.7876079069934505, 6.403313357018439, 731.6462562669515]
_TB = [math.nan, 0.6457031448458057, math.nan, 2.1350062473969214, math.nan, 7.8759173348272125, math.nan, 0.9355422380122942, math.nan]
_TI = [math.nan, 0.4499874719815851, math.nan, 2.000068541623683, math.nan, 7.999692684218543, math.nan, 0.9999877234567898, math.nan]
_CW = [0.8999999526168558, math.nan, 100.00064824653873, math.nan, 0.07999971440507851, math.nan, 5.000031275929672, math.nan, 200.00524834861903]
_CX = [0.8999999526168565, math.nan, 1.000006482465387, math.nan, 7.999971440507853, math.nan, 5.000031275929672, math.nan, 2.00005248348619]
_M = 25
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

def transistor_cd_9D():
    PB = SimpleNamespace()
    PB.name = _NAME
    PB.nvar = 9
    PB.var_types = ['C', 'D', 'C', 'D', 'C', 'D', 'C', 'D', 'C']
    PB.cardinality = [math.inf, 25.0, math.inf, 25.0, math.inf, 25.0, math.inf, 25.0, math.inf]
    PB.f_star = 1.6004932665473596e-24
    PB.f_star_justification = 'E'
    PB.n_minimizers_documented = 1
    PB.description = 'DFO Benchmark Suite v2.1.0 core (Python): CD adaptation of transistor (F = f o T).'
    PB.x0_idx = [-3.4210181832313538, 2.0, -162.99265392640007, 20.0, 0.03364993646332759, 18.0, -1.4507235938738816, 24.0, -698.5594829642966]
    PB.x0_alt_idx = [[0.9156889319419861, 4.0, -988.7757635584, 21.0, -0.02070758120817024, 15.0, -6.155900896526026, 6.0, 587.899834073344], [-7.560602247714996, 2.0, 173.7344393215999, 19.0, -0.05942264775117949, 10.0, 4.884351728651227, 14.0, 145.98632822071932], [3.084042489528656, 2.0, -664.8763691008, 20.0, 0.05828679699069712, 7.0, 0.17917442599908284, 22.0, -295.92717763190535], [-5.490810215473175, 3.0, 497.6338337791999, 21.0, 0.01957173044768787, 4.0, -4.526002876653062, 6.0, -737.8406834845299], [-1.1541031002998352, 1.0, -328.1492758528001, 18.0, -0.03620782823549074, 23.0, -9.231180179305206, 13.0, 548.6186335531106], [-9.630394279956818, 3.0, -769.1014217728, 19.0, 0.09714406763487442, 20.0, 2.7931618163613834, 21.0, 106.70512770048595], [6.533695876598358, 3.0, 8.577817395199986, 20.0, 0.042786549963376594, 17.0, -2.8961048567800978, 4.0, -335.2083781521385], [-2.090437591075897, 1.0, -432.37432852480003, 21.0, 0.004071483420367353, 13.0, -7.601282159432241, 12.0, -777.1218840047632], [2.2462695240974426, 2.0, 345.30491064319995, 19.0, -0.05028603425113047, 10.0, 4.423059836234348, 19.0, 434.1907885593878]]
    PB.suite_release = '2.1.0'
    PB.base_function = 'transistor'
    PB.m = 25
    PB.pattern = 'CD'
    PB.strategy = 'moderate'
    PB.scale_factors_C = [1.0, 100.0, 0.01, 1.0, 100.0]
    PB.lb_orig = [-10.0, -10.0, -10.0, -10.0, -7.350085376681921, -6.992661321382304, -10.0, -10.0, -10.0]
    PB.ub_orig = [6.8211669921875, 7.03312503175329, 6.7027328, 6.180008329862558, 10.0, 10.0, 6.7876079069934505, 6.403313357018439, 7.316462562669516]
    PB.lb_C_work = [-10.0, -10.0, -1000.0, -10.0, -0.07350085376681921, -6.992661321382304, -10.0, -10.0, -1000.0]
    PB.ub_C_work = [6.8211669921875, 7.03312503175329, 670.2732799999999, 6.180008329862558, 0.09999999999999999, 10.0, 6.7876079069934505, 6.403313357018439, 731.6462562669515]
    PB.ProblemData = _problem_data()
    PB.func_idx = _objective
    PB.func_F = lambda xr: _objective(_cell2idx(xr))
    PB.idx2cell = _idx2cell
    PB.cell2idx = _cell2idx
    PB.decode = lambda xr: _decode(_cell2idx(xr))
    return PB
