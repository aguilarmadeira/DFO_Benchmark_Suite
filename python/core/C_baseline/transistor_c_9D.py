"""transistor_c_9D  DFO Benchmark Suite v2.1.0 core, Python version (self-contained; standard library only).
TRANSISTOR_C_9D  DFO Benchmark Suite v2.1.0 core: transistor, dimension 9, pattern C, strategy baseline.

 Copyright (C) 2026 J. F. A. Madeira

 Objective F(z) = f(T(z)), where f is the original function (embedded below, translated from the MATLAB original) and T acts on the finite
 coordinates only: uniform grid of 0 levels on the interval of the domain, piecewise-affine map with |T - Id| <=
 half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
 Global optimum of F on this domain: f* = 1.60049326654736e-24 (E: exact by structure); the optimal solution is listed in
 reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
 Variables: x1:C x2:C x3:C x4:C x5:C x6:C x7:C x8:C x9:C


Python version of the MATLAB wrapper transistor_c_9D.m: same instance, same values (level indices 1..m as in MATLAB).
Usage:  from transistor_c_9D import transistor_c_9D;  PB = transistor_c_9D();  f = PB.func_idx(PB.x0_idx)
"""
# Copyright (C) 2026 J. F. A. Madeira. LGPL-3.0-or-later.
import math
from types import SimpleNamespace

_NAME = 'transistor_c_9D'
_N = 9
_TYPES = ['C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C']
_CARD = [math.inf, math.inf, math.inf, math.inf, math.inf, math.inf, math.inf, math.inf, math.inf]
_LB = [-10.0, -10.0, -10.0, -10.0, -7.350085376681921, -6.992661321382304, -10.0, -10.0, -10.0]
_UB = [6.8211669921875, 7.03312503175329, 6.7027328, 6.180008329862558, 10.0, 10.0, 6.7876079069934505, 6.403313357018439, 7.316462562669516]
_LBW = [-10.0, -10.0, -10.0, -10.0, -7.350085376681921, -6.992661321382304, -10.0, -10.0, -10.0]
_UBW = [6.8211669921875, 7.03312503175329, 6.7027328, 6.180008329862558, 10.0, 10.0, 6.7876079069934505, 6.403313357018439, 7.316462562669516]
_TB = [math.nan, math.nan, math.nan, math.nan, math.nan, math.nan, math.nan, math.nan, math.nan]
_TI = [math.nan, math.nan, math.nan, math.nan, math.nan, math.nan, math.nan, math.nan, math.nan]
_CW = [0.8999999526168558, 0.4499874719815846, 1.0000064824653876, 2.0000685416236834, 7.999971440507853, 7.999692684218543, 5.000031275929672, 0.9999877234567904, 2.00005248348619]
_CX = [0.8999999526168565, 0.4499874719815851, 1.000006482465387, 2.000068541623683, 7.999971440507853, 7.999692684218543, 5.000031275929672, 0.9999877234567898, 2.00005248348619]
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

def transistor_c_9D():
    PB = SimpleNamespace()
    PB.name = _NAME
    PB.nvar = 9
    PB.var_types = ['C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C']
    PB.cardinality = [math.inf, math.inf, math.inf, math.inf, math.inf, math.inf, math.inf, math.inf, math.inf]
    PB.f_star = 1.6004932665473596e-24
    PB.f_star_justification = 'E'
    PB.n_minimizers_documented = 1
    PB.description = 'DFO Benchmark Suite v2.1.0 core (Python): C adaptation of transistor (F = f o T).'
    PB.x0_idx = [-3.4210181832313538, -9.267894946051756, -1.6299265392639999, 2.8981824003985572, 3.36499364633276, 5.119540603644863, -1.4507235938738816, 5.453886997091871, -6.985594829642967]
    PB.x0_alt_idx = [[0.9156889319419861, -8.216467474955873, -9.887757635584, 3.605762939605375, -2.0707581208170245, 3.0080264749523895, -6.155900896526026, -6.4486192833981075, 5.878998340733439], [-7.560602247714996, -9.408085275531207, 1.7373443932160004, 2.0490857533503757, -5.942264775117948, -0.4028809637046873, 4.884351728651227, -1.1323099682425468, 1.459863282207193], [3.084042489528656, -9.104339561659064, -6.648763691008, 2.756666292557192, 5.828679699069712, -2.514395092397163, 0.17917442599908284, 4.183999346913016, -2.9592717763190537], [-5.490810215473175, -8.613673408480985, 4.976338337791999, 3.4642468317640116, 1.957173044768787, -4.625909221089639, -4.526002876653062, -6.857566153794689, -7.378406834845299], [-1.1541031002998352, -9.805291209056318, -3.2814927585280005, 1.5841042561573229, -3.6207828235490744, 8.855296369793091, -9.231180179305206, -1.5412568386391268, 5.486186335531107], [-9.630394279956818, -8.753863737960437, -7.691014217728, 2.291684795364141, 9.714406763487442, 6.7437822411006145, 2.7931618163613834, 3.775052476516434, 1.0670512770048592], [6.533695876598358, -8.51242483719027, 0.08577817395200071, 2.999265334570957, 4.278654996337661, 4.6322681124081395, -2.8961048567800978, -8.127453803973548, -3.3520837815213858], [-2.090437591075897, -9.704042637765603, -4.323743285248001, 3.706845873777775, 0.4071483420367352, 1.2213606737510627, -7.601282159432241, -2.8111444888179866, -7.771218840047632], [2.2462695240974426, -9.213376484587526, 3.453049106432001, 2.150168687522779, -5.028603425113047, -0.8901534549414123, 4.423059836234348, 2.505164826337573, 4.341907885593876]]
    PB.suite_release = '2.1.0'
    PB.base_function = 'transistor'
    PB.m = 0
    PB.pattern = 'C'
    PB.strategy = 'baseline'
    PB.scale_factors_C = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
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
