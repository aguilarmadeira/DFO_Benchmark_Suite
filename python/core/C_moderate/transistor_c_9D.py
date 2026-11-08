"""transistor_c_9D  DFO Benchmark Suite v2.1.0 core, Python version (self-contained; standard library only).
TRANSISTOR_C_9D  DFO Benchmark Suite v2.1.0 core: transistor, dimension 9, pattern C, strategy moderate.

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
_LBW = [-10.0, -1000.0, -0.1, -10.0, -735.008537668192, -0.06992661321382304, -10.0, -1000.0, -0.1]
_UBW = [6.8211669921875, 703.3125031753291, 0.067027328, 6.1800083298625585, 1000.0, 0.1, 6.7876079069934505, 640.331335701844, 0.07316462562669515]
_TB = [math.nan, math.nan, math.nan, math.nan, math.nan, math.nan, math.nan, math.nan, math.nan]
_TI = [math.nan, math.nan, math.nan, math.nan, math.nan, math.nan, math.nan, math.nan, math.nan]
_CW = [0.8999999526168558, 44.998747198158526, 0.010000064824653876, 2.0000685416236834, 799.9971440507852, 0.07999692684218543, 5.000031275929672, 99.99877234567907, 0.020000524834861907]
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
    PB.x0_idx = [-3.4210181832313538, -926.7894946051756, -0.016299265392639997, 2.8981824003985572, 336.49936463327595, 0.051195406036448635, -1.4507235938738816, 545.388699709187, -0.06985594829642966]
    PB.x0_alt_idx = [[0.9156889319419861, -821.6467474955874, -0.09887757635584, 3.605762939605375, -207.07581208170234, 0.030080264749523897, -6.155900896526026, -644.8619283398107, 0.058789983407334406], [-7.560602247714996, -940.8085275531207, 0.017373443932160004, 2.0490857533503757, -594.2264775117948, -0.004028809637046868, 4.884351728651227, -113.2309968242547, 0.014598632822071939], [3.084042489528656, -910.4339561659064, -0.06648763691008, 2.756666292557192, 582.8679699069712, -0.025143950923971634, 0.17917442599908284, 418.3999346913015, -0.02959271776319053], [-5.490810215473175, -861.3673408480986, 0.049763383377919995, 3.4642468317640116, 195.71730447687878, -0.046259092210896385, -4.526002876653062, -685.7566153794689, -0.073784068348453], [-1.1541031002998352, -980.5291209056318, -0.03281492758528001, 1.5841042561573229, -362.0782823549074, 0.08855296369793092, -9.231180179305206, -154.12568386391274, 0.05486186335531107], [-9.630394279956818, -875.3863737960436, -0.07691014217728001, 2.291684795364141, 971.4406763487443, 0.06743782241100615, 2.7931618163613834, 377.50524765164323, 0.010670512770048601], [6.533695876598358, -851.2424837190271, 0.0008577817395200027, 2.999265334570957, 427.86549963376615, 0.04632268112408139, -2.8961048567800978, -812.7453803973548, -0.03352083781521385], [-2.090437591075897, -970.4042637765604, -0.04323743285248001, 3.706845873777775, 40.714834203673604, 0.012213606737510624, -7.601282159432241, -281.1144488817987, -0.07771218840047632], [2.2462695240974426, -921.3376484587525, 0.03453049106431999, 2.150168687522779, -502.8603425113046, -0.008901534549414121, 4.423059836234348, 250.51648263375728, 0.043419078855938764]]
    PB.suite_release = '2.1.0'
    PB.base_function = 'transistor'
    PB.m = 0
    PB.pattern = 'C'
    PB.strategy = 'moderate'
    PB.scale_factors_C = [1.0, 100.0, 0.01, 1.0, 100.0, 0.01, 1.0, 100.0, 0.01]
    PB.lb_orig = [-10.0, -10.0, -10.0, -10.0, -7.350085376681921, -6.992661321382304, -10.0, -10.0, -10.0]
    PB.ub_orig = [6.8211669921875, 7.03312503175329, 6.7027328, 6.180008329862558, 10.0, 10.0, 6.7876079069934505, 6.403313357018439, 7.316462562669516]
    PB.lb_C_work = [-10.0, -1000.0, -0.1, -10.0, -735.008537668192, -0.06992661321382304, -10.0, -1000.0, -0.1]
    PB.ub_C_work = [6.8211669921875, 703.3125031753291, 0.067027328, 6.1800083298625585, 1000.0, 0.1, 6.7876079069934505, 640.331335701844, 0.07316462562669515]
    PB.ProblemData = _problem_data()
    PB.func_idx = _objective
    PB.func_F = lambda xr: _objective(_cell2idx(xr))
    PB.idx2cell = _idx2cell
    PB.cell2idx = _cell2idx
    PB.decode = lambda xr: _decode(_cell2idx(xr))
    return PB
