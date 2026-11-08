"""sinusoidal_c_10D  DFO Benchmark Suite v2.1.0 core, Python version (self-contained; standard library only).
SINUSOIDAL_C_10D  DFO Benchmark Suite v2.1.0 core: sinusoidal, dimension 10, pattern C, strategy moderate.

 Copyright (C) 2026 J. F. A. Madeira

 Objective F(z) = f(T(z)), where f is the original function (embedded below, translated from the MATLAB original) and T acts on the finite
 coordinates only: uniform grid of 0 levels on the interval of the domain, piecewise-affine map with |T - Id| <=
 half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
 Global optimum of F on this domain: f* = -3.5 (E: exact by structure); the optimal solution is listed in
 reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
 Variables: x1:C x2:C x3:C x4:C x5:C x6:C x7:C x8:C x9:C x10:C


Python version of the MATLAB wrapper sinusoidal_c_10D.m: same instance, same values (level indices 1..m as in MATLAB).
Usage:  from sinusoidal_c_10D import sinusoidal_c_10D;  PB = sinusoidal_c_10D();  f = PB.func_idx(PB.x0_idx)
"""
# Copyright (C) 2026 J. F. A. Madeira. LGPL-3.0-or-later.
import math
from types import SimpleNamespace

_NAME = 'sinusoidal_c_10D'
_N = 10
_TYPES = ['C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C']
_CARD = [math.inf, math.inf, math.inf, math.inf, math.inf, math.inf, math.inf, math.inf, math.inf, math.inf]
_LB = [32.1532470703125, 26.235208047553726, 0.0, 0.0, 0.0, 0.0, 34.79734198584787, 26.401758734202474, 0.0, 0.0]
_UB = [180.0, 180.0, 148.8125952, 147.71395251978342, 153.96399153063317, 149.13158502853543, 180.0, 180.0, 148.20241637215418, 153.82416663249825]
_LBW = [32.15324707031249, 2623.5208047553724, 0.0, 0.0, 0.0, 0.0, 34.79734198584788, 2640.1758734202476, 0.0, 0.0]
_UBW = [180.0, 18000.0, 1.488125952, 147.71395251978342, 15396.399153063317, 1.4913158502853543, 180.0, 18000.0, 1.4820241637215419, 153.82416663249825]
_TB = [math.nan, math.nan, math.nan, math.nan, math.nan, math.nan, math.nan, math.nan, math.nan, math.nan]
_TI = [math.nan, math.nan, math.nan, math.nan, math.nan, math.nan, math.nan, math.nan, math.nan, math.nan]
_CW = [119.99999999999999, 12000.0, 1.2, 120.0, 12000.000000000002, 1.2, 120.0, 12000.0, 1.2000000000000002, 120.0]
_CX = [120.0, 120.0, 120.0, 120.0, 120.0, 120.0, 120.0, 120.0, 120.0, 120.0]
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

def sinusoidal_c_10D():
    PB = SimpleNamespace()
    PB.name = _NAME
    PB.nvar = 10
    PB.var_types = ['C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C']
    PB.cardinality = [math.inf, math.inf, math.inf, math.inf, math.inf, math.inf, math.inf, math.inf, math.inf, math.inf]
    PB.f_star = -3.5
    PB.f_star_justification = 'E'
    PB.n_minimizers_documented = 1
    PB.description = 'DFO Benchmark Suite v2.1.0 core (Python): C adaptation of sinusoidal (F = f o T).'
    PB.x0_idx = [138.49029153585434, 2714.9219156553177, 0.9243047913062402, 72.53425656844009, 12180.622320192093, 0.030545841266654963, 176.5716449563115, 9835.26479581429, 1.22160107980302, 49.02518542106724]
    PB.x0_alt_idx = [[65.72196782827376, 3664.0872980778254, 0.53143953997824, 78.99405041041312, 7356.9570708852525, 1.3277259003906023, 135.87470620494014, 14007.207112160417, 0.7790767281304333, 75.91240598584398], [157.54866203069685, 4107.031143208329, 1.22431098322944, 64.78250395807244, 3921.3969292926104, 1.1424144633728954, 86.66598941406627, 3668.0457194765395, 0.40086640279506813, 102.79962655062073], [84.7803383231163, 2862.569864032152, 0.83144573190144, 71.24229780004548, 14366.887864841954, 0.843065218959677, 45.96905066269492, 8646.16046410211, 0.022656077459703033, 129.68684711539748], [120.2982106089592, 3305.5137091626557, 0.08428745392128001, 77.70209164201853, 10931.327723249313, 0.6577537819419701, 149.97233858286614, 13624.275208727679, 1.1236683578804325, 2.93280731002139], [47.52988690137862, 4254.679091585163, 1.1200231165132801, 84.16188548399155, 5981.47110051391, 0.35840453752875157, 109.27539983149481, 3285.113816043803, 0.7454580325450675, 29.820027874798136], [167.07784727811813, 3178.958324839655, 0.38429364584448, 66.99729041817749, 2545.9109589212685, 0.17309310051104476, 68.57846108012345, 8263.228560669373, 0.3672477072097024, 56.707248439574876], [94.30952357053755, 3396.914820062601, 1.4200293084364801, 73.45708426015051, 12991.401894470613, 1.470273159634992, 164.06997096079218, 12435.170877015498, 1.4682599876304319, 83.59446900435164], [129.39425107240675, 3839.8586651931046, 0.68429983776768, 79.91687810212353, 9555.84175287797, 1.1709239152217734, 123.37303220942084, 17413.28562164107, 1.0900496622950668, 110.48168956912838], [56.6259273648262, 2764.137898447596, 0.29143458643968, 65.70533164978288, 4732.17650357113, 0.9856124782040667, 82.67609345804948, 7074.124228957192, 0.6475253106224802, 137.3689101339051]]
    PB.suite_release = '2.1.0'
    PB.base_function = 'sinusoidal'
    PB.m = 0
    PB.pattern = 'C'
    PB.strategy = 'moderate'
    PB.scale_factors_C = [1.0, 100.0, 0.01, 1.0, 100.0, 0.01, 1.0, 100.0, 0.01, 1.0]
    PB.lb_orig = [32.1532470703125, 26.235208047553726, 0.0, 0.0, 0.0, 0.0, 34.79734198584787, 26.401758734202474, 0.0, 0.0]
    PB.ub_orig = [180.0, 180.0, 148.8125952, 147.71395251978342, 153.96399153063317, 149.13158502853543, 180.0, 180.0, 148.20241637215418, 153.82416663249825]
    PB.lb_C_work = [32.15324707031249, 2623.5208047553724, 0.0, 0.0, 0.0, 0.0, 34.79734198584788, 2640.1758734202476, 0.0, 0.0]
    PB.ub_C_work = [180.0, 18000.0, 1.488125952, 147.71395251978342, 15396.399153063317, 1.4913158502853543, 180.0, 18000.0, 1.4820241637215419, 153.82416663249825]
    PB.ProblemData = _problem_data()
    PB.func_idx = _objective
    PB.func_F = lambda xr: _objective(_cell2idx(xr))
    PB.idx2cell = _idx2cell
    PB.cell2idx = _cell2idx
    PB.decode = lambda xr: _decode(_cell2idx(xr))
    return PB
