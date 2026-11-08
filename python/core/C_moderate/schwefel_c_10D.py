"""schwefel_c_10D  DFO Benchmark Suite v2.1.0 core, Python version (self-contained; standard library only).
SCHWEFEL_C_10D  DFO Benchmark Suite v2.1.0 core: schwefel, dimension 10, pattern C, strategy moderate.

 Copyright (C) 2026 J. F. A. Madeira

 Objective F(z) = f(T(z)), where f is the original function (embedded below, translated from the MATLAB original) and T acts on the finite
 coordinates only: uniform grid of 0 levels on the interval of the domain, piecewise-affine map with |T - Id| <=
 half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
 Global optimum of F on this domain: f* = -4189.828872721624 (A: analytic reduction, 1-D root to machine precision); the optimal solution is listed in
 reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
 Variables: x1:C x2:C x3:C x4:C x5:C x6:C x7:C x8:C x9:C x10:C


Python version of the MATLAB wrapper schwefel_c_10D.m: same instance, same values (level indices 1..m as in MATLAB).
Usage:  from schwefel_c_10D import schwefel_c_10D;  PB = schwefel_c_10D();  f = PB.func_idx(PB.x0_idx)
"""
# Copyright (C) 2026 J. F. A. Madeira. LGPL-3.0-or-later.
import math
from types import SimpleNamespace

_NAME = 'schwefel_c_10D'
_N = 10
_TYPES = ['C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C']
_CARD = [math.inf, math.inf, math.inf, math.inf, math.inf, math.inf, math.inf, math.inf, math.inf, math.inf]
_LB = [-343.245849609375, -300.1884875273078, -344.09664, -322.0616409829238, -330.47947544566625, -366.20111340639335, -335.504962823721, -309.3069420891491, -344.6511054491658, -301.55848948296364]
_UB = [500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0]
_LBW = [-343.245849609375, -30018.848752730784, -3.4409664, -322.0616409829238, -33047.94754456663, -3.662011134063933, -335.5049628237209, -30930.69420891491, -3.4465110544916584, -301.55848948296364]
_UBW = [500.0, 50000.0, 5.0, 500.0, 50000.0, 5.0, 500.0, 50000.0, 5.0, 500.0]
_TB = [math.nan, math.nan, math.nan, math.nan, math.nan, math.nan, math.nan, math.nan, math.nan, math.nan]
_TI = [math.nan, math.nan, math.nan, math.nan, math.nan, math.nan, math.nan, math.nan, math.nan, math.nan]
_CW = [420.9687, 42096.869999999995, 4.209687000000001, 420.9687, 42096.87000000001, 4.209687000000001, 420.9687000000001, 42096.87000000001, 4.209687, 420.9687]
_CX = [420.9687, 420.9687, 420.9687, 420.9687, 420.9687, 420.9687, 420.9687, 420.9687, 420.9687, 420.9687]
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

def schwefel_c_10D():
    PB = SimpleNamespace()
    PB.name = _NAME
    PB.nvar = 10
    PB.var_types = ['C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C']
    PB.cardinality = [math.inf, math.inf, math.inf, math.inf, math.inf, math.inf, math.inf, math.inf, math.inf, math.inf]
    PB.f_star = -4189.828872721624
    PB.f_star_justification = 'A'
    PB.n_minimizers_documented = 1
    PB.description = 'DFO Benchmark Suite v2.1.0 core (Python): C adaptation of schwefel (F = f o T).'
    PB.x0_idx = [-0.2654820680618286, 31376.50022170097, -0.2914730168320001, 64.83117630599907, -20194.546196572097, 0.46200509089322406, 136.24158040302484, -22836.444868294464, 0.9450915534408004, -239.57401651708273]
    PB.x0_alt_idx = [[217.13383853435516, 25778.473765748615, 3.962774048768001, 100.78139384169555, 36148.27621720977, -0.6143394878957853, -97.93109083822569, -854.5403178923007, -1.5769906904129178, -99.46807364669556], [-207.78301537036896, 28083.543482905476, 1.4102258094080002, 136.73161137739186, 17616.91602131474, -1.6906840666847949, 451.5328893945125, 25375.075256742173, 4.698016576912643, 40.63786922369172], [325.83349883556366, 29510.491403050193, -2.438854868992, 57.641132798859815, -8401.862435547988, -3.4293945401131944, 217.36021815326194, -29101.818629037112, 2.5424797293580603, 180.74381209407892], [-104.0242487192154, 23912.46494709783, 3.0470980136960013, 93.5913503345563, -26933.222631443015, 4.105017511409873, -16.812453087988615, -2872.203054402642, 0.386942881803479, 320.8497549644661], [113.3750718832016, 28851.90005529109, -0.8019826647039996, 129.5415678702526, 28728.92539497815, 2.3663070379814712, -250.98512432923908, 23357.412520231832, -1.7685939657511027, 460.95569783485337], [-311.5417820215225, 31156.969772447945, 4.748796839936, 34.016704132545044, 10197.565199083117, 1.2899624591924628, 347.4561466206234, 49587.0280948663, 4.506413301574458, -199.54374712554352], [459.2376273870468, 24388.114253812735, 0.8997161615360008, 69.96692166824141, -15821.213257779604, 0.21361788040345298, 64.30618466224848, -9137.57681514529, 2.3508764540198754, -59.43780425515632], [26.909433007240295, 26693.183970969596, -3.273497626624, 105.91713920393778, 48009.02741696997, -1.5250925930249468, -169.86648657900201, 17092.038759489187, 0.1953396064652937, 80.66813861523087], [244.3087536096573, 31632.61907916285, 2.601414987776, 141.8673567396341, 21990.24896010723, -2.6014371718139566, 428.57478437086047, 43321.654334123654, -2.326742637388425, 220.77408148561813]]
    PB.suite_release = '2.1.0'
    PB.base_function = 'schwefel'
    PB.m = 0
    PB.pattern = 'C'
    PB.strategy = 'moderate'
    PB.scale_factors_C = [1.0, 100.0, 0.01, 1.0, 100.0, 0.01, 1.0, 100.0, 0.01, 1.0]
    PB.lb_orig = [-343.245849609375, -300.1884875273078, -344.09664, -322.0616409829238, -330.47947544566625, -366.20111340639335, -335.504962823721, -309.3069420891491, -344.6511054491658, -301.55848948296364]
    PB.ub_orig = [500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0]
    PB.lb_C_work = [-343.245849609375, -30018.848752730784, -3.4409664, -322.0616409829238, -33047.94754456663, -3.662011134063933, -335.5049628237209, -30930.69420891491, -3.4465110544916584, -301.55848948296364]
    PB.ub_C_work = [500.0, 50000.0, 5.0, 500.0, 50000.0, 5.0, 500.0, 50000.0, 5.0, 500.0]
    PB.ProblemData = _problem_data()
    PB.func_idx = _objective
    PB.func_F = lambda xr: _objective(_cell2idx(xr))
    PB.idx2cell = _idx2cell
    PB.cell2idx = _cell2idx
    PB.decode = lambda xr: _decode(_cell2idx(xr))
    return PB
