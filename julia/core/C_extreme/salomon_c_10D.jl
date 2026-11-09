# salomon_c_10D  DFO Benchmark Suite v2.1 core, Julia version (self-contained; Base only).
#SALOMON_C_10D  DFO Benchmark Suite v2.1.0 core: salomon, dimension 10, pattern C, strategy extreme.
#
# Copyright (C) 2026 J. F. A. Madeira
#
# Objective F(z) = f(T(z)), where f is the original function (embedded below, translated from the MATLAB original) and T acts on the finite
# coordinates only: uniform grid of 0 levels on the interval of the domain, piecewise-affine map with |T - Id| <=
# half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
# Global optimum of F on this domain: f* = 0 (E: exact by structure); the optimal solution is listed in
# reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
# Variables: x1:C x2:C x3:C x4:C x5:C x6:C x7:C x8:C x9:C x10:C
#
# Julia version of the MATLAB wrapper salomon_c_10D.m: same instance, same values (level indices 1..m as in MATLAB).
# Usage:  include("salomon_c_10D.jl");  PB = salomon_c_10D.problem();  f = PB.func_idx(PB.x0_idx)
# Copyright (C) 2026 J. F. A. Madeira. LGPL-3.0-or-later.
module salomon_c_10D

const NAME = "salomon_c_10D"
const N = 10
const TYPES = Char['C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C']
const CARD = Float64[Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf]
const LB = Float64[-100.0, -100.0, -63.779328, -70.41232819658475, -69.91407690731508, -62.55383214873429, -100.0, -100.0, -100.0, -100.0]
const UB = Float64[63.399169921875, 61.59325306101712, 100.0, 100.0, 100.0, 100.0, 69.5715808000383, 64.0719147336193, 70.75630804635489, 61.75997375866169]
const LBW = Float64[-100.0, -100.0, -63.77932799999999, -70.41232819658475, -69.91407690731506, -6255383214.873427, -10000000000.0, -10000000000.0, -10000000000.0, -10000000000.0]
const UBW = Float64[63.399169921875, 61.59325306101712, 100.0, 100.0, 100.0, 9999999999.999998, 6957158080.00383, 6407191473.361929, 7075630804.6354885, 6175997375.866169]
const TB = Float64[NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN]
const TI = Float64[NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN]
const CW = Float64[0.0, 0.0, 7.105427357601002e-15, 0.0, 1.4210854715202004e-14, 9.5367431640625e-07, 0.0, 0.0, -1.9073486328125e-06, 0.0]
const CX = Float64[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
const M = 0
const PERM = Dict{Int,Vector{Int}}()

function _linspace(a::Float64, b::Float64, m::Int)
    # MATLAB linspace: a + (k*(b-a))/(m-1), first and last values exactly a and b
    n1 = m - 1
    y = [a + (k * (b - a)) / n1 for k in 0:m-1]
    y[1] = a; y[end] = b
    return y
end

_err(msg) = throw(ArgumentError(NAME * ": " * msg))

function _vector(z)
    ok = (z isa AbstractVector || z isa Tuple) && length(z) == N && all(v -> v isa Real && !(v isa Bool) && isfinite(v), z)
    ok || _err("z must be a finite numeric vector with $N entries.")
    return Float64[v for v in z]
end

function _check_idx(v::Float64, m::Int, i::Int)
    (isfinite(v) && v == floor(v) && 1 <= v <= m) || _err("z$i = $v is not an index in 1..$m.")
    return Int(v)
end

function _orig(x::Vector{Float64})
    s = 0.0
    for v in x; s += v^2; end
    r = sqrt(s)
    return 1 - cos(2 * pi * r) + 0.1 * r
end

function _T_map(g, a, b, zs, xs)
    # piecewise-affine T with T(zs) = xs; fixes the bounds when zs is interior
    abs(g - zs) <= 1e-12 * max(1, abs(b - a)) && return xs
    if abs(zs - a) <= 1e-14 * max(1, abs(a))
        x = xs + (g - a) * (b - xs) / (b - a)
    elseif abs(zs - b) <= 1e-14 * max(1, abs(b))
        x = a + (g - a) * (xs - a) / (b - a)
    elseif g <= zs
        x = a + (g - a) * (xs - a) / (zs - a)
    else
        x = xs + (g - zs) * (b - xs) / (b - zs)
    end
    return min(max(x, min(a, xs)), max(b, xs))
end

function _map_work_to_orig(xw, lbw, ubw, lbo, ubo)
    t = abs(ubw - lbw) <= 1e-14 ? 0.0 : (xw - lbw) / (ubw - lbw)
    t = max(0.0, min(1.0, t))
    return lbo + t * (ubo - lbo)
end

const E = Dict{Int,Vector{Float64}}(i => _linspace(LB[i], UB[i], M)[p] for (i, p) in PERM)

function _decode(z)
    z = _vector(z)
    xh = zeros(N)
    for i in 1:N
        t = TYPES[i]
        if t == 'C'
            v = z[i]; tl = 1e-12 * (UBW[i] - LBW[i])
            (v < LBW[i] - tl || v > UBW[i] + tl) && _err("x$i = $v outside [$(LBW[i]), $(UBW[i])].")
            v = min(max(v, LBW[i]), UBW[i])
            xh[i] = abs(v - CW[i]) <= 1e-12 * (UBW[i] - LBW[i]) ? CX[i] : _map_work_to_orig(v, LBW[i], UBW[i], LB[i], UB[i])
        else
            k = _check_idx(z[i], M, i)
            g = t == 'K' ? E[i][k] : LB[i] + (k - 1) * (UB[i] - LB[i]) / (M - 1)
            xh[i] = _T_map(g, LB[i], UB[i], TB[i], TI[i])
        end
    end
    return xh
end

_objective(z) = _orig(_decode(z))

function _problem_data()
    pd = Any[]
    for i in 1:N
        t = TYPES[i]
        if t == 'C'; push!(pd, [LBW[i], UBW[i]])
        elseif t == 'D'; push!(pd, Any[_linspace(LB[i], UB[i], M)])
        else; push!(pd, ["k$(i)_$j" for j in 1:M]); end
    end
    return pd
end

function _idx2cell(z)
    z = _vector(z); pd = _problem_data(); xr = Any[]
    for i in 1:N
        t = TYPES[i]
        if t == 'C'; push!(xr, z[i])
        elseif t == 'D'; push!(xr, pd[i][1][_check_idx(z[i], Int(CARD[i]), i)])
        else; push!(xr, pd[i][_check_idx(z[i], Int(CARD[i]), i)]); end
    end
    return xr
end

function _cell2idx(xr)
    ((xr isa AbstractVector || xr isa Tuple) && length(xr) == N) || _err("xr must be a vector with $N entries.")
    pd = _problem_data(); z = zeros(N)
    for i in 1:N
        t = TYPES[i]; v = xr[i]
        if t == 'C'
            (v isa Real && !(v isa Bool) && isfinite(v)) || _err("x$i must be a finite scalar.")
            z[i] = Float64(v)
        elseif t == 'D'
            (v isa Real && !(v isa Bool) && isfinite(v)) || _err("x$i must be a grid value.")
            grid = pd[i][1]
            d, k = findmin(abs.(grid .- Float64(v)))
            d > 1e-10 * max(1, maximum(abs.(grid))) && _err("x$i = $v is not a grid value.")
            z[i] = k
        else
            k = v isa AbstractString ? findfirst(==(v), pd[i]) : nothing
            k === nothing && _err("x$i unknown label $v.")
            z[i] = k
        end
    end
    return z
end

"""
    problem()

Return the problem as a named tuple: `func_idx(z)` (level indices 1..m, continuous coordinates as values),
`func_F(xr)` (values, grid values, labels), `idx2cell`, `cell2idx`, `x0_idx`, `x0_alt_idx`, `f_star`, ...
"""
function problem()
    return (name = NAME, nvar = N, var_types = String[string(t) for t in TYPES],
        cardinality = Float64[Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf],
        f_star = 0.0,
        f_star_justification = "E",
        n_minimizers_documented = 1,
        description = "DFO Benchmark Suite v2.1.0 core (Julia): C adaptation of salomon (F = f o T).",
        x0_idx = Float64[27.73538625240326, 53.90889169003137, 56.29057294336002, -63.243796273279855, 78.55329457518488, -2267376135.757917, -3418013340.409669, 5412091363.591034, -3349107061.463995, 5734273063.580841],
        x0_alt_idx = Vector{Float64}[Float64[-52.68764269351959, 58.56384136668618, -18.39280062463999, -55.79136209558664, 25.31950789573304, -4287275825.1800585, -8170711625.808664, -6493228940.086225, -8447797512.128967, -7595056522.002102], Float64[-12.476128220558167, 47.25896358052452, 89.30848546816, -48.338927917893436, -12.595203480279402, 9852022000.774933, 2981044481.761532, -1175662728.4980087, 4237879059.1786804, -4767622973.5451145], Float64[-92.89915716648102, 50.14059909464416, 14.625111900160007, -64.7342831088185, -65.82899015973123, 6589107117.862243, -1771653803.6374664, 4141903483.090206, -119796098.90409851, -1940189425.0881271], Float64[53.58564555644989, 60.115491258904484, -36.15958212608, -57.28184893112528, 64.76612680208945, 4569207428.440103, -6524352089.036463, -6902272494.82378, -4477471256.986877, 887244123.3688602], Float64[-26.83738338947296, 48.8106134727428, 46.38519918592003, -49.829414753432076, 10.139696913234033, 1306292545.5274105, 5621432418.094309, -1584706283.2355633, -8835146415.06966, 3714677671.8258495], Float64[13.374131083488464, 53.46556314939761, -3.1416696012800003, -69.6315969970169, -27.775014462778415, -713607143.8947315, 868734132.6953125, 3732859928.3526516, 3850530156.2379894, -9614651913.757093], Float64[-67.04889786243439, 55.75609394267221, 79.40311171072003, -62.17916281932368, 87.50102719560289, -2733506833.316872, -4877992552.264261, -8172460375.324608, -507145001.8447914, -6787218365.300106], Float64[33.001179814338684, 44.4512161565105, 29.87624292352001, -54.72672864163047, 49.586315819590425, -5996421716.229563, -9630690837.66326, -2854894163.736394, -4864820159.92757, -3959784816.8431177], Float64[-47.42184913158417, 54.426108320770794, -44.80713064448, -47.27429446393726, -3.6474708598614, 8142876109.725429, 2515093669.467512, 2462672047.851822, -9963510610.592543, -1132351268.3861294]],
        suite_release = "2.1.0",
        base_function = "salomon",
        m = 0,
        pattern = "C",
        strategy = "extreme",
        scale_factors_C = Float64[1.0, 1.0, 1.0, 1.0, 1.0, 100000000.0, 100000000.0, 100000000.0, 100000000.0, 100000000.0],
        lb_orig = Float64[-100.0, -100.0, -63.779328, -70.41232819658475, -69.91407690731508, -62.55383214873429, -100.0, -100.0, -100.0, -100.0],
        ub_orig = Float64[63.399169921875, 61.59325306101712, 100.0, 100.0, 100.0, 100.0, 69.5715808000383, 64.0719147336193, 70.75630804635489, 61.75997375866169],
        lb_C_work = Float64[-100.0, -100.0, -63.77932799999999, -70.41232819658475, -69.91407690731506, -6255383214.873427, -10000000000.0, -10000000000.0, -10000000000.0, -10000000000.0],
        ub_C_work = Float64[63.399169921875, 61.59325306101712, 100.0, 100.0, 100.0, 9999999999.999998, 6957158080.00383, 6407191473.361929, 7075630804.6354885, 6175997375.866169],
        ProblemData = _problem_data(),
        func_idx = _objective,
        func_F = xr -> _objective(_cell2idx(xr)),
        idx2cell = _idx2cell,
        cell2idx = _cell2idx,
        decode = xr -> _decode(_cell2idx(xr)))
end

end # module
