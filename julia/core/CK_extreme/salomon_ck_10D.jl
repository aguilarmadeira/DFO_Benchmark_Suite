# salomon_ck_10D  DFO Benchmark Suite v2.1 core, Julia version (self-contained; Base only).
#SALOMON_CK_10D  DFO Benchmark Suite v2.1.0 core: salomon, dimension 10, pattern CK, strategy extreme.
#
# Copyright (C) 2026 J. F. A. Madeira
#
# Objective F(z) = f(T(z)), where f is the original function (embedded below, translated from the MATLAB original) and T acts on the finite
# coordinates only: uniform grid of 17 levels on the interval of the domain, piecewise-affine map with |T - Id| <=
# half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
# Global optimum of F on this domain: f* = 0 (E: exact by structure); the optimal solution is listed in
# reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
# Variables: x1:C x2:K x3:C x4:K x5:C x6:K x7:C x8:K x9:C x10:K
#
# Julia version of the MATLAB wrapper salomon_ck_10D.m: same instance, same values (level indices 1..m as in MATLAB).
# Usage:  include("salomon_ck_10D.jl");  PB = salomon_ck_10D.problem();  f = PB.func_idx(PB.x0_idx)
# Copyright (C) 2026 J. F. A. Madeira. LGPL-3.0-or-later.
module salomon_ck_10D

const NAME = "salomon_ck_10D"
const N = 10
const TYPES = Char['C', 'K', 'C', 'K', 'C', 'K', 'C', 'K', 'C', 'K']
const CARD = Float64[Inf, 17.0, Inf, 17.0, Inf, 17.0, Inf, 17.0, Inf, 17.0]
const LB = Float64[-100.0, -100.0, -63.779328, -70.41232819658475, -69.91407690731508, -62.55383214873429, -100.0, -100.0, -100.0, -100.0]
const UB = Float64[63.399169921875, 61.59325306101712, 100.0, 100.0, 100.0, 100.0, 69.5715808000383, 64.0719147336193, 70.75630804635489, 61.75997375866169]
const LBW = Float64[-100.0, -100.0, -63.77932799999999, -70.41232819658475, -6991407690.731507, -62.55383214873429, -10000000000.0, -100.0, -10000000000.0, -100.0]
const UBW = Float64[63.399169921875, 61.59325306101712, 100.0, 100.0, 10000000000.0, 100.0, 6957158080.00383, 64.0719147336193, 7075630804.6354885, 61.75997375866169]
const TB = Float64[NaN, 0.9957831631356981, NaN, 4.143065389421082, NaN, -1.5961450929589347, NaN, 2.544946708512043, NaN, 1.0999835991635507]
const TI = Float64[NaN, 0.0, NaN, 0.0, NaN, 0.0, NaN, 0.0, NaN, 0.0]
const CW = Float64[0.0, NaN, 7.105427357601002e-15, NaN, 0.0, NaN, 0.0, NaN, -1.9073486328125e-06, NaN]
const CX = Float64[0.0, NaN, 0.0, NaN, 0.0, NaN, 0.0, NaN, 0.0, NaN]
const M = 17
const PERM = Dict{Int,Vector{Int}}(2 => Int[14, 8, 5, 10, 2, 15, 4, 1, 7, 13, 16, 3, 17, 6, 12, 9, 11], 4 => Int[16, 5, 6, 8, 4, 13, 12, 9, 14, 1, 2, 11, 10, 7, 3, 15, 17], 6 => Int[3, 1, 9, 4, 13, 17, 7, 11, 6, 16, 5, 15, 14, 10, 8, 2, 12], 8 => Int[17, 10, 1, 9, 7, 2, 5, 3, 6, 14, 8, 4, 16, 12, 13, 11, 15], 10 => Int[7, 8, 17, 15, 2, 10, 3, 1, 5, 6, 14, 12, 11, 16, 4, 13, 9])

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
        cardinality = Float64[Inf, 17.0, Inf, 17.0, Inf, 17.0, Inf, 17.0, Inf, 17.0],
        f_star = 0.0,
        f_star_justification = "E",
        n_minimizers_documented = 1,
        description = "DFO Benchmark Suite v2.1.0 core (Julia): CK adaptation of salomon (F = f o T).",
        x0_idx = Float64[27.73538625240326, 16.0, 56.29057294336002, 2.0, 7855329457.518487, 5.0, -3418013340.409669, 16.0, -3349107061.463995, 17.0],
        x0_alt_idx = Vector{Float64}[Float64[-52.68764269351959, 17.0, -18.39280062463999, 2.0, 2531950789.5733023, 3.0, -8170711625.808664, 4.0, -8447797512.128967, 3.0], Float64[-12.476128220558167, 16.0, 89.30848546816, 3.0, -1259520348.0279417, 17.0, 2981044481.761532, 10.0, 4237879059.1786804, 6.0], Float64[-92.89915716648102, 16.0, 14.625111900160007, 2.0, -6582899015.973124, 14.0, -1771653803.6374664, 15.0, -119796098.90409851, 9.0], Float64[53.58564555644989, 17.0, -36.15958212608, 2.0, 6476612680.208945, 12.0, -6524352089.036463, 4.0, -4477471256.986877, 12.0], Float64[-26.83738338947296, 16.0, 46.38519918592003, 3.0, 1013969691.3234024, 8.0, 5621432418.094309, 9.0, -8835146415.06966, 15.0], Float64[13.374131083488464, 16.0, -3.1416696012800003, 1.0, -2777501446.277843, 6.0, 868734132.6953125, 14.0, 3850530156.2379894, 1.0], Float64[-67.04889786243439, 16.0, 79.40311171072003, 2.0, 8750102719.560287, 4.0, -4877992552.264261, 3.0, -507145001.8447914, 4.0], Float64[33.001179814338684, 15.0, 29.87624292352001, 2.0, 4958631581.959041, 1.0, -9630690837.66326, 8.0, -4864820159.92757, 7.0], Float64[-47.42184913158417, 16.0, -44.80713064448, 3.0, -364747085.9861412, 15.0, 2515093669.467512, 13.0, -9963510610.592543, 10.0]],
        suite_release = "2.1.0",
        base_function = "salomon",
        m = 17,
        pattern = "CK",
        strategy = "extreme",
        scale_factors_C = Float64[1.0, 1.0, 100000000.0, 100000000.0, 100000000.0],
        lb_orig = Float64[-100.0, -100.0, -63.779328, -70.41232819658475, -69.91407690731508, -62.55383214873429, -100.0, -100.0, -100.0, -100.0],
        ub_orig = Float64[63.399169921875, 61.59325306101712, 100.0, 100.0, 100.0, 100.0, 69.5715808000383, 64.0719147336193, 70.75630804635489, 61.75997375866169],
        lb_C_work = Float64[-100.0, -100.0, -63.77932799999999, -70.41232819658475, -6991407690.731507, -62.55383214873429, -10000000000.0, -100.0, -10000000000.0, -100.0],
        ub_C_work = Float64[63.399169921875, 61.59325306101712, 100.0, 100.0, 10000000000.0, 100.0, 6957158080.00383, 64.0719147336193, 7075630804.6354885, 61.75997375866169],
        ProblemData = _problem_data(),
        func_idx = _objective,
        func_F = xr -> _objective(_cell2idx(xr)),
        idx2cell = _idx2cell,
        cell2idx = _cell2idx,
        decode = xr -> _decode(_cell2idx(xr)))
end

end # module
