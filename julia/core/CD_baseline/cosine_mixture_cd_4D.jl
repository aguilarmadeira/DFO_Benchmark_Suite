# cosine_mixture_cd_4D  DFO Benchmark Suite v2.1 core, Julia version (self-contained; Base only).
#COSINE_MIXTURE_CD_4D  DFO Benchmark Suite v2.1.0 core: cosine_mixture, dimension 4, pattern CD, strategy baseline.
#
# Copyright (C) 2026 J. F. A. Madeira
#
# Objective F(z) = f(T(z)), where f is the original function (embedded below, translated from the MATLAB original) and T acts on the finite
# coordinates only: uniform grid of 151 levels on the interval of the domain, piecewise-affine map with |T - Id| <=
# half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
# Global optimum of F on this domain: f* = -4.4 (E: exact by structure); the optimal solution is listed in
# reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
# Variables: x1:C x2:D x3:C x4:D
#
# Julia version of the MATLAB wrapper cosine_mixture_cd_4D.m: same instance, same values (level indices 1..m as in MATLAB).
# Usage:  include("cosine_mixture_cd_4D.jl");  PB = cosine_mixture_cd_4D.problem();  f = PB.func_idx(PB.x0_idx)
# Copyright (C) 2026 J. F. A. Madeira. LGPL-3.0-or-later.
module cosine_mixture_cd_4D

const NAME = "cosine_mixture_cd_4D"
const N = 4
const TYPES = Char['C', 'D', 'C', 'D']
const CARD = Float64[Inf, 151.0, Inf, 151.0]
const LB = Float64[-0.66625732421875, -0.653957221968196, -0.62211328, -0.6984089962515618]
const UB = Float64[1.0, 1.0, 1.0, 1.0]
const LBW = Float64[-0.66625732421875, -0.653957221968196, -0.62211328, -0.6984089962515618]
const UBW = Float64[1.0, 1.0, 1.0, 1.0]
const TB = Float64[NaN, 1.0, NaN, 1.0]
const TI = Float64[NaN, 1.0, NaN, 1.0]
const CW = Float64[1.0, NaN, 1.0, NaN]
const CX = Float64[1.0, NaN, 1.0, NaN]
const M = 151
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
    s1 = 0.0
    for v in x; s1 += cos(5 * pi * v); end
    s2 = 0.0
    for v in x; s2 += v^2; end
    return 0.1 * s1 - s2
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
        cardinality = Float64[Inf, 151.0, Inf, 151.0],
        f_star = -4.4,
        f_star_justification = "E",
        n_minimizers_documented = 1,
        description = "DFO Benchmark Suite v2.1.0 core (Julia): CD adaptation of cosine_mixture (F = f o T).",
        x0_idx = Float64[0.26938521623611456, 59.0, 0.7487670951936001, 37.0],
        x0_alt_idx = Vector{Float64}[Float64[-0.5873379099369048, 63.0, 0.3205291892736002, 41.0], Float64[0.9357254254817962, 53.0, -0.4814436163584, 27.0], Float64[0.09608794569969181, 62.0, 0.6475472265216002, 33.0], Float64[0.5256699120998383, 64.0, -0.16688340910080002, 40.0], Float64[-0.3237307941913604, 54.0, 0.9621074337792002, 25.0], Float64[0.7306976687908173, 58.0, 0.16013462814720003, 29.0], Float64[-0.10893981099128724, 68.0, -0.5795490275328, 35.0], Float64[0.3206421554088593, 55.0, 0.4871526653952001, 42.0], Float64[-0.5336401641368865, 59.0, -0.25253099028479997, 27.0]],
        suite_release = "2.1.0",
        base_function = "cosine_mixture",
        m = 151,
        pattern = "CD",
        strategy = "baseline",
        scale_factors_C = Float64[1.0, 1.0],
        lb_orig = Float64[-0.66625732421875, -0.653957221968196, -0.62211328, -0.6984089962515618],
        ub_orig = Float64[1.0, 1.0, 1.0, 1.0],
        lb_C_work = Float64[-0.66625732421875, -0.653957221968196, -0.62211328, -0.6984089962515618],
        ub_C_work = Float64[1.0, 1.0, 1.0, 1.0],
        ProblemData = _problem_data(),
        func_idx = _objective,
        func_F = xr -> _objective(_cell2idx(xr)),
        idx2cell = _idx2cell,
        cell2idx = _cell2idx,
        decode = xr -> _decode(_cell2idx(xr)))
end

end # module
