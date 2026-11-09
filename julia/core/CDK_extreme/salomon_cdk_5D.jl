# salomon_cdk_5D  DFO Benchmark Suite v2.1 core, Julia version (self-contained; Base only).
#SALOMON_CDK_5D  DFO Benchmark Suite v2.1.0 core: salomon, dimension 5, pattern CDK, strategy extreme.
#
# Copyright (C) 2026 J. F. A. Madeira
#
# Objective F(z) = f(T(z)), where f is the original function (embedded below, translated from the MATLAB original) and T acts on the finite
# coordinates only: uniform grid of 41 levels on the interval of the domain, piecewise-affine map with |T - Id| <=
# half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
# Global optimum of F on this domain: f* = 0 (E: exact by structure); the optimal solution is listed in
# reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
# Variables: x1:C x2:K x3:D x4:C x5:D
#
# Julia version of the MATLAB wrapper salomon_cdk_5D.m: same instance, same values (level indices 1..m as in MATLAB).
# Usage:  include("salomon_cdk_5D.jl");  PB = salomon_cdk_5D.problem();  f = PB.func_idx(PB.x0_idx)
# Copyright (C) 2026 J. F. A. Madeira. LGPL-3.0-or-later.
module salomon_cdk_5D

const NAME = "salomon_cdk_5D"
const N = 5
const TYPES = Char['C', 'K', 'D', 'C', 'D']
const CARD = Float64[Inf, 41.0, 41.0, Inf, 41.0]
const LB = Float64[-70.399169921875, -100.0, -66.579328, -100.0, -71.18680418004234]
const UB = Float64[100.0, 66.2599197276838, 100.0, 72.41232819658475, 100.0]
const LBW = Float64[-70.399169921875, -100.0, -66.579328, -10000000000.0, -71.18680418004234]
const UBW = Float64[100.0, 66.2599197276838, 100.0, 7241232819.658475, 100.0]
const TB = Float64[NaN, -0.2440481633897349, 0.05240320000000054, NaN, 1.5675875964756472]
const TI = Float64[NaN, 0.0, 0.0, NaN, 0.0]
const CW = Float64[0.0, NaN, NaN, 0.0, NaN]
const CX = Float64[0.0, NaN, NaN, 0.0, NaN]
const M = 41
const PERM = Dict{Int,Vector{Int}}(2 => Int[11, 38, 21, 8, 10, 36, 29, 9, 23, 27, 19, 7, 6, 1, 25, 2, 31, 37, 12, 32, 40, 24, 18, 41, 5, 33, 15, 39, 30, 34, 14, 22, 17, 26, 4, 28, 35, 20, 16, 13, 3])

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
        cardinality = Float64[Inf, 41.0, 41.0, Inf, 41.0],
        f_star = 0.0,
        f_star_justification = "E",
        n_minimizers_documented = 1,
        description = "DFO Benchmark Suite v2.1.0 core (Julia): CDK adaptation of salomon (F = f o T).",
        x0_idx = Float64[-22.39120066165924, 26.0, 22.0, 5151604018.941849, 32.0],
        x0_alt_idx = Vector{Float64}[Float64[21.539835333824158, 27.0, 4.0, 5905593792.396301, 20.0], Float64[-64.32537138462067, 24.0, 30.0, 6659583565.850754, 11.0], Float64[91.4301198720932, 25.0, 12.0, 5000806064.2509575, 38.0], Float64[4.56648051738739, 27.0, 38.0, 5754795837.705412, 29.0], Float64[48.49751651287079, 25.0, 20.0, 6508785611.1598625, 16.0], Float64[-37.367690205574036, 26.0, 8.0, 6918094345.32085, 7.0], Float64[70.46303451061249, 26.0, 28.0, 5259316843.721054, 34.0], Float64[-16.899821162223816, 23.0, 16.0, 6013306617.175507, 25.0], Float64[27.031214833259583, 26.0, 36.0, 6767296390.629959, 13.0]],
        suite_release = "2.1.0",
        base_function = "salomon",
        m = 41,
        pattern = "CDK",
        strategy = "extreme",
        scale_factors_C = Float64[1.0, 100000000.0],
        lb_orig = Float64[-70.399169921875, -100.0, -66.579328, -100.0, -71.18680418004234],
        ub_orig = Float64[100.0, 66.2599197276838, 100.0, 72.41232819658475, 100.0],
        lb_C_work = Float64[-70.399169921875, -100.0, -66.579328, -10000000000.0, -71.18680418004234],
        ub_C_work = Float64[100.0, 66.2599197276838, 100.0, 7241232819.658475, 100.0],
        ProblemData = _problem_data(),
        func_idx = _objective,
        func_F = xr -> _objective(_cell2idx(xr)),
        idx2cell = _idx2cell,
        cell2idx = _cell2idx,
        decode = xr -> _decode(_cell2idx(xr)))
end

end # module
