# goldstein_price_k_2D  DFO Benchmark Suite v2.1 core, Julia version (self-contained; Base only).
#GOLDSTEIN_PRICE_K_2D  DFO Benchmark Suite v2.1.0 core: goldstein_price, dimension 2, pattern K, strategy baseline.
#
# Copyright (C) 2026 J. F. A. Madeira
#
# Objective F(z) = f(T(z)), where f is the original function (embedded below, translated from the MATLAB original) and T acts on the finite
# coordinates only: uniform grid of 151 levels on the interval of the domain, piecewise-affine map with |T - Id| <=
# half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
# Global optimum of F on this domain: f* = 3 (L: literature value at a refined point); the optimal solution is listed in
# reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
# Variables: x1:K x2:K
#
# Julia version of the MATLAB wrapper goldstein_price_k_2D.m: same instance, same values (level indices 1..m as in MATLAB).
# Usage:  include("goldstein_price_k_2D.jl");  PB = goldstein_price_k_2D.problem();  f = PB.func_idx(PB.x0_idx)
# Copyright (C) 2026 J. F. A. Madeira. LGPL-3.0-or-later.
module goldstein_price_k_2D

const NAME = "goldstein_price_k_2D"
const N = 2
const TYPES = Char['K', 'K']
const CARD = Float64[151.0, 151.0]
const LB = Float64[-2.0, -2.0]
const UB = Float64[1.2800146484375, 1.3770502464055276]
const LBW = Float64[-2.0, -2.0]
const UBW = Float64[1.2800146484375, 1.3770502464055276]
const TB = Float64[-0.01012444661458356, -1.009398594387712]
const TI = Float64[0.0, -1.0]
const CW = Float64[NaN, NaN]
const CX = Float64[NaN, NaN]
const M = 151
const PERM = Dict{Int,Vector{Int}}(1 => Int[13, 35, 137, 106, 30, 32, 67, 44, 140, 18, 48, 71, 151, 22, 41, 3, 90, 60, 108, 102, 129, 136, 123, 40, 113, 78, 66, 134, 12, 38, 58, 23, 133, 93, 34, 117, 142, 15, 125, 95, 49, 76, 126, 39, 74, 10, 72, 63, 53, 33, 143, 2, 91, 130, 47, 94, 92, 145, 147, 131, 149, 14, 68, 114, 107, 121, 4, 19, 139, 148, 88, 81, 105, 56, 9, 146, 87, 70, 17, 150, 116, 16, 101, 96, 124, 21, 122, 59, 86, 5, 83, 98, 120, 1, 20, 77, 111, 8, 25, 50, 138, 80, 6, 64, 43, 45, 112, 57, 99, 104, 42, 118, 24, 28, 103, 79, 36, 84, 54, 97, 61, 46, 82, 109, 55, 100, 85, 26, 119, 29, 73, 37, 27, 75, 62, 135, 132, 128, 127, 69, 52, 7, 65, 141, 89, 31, 110, 144, 115, 51, 11], 2 => Int[139, 64, 37, 126, 70, 145, 4, 149, 141, 39, 72, 132, 44, 43, 41, 52, 94, 76, 121, 118, 92, 14, 17, 55, 112, 60, 67, 10, 98, 15, 74, 123, 77, 113, 35, 110, 129, 148, 106, 125, 138, 32, 50, 40, 56, 75, 18, 108, 53, 34, 65, 3, 104, 59, 109, 2, 151, 22, 101, 78, 8, 16, 20, 83, 130, 116, 6, 21, 136, 84, 142, 79, 117, 124, 103, 33, 135, 73, 82, 127, 96, 19, 134, 144, 28, 24, 115, 62, 89, 146, 86, 12, 97, 5, 150, 81, 46, 13, 29, 54, 30, 85, 11, 69, 48, 51, 90, 119, 105, 1, 143, 91, 31, 36, 140, 88, 45, 131, 63, 107, 71, 58, 95, 122, 68, 114, 100, 42, 7, 47, 93, 57, 49, 99, 87, 137, 27, 26, 25, 147, 102, 61, 133, 80, 38, 9, 111, 66, 128, 23, 120])

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
    f = 1 + (x[1] + x[2] + 1)^2 * (19 - 14 * x[1] + 3 * x[1]^2 - 14 * x[2] + 6 * x[1] * x[2] + 3 * x[2]^2)
    return f * (30 + (2 * x[1] - 3 * x[2])^2 * (18 - 32 * x[1] + 12 * x[1]^2 + 48 * x[2] - 36 * x[1] * x[2] + 27 * x[2]^2))
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
        cardinality = Float64[151.0, 151.0],
        f_star = 3.0,
        f_star_justification = "L",
        n_minimizers_documented = 1,
        description = "DFO Benchmark Suite v2.1.0 core (Julia): K adaptation of goldstein_price (F = f o T).",
        x0_idx = Float64[113.0, 81.0],
        x0_alt_idx = Vector{Float64}[Float64[36.0, 70.0], Float64[131.0, 75.0], Float64[56.0, 84.0], Float64[94.0, 72.0], Float64[18.0, 76.0], Float64[141.0, 80.0], Float64[65.0, 68.0], Float64[104.0, 77.0], Float64[27.0, 81.0]],
        suite_release = "2.1.0",
        base_function = "goldstein_price",
        m = 151,
        pattern = "K",
        strategy = "baseline",
        scale_factors_C = Float64[],
        lb_orig = Float64[-2.0, -2.0],
        ub_orig = Float64[1.2800146484375, 1.3770502464055276],
        lb_C_work = Float64[-2.0, -2.0],
        ub_C_work = Float64[1.2800146484375, 1.3770502464055276],
        ProblemData = _problem_data(),
        func_idx = _objective,
        func_F = xr -> _objective(_cell2idx(xr)),
        idx2cell = _idx2cell,
        cell2idx = _cell2idx,
        decode = xr -> _decode(_cell2idx(xr)))
end

end # module
