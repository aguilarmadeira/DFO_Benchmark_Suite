# threehumpcamel_k_2D  DFO Benchmark Suite v2.1 core, Julia version (self-contained; Base only).
#THREEHUMPCAMEL_K_2D  DFO Benchmark Suite v2.1.0 core: threehumpcamel, dimension 2, pattern K, strategy baseline.
#
# Copyright (C) 2026 J. F. A. Madeira
#
# Objective F(z) = f(T(z)), where f is the original function (embedded below, translated from the MATLAB original) and T acts on the finite
# coordinates only: uniform grid of 151 levels on the interval of the domain, piecewise-affine map with |T - Id| <=
# half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
# Global optimum of F on this domain: f* = 0 (E: exact by structure); the optimal solution is listed in
# reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
# Variables: x1:K x2:K
#
# Julia version of the MATLAB wrapper threehumpcamel_k_2D.m: same instance, same values (level indices 1..m as in MATLAB).
# Usage:  include("threehumpcamel_k_2D.jl");  PB = threehumpcamel_k_2D.problem();  f = PB.func_idx(PB.x0_idx)
# Copyright (C) 2026 J. F. A. Madeira. LGPL-3.0-or-later.
module threehumpcamel_k_2D

const NAME = "threehumpcamel_k_2D"
const N = 2
const TYPES = Char['K', 'K']
const CARD = Float64[151.0, 151.0]
const LB = Float64[-5.0, -3.1276736269877556]
const UB = Float64[3.23558349609375, 5.0]
const LBW = Float64[-5.0, -3.1276736269877556]
const UBW = Float64[3.23558349609375, 5.0]
const TB = Float64[-0.0037460123697918846, 0.015026842114176375]
const TI = Float64[0.0, 0.0]
const CW = Float64[NaN, NaN]
const CX = Float64[NaN, NaN]
const M = 151
const PERM = Dict{Int,Vector{Int}}(1 => Int[32, 37, 9, 14, 81, 139, 124, 129, 90, 55, 147, 67, 89, 132, 3, 4, 47, 30, 142, 58, 40, 112, 128, 29, 149, 16, 134, 91, 5, 110, 115, 46, 79, 97, 99, 85, 12, 61, 63, 34, 104, 151, 57, 100, 136, 116, 141, 10, 21, 94, 113, 11, 48, 96, 64, 76, 72, 103, 101, 49, 69, 38, 66, 59, 138, 130, 95, 17, 144, 108, 111, 52, 19, 39, 121, 83, 65, 109, 62, 105, 78, 126, 86, 36, 2, 150, 80, 125, 54, 27, 56, 82, 43, 135, 60, 106, 107, 140, 35, 137, 119, 118, 133, 1, 145, 146, 120, 68, 148, 42, 114, 20, 50, 6, 70, 44, 41, 74, 131, 51, 7, 88, 33, 18, 13, 26, 84, 102, 87, 45, 92, 122, 15, 22, 31, 28, 73, 143, 77, 53, 75, 123, 23, 25, 127, 98, 71, 24, 93, 117, 8], 2 => Int[19, 145, 40, 14, 77, 97, 144, 70, 117, 80, 65, 137, 75, 133, 36, 4, 7, 149, 33, 15, 98, 73, 2, 129, 48, 88, 23, 94, 85, 47, 68, 87, 25, 90, 41, 3, 111, 61, 130, 34, 1, 9, 116, 106, 82, 45, 71, 11, 22, 32, 119, 12, 110, 147, 59, 108, 53, 121, 28, 51, 122, 84, 66, 91, 24, 93, 120, 18, 151, 69, 8, 54, 21, 42, 132, 43, 150, 141, 140, 38, 81, 29, 143, 127, 5, 44, 128, 16, 101, 124, 60, 86, 123, 35, 64, 107, 55, 113, 126, 115, 96, 79, 20, 6, 89, 103, 136, 74, 142, 49, 139, 27, 138, 13, 78, 52, 50, 83, 102, 62, 17, 100, 134, 31, 26, 39, 99, 118, 104, 63, 148, 10, 37, 46, 56, 57, 105, 58, 114, 95, 125, 30, 92, 109, 146, 76, 72, 112, 135, 67, 131])

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

_orig(x::Vector{Float64}) = 2 * x[1]^2 - 1.05 * x[1]^4 + 1 / 6 * x[1]^6 + x[1] * x[2] + x[2]^2

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
        f_star = 0.0,
        f_star_justification = "E",
        n_minimizers_documented = 1,
        description = "DFO Benchmark Suite v2.1.0 core (Julia): K adaptation of threehumpcamel (F = f o T).",
        x0_idx = Float64[97.0, 135.0],
        x0_alt_idx = Vector{Float64}[Float64[23.0, 144.0], Float64[117.0, 149.0], Float64[43.0, 137.0], Float64[79.0, 141.0], Float64[5.0, 150.0], Float64[147.0, 140.0], Float64[73.0, 142.0], Float64[109.0, 146.0], Float64[35.0, 136.0]],
        suite_release = "2.1.0",
        base_function = "threehumpcamel",
        m = 151,
        pattern = "K",
        strategy = "baseline",
        scale_factors_C = Float64[],
        lb_orig = Float64[-5.0, -3.1276736269877556],
        ub_orig = Float64[3.23558349609375, 5.0],
        lb_C_work = Float64[-5.0, -3.1276736269877556],
        ub_C_work = Float64[3.23558349609375, 5.0],
        ProblemData = _problem_data(),
        func_idx = _objective,
        func_F = xr -> _objective(_cell2idx(xr)),
        idx2cell = _idx2cell,
        cell2idx = _cell2idx,
        decode = xr -> _decode(_cell2idx(xr)))
end

end # module
