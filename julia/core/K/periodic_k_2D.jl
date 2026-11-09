# periodic_k_2D  DFO Benchmark Suite v2.1 core, Julia version (self-contained; Base only).
#PERIODIC_K_2D  DFO Benchmark Suite v2.1.0 core: periodic, dimension 2, pattern K, strategy baseline.
#
# Copyright (C) 2026 J. F. A. Madeira
#
# Objective F(z) = f(T(z)), where f is the original function (embedded below, translated from the MATLAB original) and T acts on the finite
# coordinates only: uniform grid of 151 levels on the interval of the domain, piecewise-affine map with |T - Id| <=
# half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
# Global optimum of F on this domain: f* = 0.9 (E: exact by structure); the optimal solution is listed in
# reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
# Variables: x1:K x2:K
#
# Julia version of the MATLAB wrapper periodic_k_2D.m: same instance, same values (level indices 1..m as in MATLAB).
# Usage:  include("periodic_k_2D.jl");  PB = periodic_k_2D.problem();  f = PB.func_idx(PB.x0_idx)
# Copyright (C) 2026 J. F. A. Madeira. LGPL-3.0-or-later.
module periodic_k_2D

const NAME = "periodic_k_2D"
const N = 2
const TYPES = Char['K', 'K']
const CARD = Float64[151.0, 151.0]
const LB = Float64[-10.0, -10.0]
const UB = Float64[6.3563232421875, 6.055621602398009]
const LBW = Float64[-10.0, -10.0]
const UBW = Float64[6.3563232421875, 6.055621602398009]
const TB = Float64[0.03187825520833343, -0.04551460651323325]
const TI = Float64[0.0, 0.0]
const CW = Float64[NaN, NaN]
const CX = Float64[NaN, NaN]
const M = 151
const PERM = Dict{Int,Vector{Int}}(1 => Int[29, 55, 117, 132, 31, 140, 43, 25, 131, 78, 18, 8, 11, 27, 115, 52, 73, 108, 75, 86, 76, 145, 126, 97, 4, 41, 46, 17, 21, 96, 1, 67, 150, 65, 49, 15, 141, 34, 51, 72, 7, 13, 151, 23, 60, 125, 109, 101, 102, 22, 136, 47, 139, 69, 74, 116, 122, 148, 39, 147, 144, 91, 114, 82, 10, 28, 40, 118, 83, 128, 38, 146, 33, 134, 9, 120, 54, 95, 35, 77, 30, 53, 66, 44, 84, 12, 42, 58, 2, 32, 36, 56, 5, 70, 119, 100, 45, 26, 68, 99, 103, 110, 92, 59, 123, 24, 71, 61, 138, 130, 149, 129, 57, 87, 90, 111, 127, 104, 93, 79, 88, 14, 112, 48, 63, 20, 85, 81, 142, 113, 107, 19, 50, 143, 105, 37, 137, 3, 106, 62, 64, 80, 133, 98, 135, 89, 121, 94, 6, 124, 16], 2 => Int[69, 22, 44, 19, 112, 71, 49, 90, 109, 120, 28, 122, 126, 140, 143, 53, 6, 33, 151, 136, 91, 46, 20, 30, 145, 92, 13, 18, 21, 76, 147, 131, 74, 75, 108, 16, 24, 94, 55, 85, 8, 14, 17, 101, 115, 63, 150, 107, 57, 142, 102, 149, 25, 66, 51, 78, 105, 70, 41, 59, 27, 138, 43, 54, 111, 106, 42, 88, 83, 81, 40, 116, 35, 52, 11, 80, 82, 7, 144, 79, 137, 56, 146, 47, 141, 15, 45, 61, 5, 36, 132, 124, 9, 99, 2, 133, 50, 31, 73, 3, 119, 62, 128, 110, 86, 29, 125, 67, 10, 118, 12, 87, 64, 95, 121, 139, 93, 113, 103, 89, 98, 26, 135, 60, 77, 34, 100, 129, 97, 1, 127, 39, 72, 23, 130, 65, 32, 38, 4, 104, 114, 134, 58, 148, 84, 68, 123, 37, 48, 96, 117])

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

_orig(x::Vector{Float64}) = 1 + sin(x[1])^2 + sin(x[2])^2 - 0.1 * exp(-x[1]^2 - x[2]^2)

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
        f_star = 0.9,
        f_star_justification = "E",
        n_minimizers_documented = 1,
        description = "DFO Benchmark Suite v2.1.0 core (Julia): K adaptation of periodic (F = f o T).",
        x0_idx = Float64[104.0, 111.0],
        x0_alt_idx = Vector{Float64}[Float64[30.0, 115.0], Float64[124.0, 104.0], Float64[50.0, 107.0], Float64[86.0, 116.0], Float64[12.0, 106.0], Float64[133.0, 110.0], Float64[59.0, 112.0], Float64[95.0, 102.0], Float64[21.0, 111.0]],
        suite_release = "2.1.0",
        base_function = "periodic",
        m = 151,
        pattern = "K",
        strategy = "baseline",
        scale_factors_C = Float64[],
        lb_orig = Float64[-10.0, -10.0],
        ub_orig = Float64[6.3563232421875, 6.055621602398009],
        lb_C_work = Float64[-10.0, -10.0],
        ub_C_work = Float64[6.3563232421875, 6.055621602398009],
        ProblemData = _problem_data(),
        func_idx = _objective,
        func_F = xr -> _objective(_cell2idx(xr)),
        idx2cell = _idx2cell,
        cell2idx = _cell2idx,
        decode = xr -> _decode(_cell2idx(xr)))
end

end # module
