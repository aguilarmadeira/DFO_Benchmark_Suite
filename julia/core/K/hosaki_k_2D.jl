# hosaki_k_2D  DFO Benchmark Suite v2.1 core, Julia version (self-contained; Base only).
#HOSAKI_K_2D  DFO Benchmark Suite v2.1.0 core: hosaki, dimension 2, pattern K, strategy baseline.
#
# Copyright (C) 2026 J. F. A. Madeira
#
# Objective F(z) = f(T(z)), where f is the original function (embedded below, translated from the MATLAB original) and T acts on the finite
# coordinates only: uniform grid of 151 levels on the interval of the domain, piecewise-affine map with |T - Id| <=
# half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
# Global optimum of F on this domain: f* = -2.345811576101292 (L: literature value at a refined point); the optimal solution is listed in
# reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
# Variables: x1:K x2:K
#
# Julia version of the MATLAB wrapper hosaki_k_2D.m: same instance, same values (level indices 1..m as in MATLAB).
# Usage:  include("hosaki_k_2D.jl");  PB = hosaki_k_2D.problem();  f = PB.func_idx(PB.x0_idx)
# Copyright (C) 2026 J. F. A. Madeira. LGPL-3.0-or-later.
module hosaki_k_2D

const NAME = "hosaki_k_2D"
const N = 2
const TYPES = Char['K', 'K']
const CARD = Float64[151.0, 151.0]
const LB = Float64[0.0, 1.027757963725042]
const UB = Float64[4.067205810546875, 6.0]
const LBW = Float64[0.0, 1.027757963725042]
const UBW = Float64[4.067205810546875, 6.0]
const TB = Float64[4.012976399739584, 1.9890580907382005]
const TI = Float64[4.0, 2.0]
const CW = Float64[NaN, NaN]
const CX = Float64[NaN, NaN]
const M = 151
const PERM = Dict{Int,Vector{Int}}(1 => Int[97, 61, 111, 68, 42, 81, 38, 59, 139, 73, 100, 117, 29, 22, 14, 41, 145, 34, 118, 52, 149, 60, 112, 135, 72, 113, 96, 127, 8, 46, 57, 19, 143, 27, 31, 87, 150, 13, 99, 75, 24, 116, 53, 39, 103, 130, 105, 20, 28, 84, 131, 86, 108, 107, 55, 121, 74, 32, 98, 56, 80, 129, 128, 45, 122, 89, 26, 43, 49, 9, 92, 146, 3, 17, 91, 82, 125, 36, 132, 21, 50, 71, 136, 54, 5, 6, 102, 47, 101, 4, 88, 94, 51, 15, 35, 2, 67, 124, 66, 114, 76, 140, 64, 141, 110, 37, 79, 58, 104, 18, 115, 134, 78, 93, 126, 138, 25, 106, 83, 30, 10, 65, 1, 40, 11, 70, 144, 147, 120, 77, 151, 63, 119, 12, 44, 133, 48, 62, 142, 33, 23, 148, 109, 85, 123, 137, 95, 16, 7, 90, 69], 2 => Int[61, 124, 140, 26, 42, 66, 128, 35, 150, 120, 63, 59, 48, 22, 118, 68, 45, 27, 4, 103, 107, 10, 46, 16, 62, 106, 79, 147, 113, 56, 108, 20, 110, 67, 31, 117, 76, 123, 15, 82, 136, 30, 114, 130, 119, 12, 126, 121, 29, 54, 72, 7, 60, 87, 91, 99, 151, 33, 41, 58, 131, 32, 37, 47, 132, 138, 148, 127, 88, 11, 17, 50, 5, 94, 149, 141, 13, 102, 65, 105, 139, 74, 2, 57, 90, 9, 134, 51, 80, 129, 109, 18, 55, 112, 39, 6, 116, 44, 135, 98, 137, 38, 142, 28, 23, 43, 85, 64, 1, 125, 122, 19, 86, 101, 78, 71, 143, 146, 93, 40, 21, 77, 14, 52, 144, 84, 100, 111, 133, 95, 89, 83, 8, 36, 70, 25, 81, 97, 53, 75, 145, 69, 34, 24, 73, 115, 96, 104, 49, 3, 92])

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

_orig(x::Vector{Float64}) = (1 - 8 * x[1] + 7 * x[1]^2 - 7 / 3 * x[1]^3 + 1 / 4 * x[1]^4) * x[2]^2 * exp(-x[2])

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
        f_star = -2.345811576101292,
        f_star_justification = "L",
        n_minimizers_documented = 1,
        description = "DFO Benchmark Suite v2.1.0 core (Julia): K adaptation of hosaki (F = f o T).",
        x0_idx = Float64[114.0, 55.0],
        x0_alt_idx = Vector{Float64}[Float64[40.0, 59.0], Float64[77.0, 64.0], Float64[3.0, 51.0], Float64[150.0, 61.0], Float64[74.0, 65.0], Float64[113.0, 55.0], Float64[37.0, 57.0], Float64[132.0, 66.0], Float64[55.0, 55.0]],
        suite_release = "2.1.0",
        base_function = "hosaki",
        m = 151,
        pattern = "K",
        strategy = "baseline",
        scale_factors_C = Float64[],
        lb_orig = Float64[0.0, 1.027757963725042],
        ub_orig = Float64[4.067205810546875, 6.0],
        lb_C_work = Float64[0.0, 1.027757963725042],
        ub_C_work = Float64[4.067205810546875, 6.0],
        ProblemData = _problem_data(),
        func_idx = _objective,
        func_F = xr -> _objective(_cell2idx(xr)),
        idx2cell = _idx2cell,
        cell2idx = _cell2idx,
        decode = xr -> _decode(_cell2idx(xr)))
end

end # module
