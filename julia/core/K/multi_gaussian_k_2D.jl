# multi_gaussian_k_2D  DFO Benchmark Suite v2.1 core, Julia version (self-contained; Base only).
#MULTI_GAUSSIAN_K_2D  DFO Benchmark Suite v2.1.0 core: multi_gaussian, dimension 2, pattern K, strategy baseline.
#
# Copyright (C) 2026 J. F. A. Madeira
#
# Objective F(z) = f(T(z)), where f is the original function (embedded below, translated from the MATLAB original) and T acts on the finite
# coordinates only: uniform grid of 151 levels on the interval of the domain, piecewise-affine map with |T - Id| <=
# half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
# Global optimum of F on this domain: f* = -1.296954045953779 (L: literature value at a refined point); the optimal solution is listed in
# reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
# Variables: x1:K x2:K
#
# Julia version of the MATLAB wrapper multi_gaussian_k_2D.m: same instance, same values (level indices 1..m as in MATLAB).
# Usage:  include("multi_gaussian_k_2D.jl");  PB = multi_gaussian_k_2D.problem();  f = PB.func_idx(PB.x0_idx)
# Copyright (C) 2026 J. F. A. Madeira. LGPL-3.0-or-later.
module multi_gaussian_k_2D

const NAME = "multi_gaussian_k_2D"
const N = 2
const TYPES = Char['K', 'K']
const CARD = Float64[151.0, 151.0]
const LB = Float64[-1.3062646484375, -1.3355687649240462]
const UB = Float64[2.0, 2.0]
const LBW = Float64[-1.3062646484375, -1.3355687649240462]
const UBW = Float64[2.0, 2.0]
const TB = Float64[-0.005800553385416851, -0.02357838405392143]
const TI = Float64[-0.01354066395522768, -0.01354066351404425]
const CW = Float64[NaN, NaN]
const CX = Float64[NaN, NaN]
const M = 151
const PERM = Dict{Int,Vector{Int}}(1 => Int[101, 72, 135, 132, 39, 110, 17, 55, 146, 54, 12, 33, 49, 23, 64, 73, 58, 126, 47, 104, 90, 5, 105, 82, 32, 120, 149, 100, 127, 7, 117, 15, 119, 151, 57, 79, 59, 62, 145, 16, 121, 28, 134, 50, 116, 11, 18, 87, 141, 112, 1, 118, 38, 42, 99, 65, 27, 111, 107, 75, 68, 46, 77, 63, 103, 80, 48, 136, 4, 138, 123, 6, 128, 95, 9, 34, 52, 91, 30, 150, 96, 114, 56, 31, 69, 133, 125, 115, 67, 3, 2, 21, 61, 106, 51, 22, 93, 71, 8, 37, 25, 88, 19, 86, 97, 142, 74, 60, 109, 41, 98, 53, 20, 139, 144, 45, 14, 143, 113, 129, 81, 83, 43, 94, 147, 40, 85, 35, 122, 124, 66, 89, 78, 84, 26, 131, 102, 140, 44, 137, 130, 29, 10, 148, 13, 108, 92, 36, 76, 70, 24], 2 => Int[27, 69, 131, 65, 78, 142, 40, 127, 122, 71, 18, 143, 62, 76, 129, 46, 31, 77, 96, 125, 119, 42, 85, 61, 15, 82, 84, 148, 5, 7, 14, 16, 56, 104, 57, 67, 128, 114, 130, 17, 149, 29, 32, 50, 1, 133, 19, 22, 101, 59, 2, 112, 88, 150, 9, 118, 113, 115, 63, 98, 108, 110, 74, 102, 21, 89, 49, 68, 90, 38, 45, 8, 37, 53, 11, 36, 116, 4, 94, 73, 146, 47, 58, 138, 72, 139, 111, 141, 70, 121, 132, 25, 64, 33, 123, 26, 97, 75, 12, 41, 30, 144, 24, 92, 103, 140, 80, 66, 6, 48, 105, 60, 28, 135, 44, 126, 136, 20, 3, 151, 91, 95, 55, 145, 106, 54, 100, 51, 10, 13, 86, 137, 99, 109, 52, 23, 134, 34, 81, 35, 39, 83, 79, 120, 107, 87, 93, 124, 117, 43, 147])

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

const _A = [0.5 0 0 0.1; 1.2 1 0 0.5; 1 0 -0.5 0.5; 1 -0.5 0 0.5; 1.2 0 1 0.5]
function _orig(x::Vector{Float64})
    s = 0.0
    for r in 1:5; s += _A[r, 1] * exp(-((x[1] - _A[r, 2])^2 + (x[2] - _A[r, 3])^2) / _A[r, 4]^2); end
    return -s
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
        f_star = -1.296954045953779,
        f_star_justification = "L",
        n_minimizers_documented = 1,
        description = "DFO Benchmark Suite v2.1.0 core (Julia): K adaptation of multi_gaussian (F = f o T).",
        x0_idx = Float64[85.0, 88.0],
        x0_alt_idx = Vector{Float64}[Float64[12.0, 93.0], Float64[133.0, 97.0], Float64[59.0, 85.0], Float64[95.0, 94.0], Float64[21.0, 98.0], Float64[114.0, 88.0], Float64[41.0, 90.0], Float64[76.0, 99.0], Float64[2.0, 89.0]],
        suite_release = "2.1.0",
        base_function = "multi_gaussian",
        m = 151,
        pattern = "K",
        strategy = "baseline",
        scale_factors_C = Float64[],
        lb_orig = Float64[-1.3062646484375, -1.3355687649240462],
        ub_orig = Float64[2.0, 2.0],
        lb_C_work = Float64[-1.3062646484375, -1.3355687649240462],
        ub_C_work = Float64[2.0, 2.0],
        ProblemData = _problem_data(),
        func_idx = _objective,
        func_F = xr -> _objective(_cell2idx(xr)),
        idx2cell = _idx2cell,
        cell2idx = _cell2idx,
        decode = xr -> _decode(_cell2idx(xr)))
end

end # module
