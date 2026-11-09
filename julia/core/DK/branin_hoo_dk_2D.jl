# branin_hoo_dk_2D  DFO Benchmark Suite v2.1 core, Julia version (self-contained; Base only).
#BRANIN_HOO_DK_2D  DFO Benchmark Suite v2.1.0 core: branin_hoo, dimension 2, pattern DK, strategy baseline.
#
# Copyright (C) 2026 J. F. A. Madeira
#
# Objective F(z) = f(T(z)), where f is the original function (embedded below, translated from the MATLAB original) and T acts on the finite
# coordinates only: uniform grid of 151 levels on the interval of the domain, piecewise-affine map with |T - Id| <=
# half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
# Global optimum of F on this domain: f* = 0.3978873577297382 (A: analytic reduction, 1-D root to machine precision); the optimal solution is listed in
# reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
# Variables: x1:D x2:K
#
# Julia version of the MATLAB wrapper branin_hoo_dk_2D.m: same instance, same values (level indices 1..m as in MATLAB).
# Usage:  include("branin_hoo_dk_2D.jl");  PB = branin_hoo_dk_2D.problem();  f = PB.func_idx(PB.x0_idx)
# Copyright (C) 2026 J. F. A. Madeira. LGPL-3.0-or-later.
module branin_hoo_dk_2D

const NAME = "branin_hoo_dk_2D"
const N = 2
const TYPES = Char['D', 'K']
const CARD = Float64[151.0, 151.0]
const LB = Float64[-5.0, 2.128654168571864]
const UB = Float64[7.283648681640625, 15.0]
const LBW = Float64[-5.0, 2.128654168571864]
const UBW = Float64[7.283648681640625, 15.0]
const TB = Float64[-3.116507202148438, 12.254112889295332]
const TI = Float64[-3.141592653589793, 12.275]
const CW = Float64[NaN, NaN]
const CX = Float64[NaN, NaN]
const M = 151
const PERM = Dict{Int,Vector{Int}}(2 => Int[110, 79, 44, 135, 119, 53, 17, 88, 144, 8, 145, 21, 36, 127, 80, 77, 128, 147, 96, 43, 87, 32, 14, 16, 38, 23, 61, 2, 51, 70, 71, 132, 15, 50, 20, 22, 107, 35, 136, 40, 57, 7, 151, 108, 140, 26, 62, 117, 134, 141, 98, 10, 66, 148, 46, 69, 3, 138, 19, 30, 68, 92, 78, 4, 18, 76, 37, 81, 83, 6, 105, 34, 54, 101, 12, 89, 41, 60, 1, 28, 49, 64, 104, 29, 116, 33, 52, 63, 84, 106, 146, 56, 82, 31, 121, 109, 42, 93, 85, 102, 123, 120, 133, 73, 11, 72, 103, 131, 122, 129, 86, 142, 111, 90, 9, 55, 113, 74, 24, 5, 139, 91, 143, 125, 59, 149, 99, 114, 47, 94, 58, 75, 27, 39, 48, 126, 25, 13, 45, 118, 95, 130, 65, 124, 137, 97, 100, 112, 67, 150, 115])

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

_orig(x::Vector{Float64}) = (x[2] - 1.275 * (x[1] / pi)^2 + 5 * x[1] / pi - 6)^2 + 10 * (1 - 1 / (8 * pi)) * cos(x[1]) + 10

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
        f_star = 0.39788735772973816,
        f_star_justification = "A",
        n_minimizers_documented = 2,
        description = "DFO Benchmark Suite v2.1.0 core (Julia): DK adaptation of branin_hoo (F = f o T).",
        x0_idx = Float64[109.0, 36.0],
        x0_alt_idx = Vector{Float64}[Float64[32.0, 41.0], Float64[126.0, 50.0], Float64[51.0, 40.0], Float64[90.0, 42.0], Float64[13.0, 47.0], Float64[136.0, 36.0], Float64[61.0, 45.0], Float64[99.0, 47.0], Float64[22.0, 37.0]],
        suite_release = "2.1.0",
        base_function = "branin_hoo",
        m = 151,
        pattern = "DK",
        strategy = "baseline",
        scale_factors_C = Float64[],
        lb_orig = Float64[-5.0, 2.128654168571864],
        ub_orig = Float64[7.283648681640625, 15.0],
        lb_C_work = Float64[-5.0, 2.128654168571864],
        ub_C_work = Float64[7.283648681640625, 15.0],
        ProblemData = _problem_data(),
        func_idx = _objective,
        func_F = xr -> _objective(_cell2idx(xr)),
        idx2cell = _idx2cell,
        cell2idx = _cell2idx,
        decode = xr -> _decode(_cell2idx(xr)))
end

end # module
