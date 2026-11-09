# branin_hoo_ck_2D  DFO Benchmark Suite v2.1 core, Julia version (self-contained; Base only).
#BRANIN_HOO_CK_2D  DFO Benchmark Suite v2.1.0 core: branin_hoo, dimension 2, pattern CK, strategy baseline.
#
# Copyright (C) 2026 J. F. A. Madeira
#
# Objective F(z) = f(T(z)), where f is the original function (embedded below, translated from the MATLAB original) and T acts on the finite
# coordinates only: uniform grid of 65 levels on the interval of the domain, piecewise-affine map with |T - Id| <=
# half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
# Global optimum of F on this domain: f* = 0.3978873577297382 (A: analytic reduction, 1-D root to machine precision); the optimal solution is listed in
# reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
# Variables: x1:C x2:K
#
# Julia version of the MATLAB wrapper branin_hoo_ck_2D.m: same instance, same values (level indices 1..m as in MATLAB).
# Usage:  include("branin_hoo_ck_2D.jl");  PB = branin_hoo_ck_2D.problem();  f = PB.func_idx(PB.x0_idx)
# Copyright (C) 2026 J. F. A. Madeira. LGPL-3.0-or-later.
module branin_hoo_ck_2D

const NAME = "branin_hoo_ck_2D"
const N = 2
const TYPES = Char['C', 'K']
const CARD = Float64[Inf, 65.0]
const LB = Float64[-5.0, 2.128654168571864]
const UB = Float64[7.283648681640625, 15.0]
const LBW = Float64[-5.0, 2.128654168571864]
const UBW = Float64[7.283648681640625, 15.0]
const TB = Float64[NaN, 12.184393099375097]
const TI = Float64[NaN, 12.275]
const CW = Float64[-3.141592653589793, NaN]
const CX = Float64[-3.141592653589793, NaN]
const M = 65
const PERM = Dict{Int,Vector{Int}}(2 => Int[28, 39, 9, 6, 8, 4, 33, 63, 24, 31, 17, 61, 54, 43, 58, 1, 59, 35, 34, 48, 21, 12, 25, 11, 41, 60, 36, 23, 3, 15, 30, 55, 7, 2, 13, 27, 45, 49, 19, 22, 44, 38, 16, 32, 20, 57, 51, 14, 53, 62, 10, 5, 18, 46, 37, 52, 26, 65, 56, 40, 42, 47, 29, 64, 50])

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
        cardinality = Float64[Inf, 65.0],
        f_star = 0.39788735772973816,
        f_star_justification = "A",
        n_minimizers_documented = 2,
        description = "DFO Benchmark Suite v2.1.0 core (Julia): CK adaptation of branin_hoo (F = f o T).",
        x0_idx = Float64[3.8168767392635345, 16.0],
        x0_alt_idx = Vector{Float64}[Float64[-2.498885986208916, 18.0], Float64[5.274360445141792, 22.0], Float64[-0.9154468983411794, 17.0], Float64[2.2514312773942944, 19.0], Float64[-4.0103505700826645, 20.0], Float64[6.066079989075661, 16.0], Float64[-0.1237273544073112, 20.0], Float64[3.0431508213281635, 21.0], Float64[-3.25461827814579, 16.0]],
        suite_release = "2.1.0",
        base_function = "branin_hoo",
        m = 65,
        pattern = "CK",
        strategy = "baseline",
        scale_factors_C = Float64[1.0],
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
