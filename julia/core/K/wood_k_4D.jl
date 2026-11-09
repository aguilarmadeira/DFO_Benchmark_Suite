# wood_k_4D  DFO Benchmark Suite v2.1 core, Julia version (self-contained; Base only).
#WOOD_K_4D  DFO Benchmark Suite v2.1.0 core: wood, dimension 4, pattern K, strategy baseline.
#
# Copyright (C) 2026 J. F. A. Madeira
#
# Objective F(z) = f(T(z)), where f is the original function (embedded below, translated from the MATLAB original) and T acts on the finite
# coordinates only: uniform grid of 25 levels on the interval of the domain, piecewise-affine map with |T - Id| <=
# half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
# Global optimum of F on this domain: f* = 0 (E: exact by structure); the optimal solution is listed in
# reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
# Variables: x1:K x2:K x3:K x4:K
#
# Julia version of the MATLAB wrapper wood_k_4D.m: same instance, same values (level indices 1..m as in MATLAB).
# Usage:  include("wood_k_4D.jl");  PB = wood_k_4D.problem();  f = PB.func_idx(PB.x0_idx)
# Copyright (C) 2026 J. F. A. Madeira. LGPL-3.0-or-later.
module wood_k_4D

const NAME = "wood_k_4D"
const N = 4
const TYPES = Char['K', 'K', 'K', 'K']
const CARD = Float64[25.0, 25.0, 25.0, 25.0]
const LB = Float64[-10.0, -10.0, -10.0, -10.0]
const UB = Float64[6.1211669921875, 6.566458365086623, 6.4227327999999995, 7.3514369012911285]
const LBW = Float64[-10.0, -10.0, -10.0, -10.0]
const UBW = Float64[6.1211669921875, 6.566458365086623, 6.4227327999999995, 7.3514369012911285]
const TB = Float64[0.7474446614583332, 1.0443055767244172, 0.9484885333333324, 0.8446480633069555]
const TI = Float64[1.0, 1.0, 1.0, 1.0]
const CW = Float64[NaN, NaN, NaN, NaN]
const CX = Float64[NaN, NaN, NaN, NaN]
const M = 25
const PERM = Dict{Int,Vector{Int}}(1 => Int[11, 17, 7, 19, 2, 13, 1, 3, 5, 23, 15, 18, 21, 24, 20, 9, 10, 25, 6, 12, 4, 22, 14, 16, 8], 2 => Int[20, 6, 12, 4, 16, 14, 21, 5, 7, 17, 10, 22, 8, 13, 3, 23, 1, 25, 18, 9, 24, 11, 15, 2, 19], 3 => Int[13, 20, 16, 1, 4, 17, 10, 6, 15, 14, 12, 23, 18, 21, 8, 5, 9, 24, 11, 7, 3, 25, 22, 19, 2], 4 => Int[18, 14, 22, 24, 10, 9, 16, 8, 17, 19, 4, 6, 25, 7, 13, 11, 1, 21, 5, 23, 3, 15, 2, 12, 20])

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

_orig(x::Vector{Float64}) = 100 * (x[1]^2 - x[2])^2 + (1 - x[1])^2 + 90 * (x[3]^2 - x[4])^2 + (1 - x[3])^2 + 10.1 * ((1 - x[2])^2 + (1 - x[4])^2) + 19.8 * (1 - x[2]) * (1 - x[4])

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
        cardinality = Float64[25.0, 25.0, 25.0, 25.0],
        f_star = 0.0,
        f_star_justification = "E",
        n_minimizers_documented = 1,
        description = "DFO Benchmark Suite v2.1.0 core (Julia): K adaptation of wood (F = f o T).",
        x0_idx = Float64[22.0, 10.0, 18.0, 24.0],
        x0_alt_idx = Vector{Float64}[Float64[11.0, 12.0, 6.0, 25.0], Float64[16.0, 10.0, 23.0, 22.0], Float64[5.0, 10.0, 11.0, 23.0], Float64[19.0, 11.0, 4.0, 24.0], Float64[8.0, 9.0, 15.0, 22.0], Float64[14.0, 11.0, 9.0, 23.0], Float64[2.0, 11.0, 20.0, 24.0], Float64[24.0, 9.0, 14.0, 25.0], Float64[12.0, 10.0, 2.0, 22.0]],
        suite_release = "2.1.0",
        base_function = "wood",
        m = 25,
        pattern = "K",
        strategy = "baseline",
        scale_factors_C = Float64[],
        lb_orig = Float64[-10.0, -10.0, -10.0, -10.0],
        ub_orig = Float64[6.1211669921875, 6.566458365086623, 6.4227327999999995, 7.3514369012911285],
        lb_C_work = Float64[-10.0, -10.0, -10.0, -10.0],
        ub_C_work = Float64[6.1211669921875, 6.566458365086623, 6.4227327999999995, 7.3514369012911285],
        ProblemData = _problem_data(),
        func_idx = _objective,
        func_F = xr -> _objective(_cell2idx(xr)),
        idx2cell = _idx2cell,
        cell2idx = _cell2idx,
        decode = xr -> _decode(_cell2idx(xr)))
end

end # module
