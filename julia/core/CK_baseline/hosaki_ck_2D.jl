# hosaki_ck_2D  DFO Benchmark Suite v2.1 core, Julia version (self-contained; Base only).
#HOSAKI_CK_2D  DFO Benchmark Suite v2.1.0 core: hosaki, dimension 2, pattern CK, strategy baseline.
#
# Copyright (C) 2026 J. F. A. Madeira
#
# Objective F(z) = f(T(z)), where f is the original function (embedded below, translated from the MATLAB original) and T acts on the finite
# coordinates only: uniform grid of 65 levels on the interval of the domain, piecewise-affine map with |T - Id| <=
# half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
# Global optimum of F on this domain: f* = -2.345811576101292 (L: literature value at a refined point); the optimal solution is listed in
# reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
# Variables: x1:C x2:K
#
# Julia version of the MATLAB wrapper hosaki_ck_2D.m: same instance, same values (level indices 1..m as in MATLAB).
# Usage:  include("hosaki_ck_2D.jl");  PB = hosaki_ck_2D.problem();  f = PB.func_idx(PB.x0_idx)
# Copyright (C) 2026 J. F. A. Madeira. LGPL-3.0-or-later.
module hosaki_ck_2D

const NAME = "hosaki_ck_2D"
const N = 2
const TYPES = Char['C', 'K']
const CARD = Float64[Inf, 65.0]
const LB = Float64[0.0, 1.027757963725042]
const UB = Float64[4.067205810546875, 6.0]
const LBW = Float64[0.0, 1.027757963725042]
const UBW = Float64[4.067205810546875, 6.0]
const TB = Float64[NaN, 2.037744627343393]
const TI = Float64[NaN, 2.0]
const CW = Float64[4.0, NaN]
const CX = Float64[4.0, NaN]
const M = 65
const PERM = Dict{Int,Vector{Int}}(2 => Int[11, 15, 43, 7, 60, 46, 56, 64, 42, 36, 45, 38, 50, 22, 61, 24, 34, 51, 41, 39, 17, 57, 1, 35, 18, 5, 53, 25, 20, 19, 9, 62, 49, 12, 6, 23, 4, 16, 8, 27, 33, 52, 47, 32, 31, 54, 3, 13, 26, 58, 59, 37, 63, 55, 29, 28, 14, 10, 30, 48, 65, 44, 21, 2, 40])

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
        cardinality = Float64[Inf, 65.0],
        f_star = -2.345811576101292,
        f_star_justification = "L",
        n_minimizers_documented = 1,
        description = "DFO Benchmark Suite v2.1.0 core (Julia): CK adaptation of hosaki (F = f o T).",
        x0_idx = Float64[3.0523902982473374, 24.0],
        x0_alt_idx = Vector{Float64}[Float64[1.0505624383687975, 26.0], Float64[2.0514763683080677, 28.0], Float64[0.04964850842952729, 23.0], Float64[4.053304228186608, 26.0], Float64[1.9799825161695481, 28.0], Float64[3.0285590142011647, 24.0], Float64[0.9790685862302781, 25.0], Float64[3.5528472632169725, 29.0], Float64[1.4676099091768267, 24.0]],
        suite_release = "2.1.0",
        base_function = "hosaki",
        m = 65,
        pattern = "CK",
        strategy = "baseline",
        scale_factors_C = Float64[1.0],
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
