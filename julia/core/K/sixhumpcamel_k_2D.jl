# sixhumpcamel_k_2D  DFO Benchmark Suite v2.1 core, Julia version (self-contained; Base only).
#SIXHUMPCAMEL_K_2D  DFO Benchmark Suite v2.1.0 core: sixhumpcamel, dimension 2, pattern K, strategy baseline.
#
# Copyright (C) 2026 J. F. A. Madeira
#
# Objective F(z) = f(T(z)), where f is the original function (embedded below, translated from the MATLAB original) and T acts on the finite
# coordinates only: uniform grid of 151 levels on the interval of the domain, piecewise-affine map with |T - Id| <=
# half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
# Global optimum of F on this domain: f* = -1.031628453488552 (A: analytic reduction, 1-D root to machine precision); the optimal solution is listed in
# reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
# Variables: x1:K x2:K
#
# Julia version of the MATLAB wrapper sixhumpcamel_k_2D.m: same instance, same values (level indices 1..m as in MATLAB).
# Usage:  include("sixhumpcamel_k_2D.jl");  PB = sixhumpcamel_k_2D.problem();  f = PB.func_idx(PB.x0_idx)
# Copyright (C) 2026 J. F. A. Madeira. LGPL-3.0-or-later.
module sixhumpcamel_k_2D

const NAME = "sixhumpcamel_k_2D"
const N = 2
const TYPES = Char['K', 'K']
const CARD = Float64[151.0, 151.0]
const LB = Float64[-2.03322509765625, -1.3236620433876949]
const UB = Float64[3.0, 2.0]
const LBW = Float64[-2.03322509765625, -1.3236620433876949]
const UBW = Float64[3.0, 2.0]
const TB = Float64[0.08072944335937482, -0.7032451286219918]
const TI = Float64[0.089842, -0.712656]
const CW = Float64[NaN, NaN]
const CX = Float64[NaN, NaN]
const M = 151
const PERM = Dict{Int,Vector{Int}}(1 => Int[145, 150, 39, 49, 14, 140, 75, 129, 124, 34, 26, 109, 70, 144, 77, 6, 55, 95, 80, 57, 19, 79, 42, 30, 108, 17, 127, 48, 113, 7, 103, 18, 23, 40, 71, 68, 62, 56, 126, 119, 122, 59, 72, 13, 52, 36, 83, 147, 137, 128, 96, 138, 32, 37, 8, 22, 118, 51, 81, 67, 76, 66, 15, 112, 110, 9, 24, 46, 106, 25, 115, 123, 41, 105, 148, 58, 134, 38, 74, 92, 107, 31, 45, 141, 98, 111, 125, 151, 35, 114, 20, 44, 1, 116, 10, 94, 121, 82, 139, 89, 104, 78, 16, 100, 61, 131, 3, 12, 73, 69, 120, 28, 43, 101, 33, 102, 93, 88, 133, 29, 86, 87, 5, 135, 65, 142, 84, 21, 63, 97, 54, 47, 130, 4, 99, 53, 11, 50, 136, 2, 143, 60, 90, 146, 132, 64, 117, 27, 149, 91, 85], 2 => Int[103, 143, 53, 49, 141, 87, 27, 42, 142, 15, 71, 20, 58, 28, 64, 73, 55, 128, 149, 61, 122, 130, 3, 72, 12, 31, 65, 80, 59, 8, 126, 133, 91, 125, 34, 1, 62, 45, 132, 30, 54, 100, 77, 95, 113, 116, 46, 147, 105, 22, 39, 51, 136, 84, 119, 23, 98, 52, 5, 70, 109, 10, 123, 75, 69, 11, 25, 102, 56, 26, 135, 7, 137, 63, 134, 106, 94, 120, 110, 47, 6, 33, 92, 89, 85, 32, 68, 88, 38, 121, 151, 48, 117, 129, 14, 104, 90, 138, 57, 93, 144, 83, 21, 2, 66, 60, 9, 18, 79, 76, 36, 139, 50, 108, 41, 111, 101, 97, 107, 40, 96, 127, 17, 112, 78, 131, 99, 37, 145, 148, 74, 67, 19, 29, 124, 82, 146, 86, 35, 44, 140, 114, 16, 150, 81, 43, 118, 115, 13, 4, 24])

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

_orig(x::Vector{Float64}) = (4 - 2.1 * x[1]^2 + x[1]^4 / 3) * x[1]^2 + x[1] * x[2] + (4 * x[2]^2 - 4) * x[2]^2

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
        f_star = -1.0316284534885518,
        f_star_justification = "A",
        n_minimizers_documented = 2,
        description = "DFO Benchmark Suite v2.1.0 core (Julia): K adaptation of sixhumpcamel (F = f o T).",
        x0_idx = Float64[71.0, 52.0],
        x0_alt_idx = Vector{Float64}[Float64[110.0, 61.0], Float64[34.0, 65.0], Float64[129.0, 53.0], Float64[53.0, 58.0], Float64[92.0, 67.0], Float64[16.0, 56.0], Float64[139.0, 59.0], Float64[62.0, 63.0], Float64[101.0, 52.0]],
        suite_release = "2.1.0",
        base_function = "sixhumpcamel",
        m = 151,
        pattern = "K",
        strategy = "baseline",
        scale_factors_C = Float64[],
        lb_orig = Float64[-2.03322509765625, -1.3236620433876949],
        ub_orig = Float64[3.0, 2.0],
        lb_C_work = Float64[-2.03322509765625, -1.3236620433876949],
        ub_C_work = Float64[3.0, 2.0],
        ProblemData = _problem_data(),
        func_idx = _objective,
        func_F = xr -> _objective(_cell2idx(xr)),
        idx2cell = _idx2cell,
        cell2idx = _cell2idx,
        decode = xr -> _decode(_cell2idx(xr)))
end

end # module
