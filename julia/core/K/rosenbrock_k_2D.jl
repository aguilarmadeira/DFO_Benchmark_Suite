# rosenbrock_k_2D  DFO Benchmark Suite v2.1 core, Julia version (self-contained; Base only).
#ROSENBROCK_K_2D  DFO Benchmark Suite v2.1.0 core: rosenbrock, dimension 2, pattern K, strategy baseline.
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
# Julia version of the MATLAB wrapper rosenbrock_k_2D.m: same instance, same values (level indices 1..m as in MATLAB).
# Usage:  include("rosenbrock_k_2D.jl");  PB = rosenbrock_k_2D.problem();  f = PB.func_idx(PB.x0_idx)
# Copyright (C) 2026 J. F. A. Madeira. LGPL-3.0-or-later.
module rosenbrock_k_2D

const NAME = "rosenbrock_k_2D"
const N = 2
const TYPES = Char['K', 'K']
const CARD = Float64[151.0, 151.0]
const LB = Float64[-5.12, -5.12]
const UB = Float64[3.4252375, 3.6314412233907434]
const LBW = Float64[-5.12, -5.12]
const UBW = Float64[3.4252375, 3.6314412233907434]
const TB = Float64[0.9756027500000002, 1.00600885637352]
const TI = Float64[1.0, 1.0]
const CW = Float64[NaN, NaN]
const CX = Float64[NaN, NaN]
const M = 151
const PERM = Dict{Int,Vector{Int}}(1 => Int[83, 99, 111, 123, 88, 28, 92, 44, 77, 45, 97, 142, 68, 72, 55, 12, 54, 90, 113, 65, 4, 138, 115, 13, 96, 22, 76, 5, 9, 11, 19, 21, 82, 110, 47, 64, 53, 66, 51, 60, 125, 41, 84, 150, 144, 69, 35, 118, 120, 109, 23, 87, 148, 85, 37, 49, 124, 106, 27, 46, 121, 130, 149, 43, 2, 133, 145, 52, 101, 89, 137, 71, 128, 141, 78, 36, 56, 8, 117, 3, 129, 103, 75, 57, 14, 31, 105, 98, 25, 100, 70, 91, 42, 16, 40, 114, 95, 132, 30, 116, 59, 26, 112, 139, 63, 140, 122, 62, 7, 143, 24, 108, 93, 17, 29, 61, 39, 79, 73, 74, 94, 67, 146, 1, 151, 126, 86, 127, 38, 102, 147, 80, 20, 135, 48, 81, 33, 6, 32, 131, 134, 136, 119, 15, 104, 58, 18, 10, 107, 50, 34], 2 => Int[51, 116, 145, 124, 122, 144, 110, 95, 39, 121, 114, 117, 32, 52, 104, 30, 64, 85, 19, 82, 4, 7, 59, 136, 18, 22, 150, 5, 9, 12, 77, 129, 26, 102, 35, 141, 118, 87, 23, 93, 103, 41, 57, 56, 108, 62, 36, 81, 49, 8, 106, 113, 75, 128, 76, 74, 97, 94, 29, 125, 17, 89, 91, 45, 92, 96, 71, 54, 61, 27, 86, 1, 140, 50, 11, 131, 58, 10, 133, 6, 40, 90, 78, 60, 16, 99, 127, 119, 28, 137, 107, 112, 46, 147, 44, 21, 2, 130, 34, 66, 63, 111, 67, 115, 139, 24, 151, 68, 13, 80, 31, 3, 101, 25, 146, 123, 47, 88, 83, 84, 105, 79, 70, 126, 142, 14, 132, 15, 55, 120, 143, 100, 138, 72, 73, 109, 65, 42, 69, 33, 43, 48, 38, 148, 53, 37, 20, 98, 149, 134, 135])

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
    n = length(x)
    s = 0.0
    for i in 1:n-1; s += (x[i] - 1)^2 + 100 * (x[i]^2 - x[i+1])^2; end
    return s
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
        f_star = 0.0,
        f_star_justification = "E",
        n_minimizers_documented = 1,
        description = "DFO Benchmark Suite v2.1.0 core (Julia): K adaptation of rosenbrock (F = f o T).",
        x0_idx = Float64[81.0, 44.0],
        x0_alt_idx = Vector{Float64}[Float64[7.0, 48.0], Float64[142.0, 38.0], Float64[68.0, 40.0], Float64[104.0, 50.0], Float64[31.0, 39.0], Float64[124.0, 43.0], Float64[50.0, 46.0], Float64[86.0, 35.0], Float64[12.0, 44.0]],
        suite_release = "2.1.0",
        base_function = "rosenbrock",
        m = 151,
        pattern = "K",
        strategy = "baseline",
        scale_factors_C = Float64[],
        lb_orig = Float64[-5.12, -5.12],
        ub_orig = Float64[3.4252375, 3.6314412233907434],
        lb_C_work = Float64[-5.12, -5.12],
        ub_C_work = Float64[3.4252375, 3.6314412233907434],
        ProblemData = _problem_data(),
        func_idx = _objective,
        func_F = xr -> _objective(_cell2idx(xr)),
        idx2cell = _idx2cell,
        cell2idx = _cell2idx,
        decode = xr -> _decode(_cell2idx(xr)))
end

end # module
