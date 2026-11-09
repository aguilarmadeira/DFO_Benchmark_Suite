# dekkers_aarts_k_2D  DFO Benchmark Suite v2.1 core, Julia version (self-contained; Base only).
#DEKKERS_AARTS_K_2D  DFO Benchmark Suite v2.1.0 core: dekkers_aarts, dimension 2, pattern K, strategy baseline.
#
# Copyright (C) 2026 J. F. A. Madeira
#
# Objective F(z) = f(T(z)), where f is the original function (embedded below, translated from the MATLAB original) and T acts on the finite
# coordinates only: uniform grid of 151 levels on the interval of the domain, piecewise-affine map with |T - Id| <=
# half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
# Global optimum of F on this domain: f* = -24776.51834231769 (A: analytic reduction, 1-D root to machine precision); the optimal solution is listed in
# reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
# Variables: x1:K x2:K
#
# Julia version of the MATLAB wrapper dekkers_aarts_k_2D.m: same instance, same values (level indices 1..m as in MATLAB).
# Usage:  include("dekkers_aarts_k_2D.jl");  PB = dekkers_aarts_k_2D.problem();  f = PB.func_idx(PB.x0_idx)
# Copyright (C) 2026 J. F. A. Madeira. LGPL-3.0-or-later.
module dekkers_aarts_k_2D

const NAME = "dekkers_aarts_k_2D"
const N = 2
const TYPES = Char['K', 'K']
const CARD = Float64[151.0, 151.0]
const LB = Float64[-20.0, -12.145811106030584]
const UB = Float64[14.025146484375, 20.0]
const LBW = Float64[-20.0, -12.145811106030584]
const UBW = Float64[14.025146484375, 20.0]
const TB = Float64[-0.03858072916666799, 14.856670223035108]
const TI = Float64[0.0, 14.94511215764898]
const CW = Float64[NaN, NaN]
const CX = Float64[NaN, NaN]
const M = 151
const PERM = Dict{Int,Vector{Int}}(1 => Int[104, 16, 90, 131, 76, 136, 122, 139, 87, 132, 60, 74, 130, 66, 103, 91, 70, 46, 19, 1, 55, 44, 33, 134, 41, 138, 101, 22, 51, 92, 2, 135, 64, 96, 58, 98, 18, 54, 28, 49, 109, 40, 89, 118, 143, 59, 95, 25, 31, 105, 48, 78, 144, 26, 73, 69, 84, 11, 117, 30, 50, 57, 62, 150, 88, 34, 149, 147, 99, 15, 83, 45, 67, 110, 145, 151, 61, 100, 24, 53, 77, 115, 20, 63, 106, 75, 125, 21, 43, 137, 102, 97, 56, 7, 23, 72, 128, 82, 81, 3, 68, 13, 36, 86, 141, 6, 39, 107, 142, 5, 47, 148, 140, 120, 121, 71, 146, 4, 79, 80, 12, 114, 108, 112, 129, 29, 85, 27, 10, 94, 123, 32, 37, 8, 52, 116, 38, 14, 124, 93, 133, 35, 9, 119, 126, 17, 113, 65, 42, 127, 111], 2 => Int[94, 104, 82, 92, 65, 55, 64, 115, 107, 21, 46, 67, 74, 130, 29, 105, 96, 137, 95, 2, 40, 116, 146, 34, 84, 15, 37, 71, 111, 39, 61, 3, 98, 85, 68, 147, 19, 48, 54, 60, 36, 120, 57, 5, 10, 62, 88, 81, 32, 38, 138, 122, 25, 28, 128, 52, 124, 106, 18, 145, 51, 93, 101, 113, 149, 35, 44, 144, 143, 102, 22, 119, 69, 7, 97, 50, 63, 108, 117, 56, 79, 151, 123, 66, 136, 135, 118, 24, 47, 103, 16, 30, 133, 131, 129, 76, 31, 150, 86, 148, 73, 17, 109, 91, 126, 110, 45, 6, 41, 12, 132, 72, 127, 87, 20, 80, 26, 13, 142, 90, 125, 4, 121, 1, 139, 141, 100, 43, 27, 112, 11, 53, 59, 33, 78, 9, 70, 49, 23, 134, 42, 89, 77, 58, 75, 140, 114, 8, 83, 99, 14])

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

_orig(x::Vector{Float64}) = 10.0^5 * x[1]^2 + x[2]^2 - (x[1]^2 + x[2]^2)^2 + 1 / 10.0^5 * (x[1]^2 + x[2]^2)^4

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
        f_star = -24776.51834231769,
        f_star_justification = "A",
        n_minimizers_documented = 1,
        description = "DFO Benchmark Suite v2.1.0 core (Julia): K adaptation of dekkers_aarts (F = f o T).",
        x0_idx = Float64[48.0, 109.0],
        x0_alt_idx = Vector{Float64}[Float64[83.0, 113.0], Float64[9.0, 103.0], Float64[145.0, 112.0], Float64[71.0, 114.0], Float64[107.0, 104.0], Float64[33.0, 108.0], Float64[126.0, 118.0], Float64[52.0, 105.0], Float64[88.0, 109.0]],
        suite_release = "2.1.0",
        base_function = "dekkers_aarts",
        m = 151,
        pattern = "K",
        strategy = "baseline",
        scale_factors_C = Float64[],
        lb_orig = Float64[-20.0, -12.145811106030584],
        ub_orig = Float64[14.025146484375, 20.0],
        lb_C_work = Float64[-20.0, -12.145811106030584],
        ub_C_work = Float64[14.025146484375, 20.0],
        ProblemData = _problem_data(),
        func_idx = _objective,
        func_F = xr -> _objective(_cell2idx(xr)),
        idx2cell = _idx2cell,
        cell2idx = _cell2idx,
        decode = xr -> _decode(_cell2idx(xr)))
end

end # module
