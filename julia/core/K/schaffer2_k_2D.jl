# schaffer2_k_2D  DFO Benchmark Suite v2.1 core, Julia version (self-contained; Base only).
#SCHAFFER2_K_2D  DFO Benchmark Suite v2.1.0 core: schaffer2, dimension 2, pattern K, strategy baseline.
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
# Julia version of the MATLAB wrapper schaffer2_k_2D.m: same instance, same values (level indices 1..m as in MATLAB).
# Usage:  include("schaffer2_k_2D.jl");  PB = schaffer2_k_2D.problem();  f = PB.func_idx(PB.x0_idx)
# Copyright (C) 2026 J. F. A. Madeira. LGPL-3.0-or-later.
module schaffer2_k_2D

const NAME = "schaffer2_k_2D"
const N = 2
const TYPES = Char['K', 'K']
const CARD = Float64[151.0, 151.0]
const LB = Float64[-100.0, -64.70436417212824]
const UB = Float64[65.149169921875, 100.0]
const LBW = Float64[-100.0, -64.70436417212824]
const UBW = Float64[65.149169921875, 100.0]
const TB = Float64[0.19049641927082916, 0.07935240224219342]
const TI = Float64[0.0, 0.0]
const CW = Float64[NaN, NaN]
const CX = Float64[NaN, NaN]
const M = 151
const PERM = Dict{Int,Vector{Int}}(1 => Int[55, 107, 66, 9, 88, 137, 110, 3, 95, 61, 134, 144, 41, 102, 83, 123, 111, 138, 69, 133, 131, 52, 149, 99, 31, 80, 62, 114, 121, 73, 38, 39, 146, 113, 21, 48, 103, 65, 51, 77, 87, 26, 29, 56, 93, 91, 20, 112, 33, 42, 75, 18, 127, 118, 72, 36, 74, 6, 15, 143, 104, 16, 22, 92, 58, 78, 100, 128, 119, 19, 27, 67, 85, 145, 8, 37, 94, 11, 44, 117, 105, 148, 2, 122, 28, 50, 115, 14, 49, 1, 116, 34, 136, 57, 82, 108, 47, 35, 130, 24, 25, 139, 32, 129, 46, 13, 64, 63, 132, 150, 43, 23, 17, 140, 84, 12, 109, 54, 53, 70, 101, 120, 151, 79, 5, 106, 86, 142, 125, 90, 59, 10, 135, 96, 126, 124, 141, 7, 97, 147, 68, 45, 98, 76, 81, 40, 60, 71, 30, 89, 4], 2 => Int[89, 65, 148, 76, 84, 77, 151, 119, 27, 93, 85, 22, 73, 41, 3, 90, 55, 32, 75, 44, 133, 34, 127, 2, 6, 11, 91, 141, 138, 111, 56, 144, 98, 62, 57, 82, 49, 58, 142, 131, 116, 96, 71, 103, 45, 66, 107, 113, 33, 43, 9, 106, 51, 117, 25, 115, 114, 94, 88, 36, 7, 146, 23, 35, 59, 78, 26, 47, 54, 21, 29, 130, 101, 81, 104, 39, 60, 129, 46, 12, 48, 74, 5, 149, 145, 53, 87, 109, 97, 4, 123, 38, 83, 128, 86, 67, 105, 40, 110, 28, 121, 140, 37, 122, 52, 125, 70, 69, 136, 95, 50, 132, 143, 68, 124, 20, 1, 64, 63, 80, 112, 10, 139, 92, 19, 134, 137, 147, 13, 108, 79, 30, 15, 120, 150, 17, 100, 42, 135, 102, 118, 99, 24, 8, 31, 18, 61, 14, 72, 16, 126])

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

_orig(x::Vector{Float64}) = (x[1]^2 + x[2]^2)^0.25 * (sin(50 * (x[1]^2 + x[2]^2)^0.1)^2 + 1)

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
        description = "DFO Benchmark Suite v2.1.0 core (Julia): K adaptation of schaffer2 (F = f o T).",
        x0_idx = Float64[100.0, 66.0],
        x0_alt_idx = Vector{Float64}[Float64[26.0, 56.0], Float64[119.0, 60.0], Float64[45.0, 63.0], Float64[81.0, 52.0], Float64[7.0, 61.0], Float64[143.0, 66.0], Float64[69.0, 53.0], Float64[104.0, 57.0], Float64[31.0, 67.0]],
        suite_release = "2.1.0",
        base_function = "schaffer2",
        m = 151,
        pattern = "K",
        strategy = "baseline",
        scale_factors_C = Float64[],
        lb_orig = Float64[-100.0, -64.70436417212824],
        ub_orig = Float64[65.149169921875, 100.0],
        lb_C_work = Float64[-100.0, -64.70436417212824],
        ub_C_work = Float64[65.149169921875, 100.0],
        ProblemData = _problem_data(),
        func_idx = _objective,
        func_F = xr -> _objective(_cell2idx(xr)),
        idx2cell = _idx2cell,
        cell2idx = _cell2idx,
        decode = xr -> _decode(_cell2idx(xr)))
end

end # module
