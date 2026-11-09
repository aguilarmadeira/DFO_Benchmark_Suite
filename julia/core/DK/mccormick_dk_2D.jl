# mccormick_dk_2D  DFO Benchmark Suite v2.1 core, Julia version (self-contained; Base only).
#MCCORMICK_DK_2D  DFO Benchmark Suite v2.1.0 core: mccormick, dimension 2, pattern DK, strategy baseline.
#
# Copyright (C) 2026 J. F. A. Madeira
#
# Objective F(z) = f(T(z)), where f is the original function (embedded below, translated from the MATLAB original) and T acts on the finite
# coordinates only: uniform grid of 151 levels on the interval of the domain, piecewise-affine map with |T - Id| <=
# half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
# Global optimum of F on this domain: f* = -1.913222954981037 (L: literature value at a refined point); the optimal solution is listed in
# reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
# Variables: x1:D x2:K
#
# Julia version of the MATLAB wrapper mccormick_dk_2D.m: same instance, same values (level indices 1..m as in MATLAB).
# Usage:  include("mccormick_dk_2D.jl");  PB = mccormick_dk_2D.problem();  f = PB.func_idx(PB.x0_idx)
# Copyright (C) 2026 J. F. A. Madeira. LGPL-3.0-or-later.
module mccormick_dk_2D

const NAME = "mccormick_dk_2D"
const N = 2
const TYPES = Char['D', 'K']
const CARD = Float64[151.0, 151.0]
const LB = Float64[-1.5, -3.0]
const UB = Float64[3.1183013916015625, 2.0500198140527357]
const LBW = Float64[-1.5, -3.0]
const UBW = Float64[3.1183013916015625, 2.0500198140527357]
const TB = Float64[-0.545551045735677, -1.5523276533048824]
const TI = Float64[-0.5471975537167932, -1.547197552675593]
const CW = Float64[NaN, NaN]
const CX = Float64[NaN, NaN]
const M = 151
const PERM = Dict{Int,Vector{Int}}(2 => Int[48, 37, 89, 121, 20, 41, 15, 122, 64, 145, 126, 10, 27, 32, 4, 6, 80, 51, 43, 36, 112, 127, 77, 148, 55, 81, 50, 141, 78, 88, 1, 86, 79, 95, 12, 14, 68, 31, 102, 151, 99, 53, 119, 70, 24, 129, 42, 97, 71, 114, 94, 111, 109, 105, 150, 65, 52, 116, 17, 146, 115, 9, 13, 23, 44, 85, 101, 147, 117, 59, 120, 137, 113, 73, 143, 130, 83, 106, 49, 118, 39, 93, 72, 46, 82, 2, 30, 45, 75, 7, 125, 21, 103, 26, 47, 128, 84, 104, 110, 19, 134, 138, 98, 60, 63, 136, 29, 11, 58, 91, 35, 123, 57, 142, 135, 62, 28, 8, 108, 139, 140, 25, 132, 18, 3, 149, 100, 16, 74, 133, 87, 96, 54, 5, 61, 92, 107, 66, 67, 56, 124, 144, 131, 76, 90, 34, 22, 40, 69, 33, 38])

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

_orig(x::Vector{Float64}) = sin(x[1] + x[2]) + (x[1] - x[2])^2 - 3 / 2 * x[1] + 5 / 2 * x[2] + 1

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
        f_star = -1.9132229549810367,
        f_star_justification = "L",
        n_minimizers_documented = 1,
        description = "DFO Benchmark Suite v2.1.0 core (Julia): DK adaptation of mccormick (F = f o T).",
        x0_idx = Float64[57.0, 72.0],
        x0_alt_idx = Vector{Float64}[Float64[96.0, 76.0], Float64[20.0, 80.0], Float64[115.0, 68.0], Float64[39.0, 77.0], Float64[78.0, 82.0], Float64[2.0, 71.0], Float64[150.0, 73.0], Float64[76.0, 83.0], Float64[111.0, 72.0]],
        suite_release = "2.1.0",
        base_function = "mccormick",
        m = 151,
        pattern = "DK",
        strategy = "baseline",
        scale_factors_C = Float64[],
        lb_orig = Float64[-1.5, -3.0],
        ub_orig = Float64[3.1183013916015625, 2.0500198140527357],
        lb_C_work = Float64[-1.5, -3.0],
        ub_C_work = Float64[3.1183013916015625, 2.0500198140527357],
        ProblemData = _problem_data(),
        func_idx = _objective,
        func_F = xr -> _objective(_cell2idx(xr)),
        idx2cell = _idx2cell,
        cell2idx = _cell2idx,
        decode = xr -> _decode(_cell2idx(xr)))
end

end # module
