# shekel_410_ck_4D  DFO Benchmark Suite v2.1 core, Julia version (self-contained; Base only).
#SHEKEL_410_CK_4D  DFO Benchmark Suite v2.1.0 core: shekel_410, dimension 4, pattern CK, strategy moderate.
#
# Copyright (C) 2026 J. F. A. Madeira
#
# Objective F(z) = f(T(z)), where f is the original function (embedded below, translated from the MATLAB original) and T acts on the finite
# coordinates only: uniform grid of 151 levels on the interval of the domain, piecewise-affine map with |T - Id| <=
# half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
# Global optimum of F on this domain: f* = -10.53640981669204 (L: literature value at a refined point); the optimal solution is listed in
# reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
# Variables: x1:C x2:K x3:C x4:K
#
# Julia version of the MATLAB wrapper shekel_410_ck_4D.m: same instance, same values (level indices 1..m as in MATLAB).
# Usage:  include("shekel_410_ck_4D.jl");  PB = shekel_410_ck_4D.problem();  f = PB.func_idx(PB.x0_idx)
# Copyright (C) 2026 J. F. A. Madeira. LGPL-3.0-or-later.
module shekel_410_ck_4D

const NAME = "shekel_410_ck_4D"
const N = 4
const TYPES = Char['C', 'K', 'C', 'K']
const CARD = Float64[Inf, 151.0, Inf, 151.0]
const LB = Float64[1.69879150390625, 1.7686226693085405, 1.9790336000000002, 1.3936693044564765]
const UB = Float64[10.0, 10.0, 10.0, 10.0]
const LBW = Float64[1.69879150390625, 1.7686226693085405, 197.90336000000008, 1.3936693044564765]
const UBW = Float64[10.0, 10.0, 1000.0, 10.0]
const TB = Float64[NaN, 4.018532473030873, NaN, 3.9755685131195335]
const TI = Float64[NaN, 4.000592937456021, NaN, 3.999509802482246]
const CW = Float64[4.000746531113532, NaN, 399.9663397708698, NaN]
const CX = Float64[4.000746531113532, NaN, 3.999663397708697, NaN]
const M = 151
const PERM = Dict{Int,Vector{Int}}(2 => Int[147, 23, 126, 149, 3, 29, 119, 62, 123, 49, 116, 58, 26, 146, 89, 54, 102, 78, 82, 121, 41, 105, 11, 83, 109, 125, 77, 36, 122, 70, 127, 98, 32, 52, 151, 48, 80, 115, 90, 88, 59, 150, 14, 22, 30, 44, 6, 9, 128, 92, 45, 129, 37, 42, 12, 60, 34, 53, 66, 64, 81, 135, 106, 25, 79, 5, 18, 76, 148, 16, 24, 65, 140, 141, 10, 40, 113, 96, 51, 19, 112, 84, 15, 110, 124, 117, 21, 38, 131, 95, 50, 73, 27, 103, 33, 17, 143, 93, 56, 132, 130, 108, 7, 87, 134, 142, 67, 71, 20, 1, 68, 57, 137, 111, 31, 86, 72, 39, 43, 75, 118, 63, 85, 138, 114, 100, 99, 107, 101, 97, 133, 91, 120, 35, 94, 139, 136, 46, 61, 28, 47, 13, 4, 69, 8, 144, 104, 74, 2, 55, 145], 4 => Int[72, 142, 42, 144, 130, 118, 151, 64, 67, 38, 70, 126, 91, 59, 99, 97, 94, 50, 109, 104, 41, 53, 36, 68, 149, 133, 34, 81, 147, 116, 26, 61, 127, 146, 139, 48, 49, 113, 93, 90, 52, 114, 137, 23, 73, 148, 77, 128, 108, 30, 115, 136, 39, 44, 75, 122, 37, 56, 84, 112, 60, 10, 125, 131, 54, 8, 22, 43, 74, 19, 28, 69, 32, 85, 14, 88, 66, 21, 55, 24, 119, 129, 20, 7, 145, 78, 27, 45, 143, 120, 57, 102, 35, 123, 40, 106, 12, 4, 65, 100, 5, 80, 17, 98, 46, 25, 79, 134, 33, 15, 82, 124, 76, 13, 47, 103, 89, 58, 62, 95, 18, 87, 110, 71, 16, 2, 1, 11, 141, 3, 29, 150, 31, 83, 9, 138, 63, 117, 135, 111, 6, 121, 140, 92, 51, 101, 107, 105, 86, 132, 96])

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

const _A = [4 4 4 4 0.1; 1 1 1 1 0.2; 8 8 8 8 0.2; 6 6 6 6 0.4; 3 7 3 7 0.4;
            2 9 2 9 0.6; 5 5 3 3 0.3; 8 1 8 1 0.7; 6 2 6 2 0.5; 7 3.6 7 3.6 0.5]
function _orig(x::Vector{Float64})
    s = 0.0
    for r in 1:10
        e = 0.0
        for j in 1:4; e += (x[j] - _A[r, j])^2; end
        s += 1 / (e + _A[r, 5])
    end
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
        cardinality = Float64[Inf, 151.0, Inf, 151.0],
        f_star = -10.536409816692045,
        f_star_justification = "L",
        n_minimizers_documented = 1,
        description = "DFO Benchmark Suite v2.1.0 core (Julia): CK adaptation of shekel_410 (F = f o T).",
        x0_idx = Float64[6.631687372922897, 133.0, 978.4396423168002, 136.0],
        x0_alt_idx = Vector{Float64}[Float64[2.5459363162517548, 122.0, 612.6835744768, 142.0], Float64[9.25824162364006, 127.0, 216.12699566080008, 149.0], Float64[5.172490566968918, 129.0, 774.3862571008001, 135.0], Float64[7.16672620177269, 119.0, 371.6695760896001, 141.0], Float64[3.0809751451015472, 128.0, 929.9288375296001, 148.0], Float64[8.236803859472275, 132.0, 533.3722587136001, 130.0], Float64[4.151052802801132, 120.0, 321.6187457536001, 137.0], Float64[6.1209684908390045, 124.0, 695.0749413376002, 143.0], Float64[2.035217434167862, 133.0, 483.3214283776, 150.0]],
        suite_release = "2.1.0",
        base_function = "shekel_410",
        m = 151,
        pattern = "CK",
        strategy = "moderate",
        scale_factors_C = Float64[1.0, 100.0],
        lb_orig = Float64[1.69879150390625, 1.7686226693085405, 1.9790336000000002, 1.3936693044564765],
        ub_orig = Float64[10.0, 10.0, 10.0, 10.0],
        lb_C_work = Float64[1.69879150390625, 1.7686226693085405, 197.90336000000008, 1.3936693044564765],
        ub_C_work = Float64[10.0, 10.0, 1000.0, 10.0],
        ProblemData = _problem_data(),
        func_idx = _objective,
        func_F = xr -> _objective(_cell2idx(xr)),
        idx2cell = _idx2cell,
        cell2idx = _cell2idx,
        decode = xr -> _decode(_cell2idx(xr)))
end

end # module
