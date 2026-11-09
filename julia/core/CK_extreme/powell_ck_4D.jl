# powell_ck_4D  DFO Benchmark Suite v2.1 core, Julia version (self-contained; Base only).
#POWELL_CK_4D  DFO Benchmark Suite v2.1.0 core: powell, dimension 4, pattern CK, strategy extreme.
#
# Copyright (C) 2026 J. F. A. Madeira
#
# Objective F(z) = f(T(z)), where f is the original function (embedded below, translated from the MATLAB original) and T acts on the finite
# coordinates only: uniform grid of 151 levels on the interval of the domain, piecewise-affine map with |T - Id| <=
# half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
# Global optimum of F on this domain: f* = 0 (E: exact by structure); the optimal solution is listed in
# reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
# Variables: x1:C x2:K x3:C x4:K
#
# Julia version of the MATLAB wrapper powell_ck_4D.m: same instance, same values (level indices 1..m as in MATLAB).
# Usage:  include("powell_ck_4D.jl");  PB = powell_ck_4D.problem();  f = PB.func_idx(PB.x0_idx)
# Copyright (C) 2026 J. F. A. Madeira. LGPL-3.0-or-later.
module powell_ck_4D

const NAME = "powell_ck_4D"
const N = 4
const TYPES = Char['C', 'K', 'C', 'K']
const CARD = Float64[Inf, 151.0, Inf, 151.0]
const LB = Float64[-6.0063232421875, -6.781547528323934, -10.0, -10.0]
const UB = Float64[10.0, 10.0, 6.1539328, 6.469804248229904]
const LBW = Float64[-6.0063232421874995, -6.781547528323934, -1000000000.0, -10.0]
const UBW = Float64[10.0, 10.0, 615393280.0, 6.469804248229904]
const TB = Float64[NaN, 0.04294846652779949, NaN, -0.008318756073856903]
const TI = Float64[NaN, 0.0, NaN, 0.0]
const CW = Float64[8.881784197001252e-16, NaN, 0.0, NaN]
const CX = Float64[0.0, NaN, 0.0, NaN]
const M = 151
const PERM = Dict{Int,Vector{Int}}(2 => Int[86, 57, 13, 44, 125, 103, 104, 143, 63, 119, 66, 118, 111, 107, 126, 98, 112, 151, 93, 144, 94, 56, 69, 145, 110, 73, 135, 35, 89, 122, 115, 102, 16, 137, 25, 139, 36, 149, 84, 127, 150, 88, 92, 38, 1, 87, 20, 130, 113, 121, 59, 15, 45, 49, 17, 29, 116, 55, 90, 108, 53, 71, 7, 81, 40, 133, 132, 23, 30, 65, 2, 39, 95, 14, 147, 76, 18, 48, 79, 43, 105, 136, 32, 142, 52, 67, 117, 28, 61, 5, 140, 31, 74, 96, 70, 134, 26, 8, 54, 85, 80, 146, 131, 51, 77, 128, 75, 68, 11, 72, 22, 101, 83, 3, 10, 124, 129, 42, 33, 27, 41, 114, 97, 37, 62, 24, 100, 120, 12, 60, 82, 4, 58, 47, 34, 19, 50, 109, 106, 138, 6, 91, 123, 64, 148, 46, 78, 141, 9, 21, 99], 4 => Int[70, 112, 50, 9, 72, 129, 66, 148, 104, 36, 141, 136, 49, 138, 140, 59, 133, 108, 123, 33, 20, 86, 30, 100, 44, 58, 149, 111, 151, 76, 11, 13, 74, 53, 118, 28, 83, 150, 101, 16, 77, 143, 32, 103, 127, 14, 22, 25, 69, 65, 7, 135, 120, 126, 19, 31, 145, 110, 125, 115, 134, 71, 10, 21, 43, 121, 8, 27, 105, 68, 79, 139, 93, 107, 95, 4, 23, 119, 5, 48, 1, 109, 37, 98, 142, 73, 91, 113, 67, 12, 137, 39, 81, 116, 78, 55, 34, 17, 124, 94, 90, 56, 92, 61, 88, 38, 144, 147, 128, 85, 35, 146, 130, 18, 26, 54, 29, 60, 52, 47, 62, 15, 131, 63, 89, 51, 2, 24, 46, 97, 122, 45, 102, 96, 84, 75, 114, 41, 42, 82, 106, 57, 117, 87, 40, 3, 80, 6, 64, 99, 132])

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

_orig(x::Vector{Float64}) = (x[1] + 10 * x[2])^2 + 5 * (x[3] - x[4])^2 + (x[2] - 2 * x[3])^4 + 10 * (x[1] - x[4])^4

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
        f_star = 0.0,
        f_star_justification = "E",
        n_minimizers_documented = 1,
        description = "DFO Benchmark Suite v2.1.0 core (Julia): CK adaptation of powell (F = f o T).",
        x0_idx = Float64[9.00742038488388, 77.0, 442740046.23360014, 69.0],
        x0_alt_idx = Vector{Float64}[Float64[1.1293081641197205, 82.0, -293879289.4463999, 75.0], Float64[5.068364274501801, 71.0, -782374217.3184, 82.0], Float64[-2.8097479462623594, 74.0, 31783995.8016001, 67.0], Float64[7.03789232969284, 83.0, -469117152.4608, 74.0], Float64[-0.8402198910713192, 72.0, 345041060.6592002, 80.0], Float64[3.098836219310761, 77.0, -143453867.2127999, 84.0], Float64[-4.7792760014533995, 79.0, -880073202.8928, 69.0], Float64[8.022656357288358, 68.0, 182209418.03520012, 76.0], Float64[0.14454413652420062, 78.0, -554409917.6448001, 83.0]],
        suite_release = "2.1.0",
        base_function = "powell",
        m = 151,
        pattern = "CK",
        strategy = "extreme",
        scale_factors_C = Float64[1.0, 100000000.0],
        lb_orig = Float64[-6.0063232421875, -6.781547528323934, -10.0, -10.0],
        ub_orig = Float64[10.0, 10.0, 6.1539328, 6.469804248229904],
        lb_C_work = Float64[-6.0063232421874995, -6.781547528323934, -1000000000.0, -10.0],
        ub_C_work = Float64[10.0, 10.0, 615393280.0, 6.469804248229904],
        ProblemData = _problem_data(),
        func_idx = _objective,
        func_F = xr -> _objective(_cell2idx(xr)),
        idx2cell = _idx2cell,
        cell2idx = _cell2idx,
        decode = xr -> _decode(_cell2idx(xr)))
end

end # module
