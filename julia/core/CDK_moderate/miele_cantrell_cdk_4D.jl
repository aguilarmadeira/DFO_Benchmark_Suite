# miele_cantrell_cdk_4D  DFO Benchmark Suite v2.1 core, Julia version (self-contained; Base only).
#MIELE_CANTRELL_CDK_4D  DFO Benchmark Suite v2.1.0 core: miele_cantrell, dimension 4, pattern CDK, strategy moderate.
#
# Copyright (C) 2026 J. F. A. Madeira
#
# Objective F(z) = f(T(z)), where f is the original function (embedded below, translated from the MATLAB original) and T acts on the finite
# coordinates only: uniform grid of 151 levels on the interval of the domain, piecewise-affine map with |T - Id| <=
# half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
# Global optimum of F on this domain: f* = 0 (E: exact by structure); the optimal solution is listed in
# reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
# Variables: x1:C x2:K x3:D x4:C
#
# Julia version of the MATLAB wrapper miele_cantrell_cdk_4D.m: same instance, same values (level indices 1..m as in MATLAB).
# Usage:  include("miele_cantrell_cdk_4D.jl");  PB = miele_cantrell_cdk_4D.problem();  f = PB.func_idx(PB.x0_idx)
# Copyright (C) 2026 J. F. A. Madeira. LGPL-3.0-or-later.
module miele_cantrell_cdk_4D

const NAME = "miele_cantrell_cdk_4D"
const N = 4
const TYPES = Char['C', 'K', 'D', 'C']
const CARD = Float64[Inf, 151.0, 151.0, Inf]
const LB = Float64[-10.0, -10.0, -10.0, -10.0]
const UB = Float64[7.2313232421875, 7.144510491286898, 6.769932799999999, 6.4983756768013325]
const LBW = Float64[-10.0, -10.0, -10.0, -1000.0]
const UBW = Float64[7.2313232421875, 7.144510491286898, 6.769932799999999, 649.8375676801331]
const TB = Float64[NaN, 0.9724867144236153, 0.9563560960000004, NaN]
const TI = Float64[NaN, 1.0, 1.0, NaN]
const CW = Float64[0.0, NaN, NaN, 100.0]
const CX = Float64[0.0, NaN, NaN, 1.0]
const M = 151
const PERM = Dict{Int,Vector{Int}}(2 => Int[77, 122, 104, 101, 150, 18, 6, 128, 114, 56, 57, 68, 102, 78, 82, 129, 36, 48, 71, 19, 107, 52, 111, 8, 79, 123, 141, 26, 127, 3, 96, 110, 15, 73, 66, 131, 86, 121, 140, 92, 16, 49, 97, 31, 119, 4, 11, 14, 70, 139, 44, 81, 134, 142, 105, 9, 76, 94, 41, 58, 25, 137, 133, 80, 7, 23, 74, 53, 126, 22, 28, 63, 116, 148, 90, 17, 34, 62, 115, 51, 145, 24, 35, 10, 47, 124, 84, 108, 42, 67, 65, 85, 32, 93, 20, 87, 132, 37, 72, 1, 106, 50, 83, 138, 54, 89, 27, 12, 60, 99, 45, 147, 151, 98, 88, 95, 64, 55, 38, 5, 120, 112, 69, 118, 130, 59, 135, 39, 117, 113, 43, 61, 40, 29, 100, 143, 125, 136, 149, 91, 21, 33, 46, 103, 2, 13, 109, 144, 30, 146, 75])

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

_orig(x::Vector{Float64}) = (exp(x[1]) - x[2])^4 + 100 * (x[2] - x[3])^6 + tan(x[3] - x[4])^4 + x[1]^8

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
        cardinality = Float64[Inf, 151.0, 151.0, Inf],
        f_star = 0.0,
        f_star_justification = "E",
        n_minimizers_documented = 1,
        description = "DFO Benchmark Suite v2.1.0 core (Julia): CDK adaptation of miele_cantrell (F = f o T).",
        x0_idx = Float64[-8.914628565311432, 38.0, 69.0, -290.8653186814256],
        x0_alt_idx = Vector{Float64}[Float64[5.220441281795502, 43.0, 29.0, -218.71498773331484], Float64[-3.462530195713043, 47.0, 99.0, -146.56465678520385], Float64[0.9799203276634216, 35.0, 60.0, -74.4143258370932], Float64[-7.804015934467316, 44.0, 128.0, -233.14505392293688], Float64[3.1001808047294617, 48.0, 89.0, -160.994722974826], Float64[-5.582790672779083, 38.0, 15.0, -121.8274004601376], Float64[-1.1403401494026184, 40.0, 119.0, -280.55812854598116], Float64[-9.974758803844452, 49.0, 45.0, -208.4077975978704], Float64[6.987325012683868, 39.0, 149.0, -136.25746664975975]],
        suite_release = "2.1.0",
        base_function = "miele_cantrell",
        m = 151,
        pattern = "CDK",
        strategy = "moderate",
        scale_factors_C = Float64[1.0, 100.0],
        lb_orig = Float64[-10.0, -10.0, -10.0, -10.0],
        ub_orig = Float64[7.2313232421875, 7.144510491286898, 6.769932799999999, 6.4983756768013325],
        lb_C_work = Float64[-10.0, -10.0, -10.0, -1000.0],
        ub_C_work = Float64[7.2313232421875, 7.144510491286898, 6.769932799999999, 649.8375676801331],
        ProblemData = _problem_data(),
        func_idx = _objective,
        func_F = xr -> _objective(_cell2idx(xr)),
        idx2cell = _idx2cell,
        cell2idx = _cell2idx,
        decode = xr -> _decode(_cell2idx(xr)))
end

end # module
