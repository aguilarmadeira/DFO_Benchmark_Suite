# exponencial_ck_4D  DFO Benchmark Suite v2.1 core, Julia version (self-contained; Base only).
#EXPONENCIAL_CK_4D  DFO Benchmark Suite v2.1.0 core: exponencial, dimension 4, pattern CK, strategy baseline.
#
# Copyright (C) 2026 J. F. A. Madeira
#
# Objective F(z) = f(T(z)), where f is the original function (embedded below, translated from the MATLAB original) and T acts on the finite
# coordinates only: uniform grid of 151 levels on the interval of the domain, piecewise-affine map with |T - Id| <=
# half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
# Global optimum of F on this domain: f* = -1 (E: exact by structure); the optimal solution is listed in
# reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
# Variables: x1:C x2:K x3:C x4:K
#
# Julia version of the MATLAB wrapper exponencial_ck_4D.m: same instance, same values (level indices 1..m as in MATLAB).
# Usage:  include("exponencial_ck_4D.jl");  PB = exponencial_ck_4D.problem();  f = PB.func_idx(PB.x0_idx)
# Copyright (C) 2026 J. F. A. Madeira. LGPL-3.0-or-later.
module exponencial_ck_4D

const NAME = "exponencial_ck_4D"
const N = 4
const TYPES = Char['C', 'K', 'C', 'K']
const CARD = Float64[Inf, 151.0, Inf, 151.0]
const LB = Float64[-0.68375732421875, -0.7109942590052329, -0.6165132799999999, -1.0]
const UB = Float64[1.0, 1.0, 1.0, 0.735551853394419]
const LBW = Float64[-0.68375732421875, -0.7109942590052329, -0.6165132799999999, -1.0]
const UBW = Float64[1.0, 1.0, 1.0, 0.735551853394419]
const TB = Float64[NaN, -0.0037832986164033633, NaN, -0.004950270720533112]
const TI = Float64[NaN, 0.0, NaN, 0.0]
const CW = Float64[0.0, NaN, 0.0, NaN]
const CX = Float64[0.0, NaN, 0.0, NaN]
const M = 151
const PERM = Dict{Int,Vector{Int}}(2 => Int[112, 146, 60, 82, 73, 111, 32, 59, 122, 68, 19, 48, 1, 2, 93, 37, 66, 8, 43, 38, 95, 16, 18, 75, 96, 102, 106, 10, 13, 15, 61, 126, 94, 55, 83, 35, 5, 12, 54, 104, 135, 26, 110, 76, 40, 142, 78, 63, 71, 23, 87, 98, 128, 77, 116, 101, 109, 58, 115, 20, 41, 49, 131, 80, 14, 27, 36, 91, 56, 86, 129, 145, 123, 141, 84, 51, 92, 9, 133, 62, 6, 120, 151, 124, 113, 103, 34, 46, 90, 127, 50, 65, 4, 52, 125, 137, 134, 118, 47, 70, 39, 139, 147, 140, 25, 108, 42, 149, 45, 28, 72, 130, 24, 17, 74, 30, 97, 107, 148, 85, 29, 143, 67, 79, 144, 22, 100, 89, 132, 99, 136, 114, 31, 81, 21, 150, 105, 53, 88, 64, 117, 3, 57, 119, 44, 121, 11, 69, 7, 138, 33], 4 => Int[113, 131, 140, 103, 58, 26, 80, 43, 122, 139, 97, 42, 102, 2, 46, 37, 96, 115, 116, 123, 82, 65, 101, 133, 64, 4, 118, 93, 83, 68, 60, 45, 27, 31, 34, 149, 6, 147, 70, 150, 132, 28, 66, 143, 63, 130, 138, 53, 106, 25, 100, 61, 146, 16, 117, 110, 72, 3, 9, 22, 44, 119, 148, 112, 144, 85, 145, 107, 94, 15, 20, 48, 135, 10, 35, 55, 69, 13, 120, 67, 12, 29, 38, 84, 30, 17, 40, 52, 95, 77, 136, 73, 11, 59, 141, 33, 89, 128, 56, 79, 49, 98, 125, 74, 36, 108, 54, 134, 57, 41, 86, 127, 151, 137, 90, 47, 114, 8, 88, 105, 126, 81, 91, 104, 21, 50, 1, 121, 19, 5, 92, 23, 75, 129, 71, 109, 32, 124, 24, 7, 76, 111, 51, 142, 87, 78, 14, 99, 62, 18, 39])

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
    s = 0.0
    for v in x; s += v^2; end
    return -exp(-0.5 * s)
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
        f_star = -1.0,
        f_star_justification = "E",
        n_minimizers_documented = 1,
        description = "DFO Benchmark Suite v2.1.0 core (Julia): CK adaptation of exponencial (F = f o T).",
        x0_idx = Float64[0.05124221086502079, 42.0, 0.8142949543936, 145.0],
        x0_alt_idx = Vector{Float64}[Float64[0.4483392751216888, 46.0, 0.0771648987136001, 148.0], Float64[-0.3803850328922272, 36.0, -0.41166871715840003, 134.0], Float64[0.6653861176967621, 45.0, 0.4030539759616001, 140.0], Float64[-0.1633381903171539, 48.0, -0.09819446190079995, 147.0], Float64[0.24115819811820982, 37.0, 0.7165282312192001, 132.0], Float64[-0.5875661098957062, 42.0, 0.22769461534720004, 136.0], Float64[0.9317617881298066, 51.0, -0.5094354403328, 143.0], Float64[0.10303748011589053, 38.0, 0.5535836925952, 149.0], Float64[0.5026009857654572, 42.0, -0.1835463630848, 135.0]],
        suite_release = "2.1.0",
        base_function = "exponencial",
        m = 151,
        pattern = "CK",
        strategy = "baseline",
        scale_factors_C = Float64[1.0, 1.0],
        lb_orig = Float64[-0.68375732421875, -0.7109942590052329, -0.6165132799999999, -1.0],
        ub_orig = Float64[1.0, 1.0, 1.0, 0.735551853394419],
        lb_C_work = Float64[-0.68375732421875, -0.7109942590052329, -0.6165132799999999, -1.0],
        ub_C_work = Float64[1.0, 1.0, 1.0, 0.735551853394419],
        ProblemData = _problem_data(),
        func_idx = _objective,
        func_F = xr -> _objective(_cell2idx(xr)),
        idx2cell = _idx2cell,
        cell2idx = _cell2idx,
        decode = xr -> _decode(_cell2idx(xr)))
end

end # module
