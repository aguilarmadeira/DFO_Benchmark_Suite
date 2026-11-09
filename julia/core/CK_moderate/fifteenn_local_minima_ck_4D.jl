# fifteenn_local_minima_ck_4D  DFO Benchmark Suite v2.1 core, Julia version (self-contained; Base only).
#FIFTEENN_LOCAL_MINIMA_CK_4D  DFO Benchmark Suite v2.1.0 core: fifteenn_local_minima, dimension 4, pattern CK, strategy moderate.
#
# Copyright (C) 2026 J. F. A. Madeira
#
# Objective F(z) = f(T(z)), where f is the original function (embedded below, translated from the MATLAB original) and T acts on the finite
# coordinates only: uniform grid of 151 levels on the interval of the domain, piecewise-affine map with |T - Id| <=
# half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
# Global optimum of F on this domain: f* = 1.349783804395672e-32 (E: exact by structure); the optimal solution is listed in
# reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
# Variables: x1:C x2:K x3:C x4:K
#
# Julia version of the MATLAB wrapper fifteenn_local_minima_ck_4D.m: same instance, same values (level indices 1..m as in MATLAB).
# Usage:  include("fifteenn_local_minima_ck_4D.jl");  PB = fifteenn_local_minima_ck_4D.problem();  f = PB.func_idx(PB.x0_idx)
# Copyright (C) 2026 J. F. A. Madeira. LGPL-3.0-or-later.
module fifteenn_local_minima_ck_4D

const NAME = "fifteenn_local_minima_ck_4D"
const N = 4
const TYPES = Char['C', 'K', 'C', 'K']
const CARD = Float64[Inf, 151.0, Inf, 151.0]
const LB = Float64[-10.0, -10.0, -10.0, -10.0]
const UB = Float64[7.275073242187499, 6.176609256718996, 6.949132799999999, 6.95551853394419]
const LBW = Float64[-10.0, -10.0, -1000.0000000000001, -10.0]
const UBW = Float64[7.275073242187499, 6.176609256718996, 694.9132800000001, 6.95551853394419]
const TB = Float64[NaN, 1.0000942945689175, NaN, 0.964568651950577]
const TI = Float64[NaN, 1.0, NaN, 1.0]
const CW = Float64[1.0, NaN, 100.00000000000011, NaN]
const CX = Float64[1.0, NaN, 1.0, NaN]
const M = 151
const PERM = Dict{Int,Vector{Int}}(2 => Int[137, 48, 20, 33, 103, 71, 47, 84, 89, 58, 12, 144, 53, 61, 125, 45, 11, 62, 117, 116, 145, 115, 55, 135, 108, 110, 14, 82, 91, 114, 59, 32, 3, 72, 52, 97, 16, 24, 26, 83, 35, 90, 65, 67, 9, 140, 86, 77, 54, 41, 87, 131, 31, 34, 111, 106, 10, 120, 30, 151, 5, 75, 134, 80, 85, 138, 66, 76, 18, 44, 78, 6, 130, 95, 68, 88, 27, 49, 142, 94, 127, 143, 132, 43, 74, 64, 128, 126, 92, 107, 21, 37, 69, 148, 113, 96, 57, 19, 129, 139, 22, 150, 147, 50, 133, 4, 40, 136, 46, 119, 81, 1, 42, 38, 101, 63, 15, 36, 118, 122, 93, 23, 39, 56, 104, 141, 100, 149, 28, 29, 8, 79, 17, 105, 60, 7, 123, 124, 51, 98, 13, 102, 121, 2, 109, 112, 99, 25, 146, 73, 70], 4 => Int[9, 100, 24, 51, 69, 14, 34, 123, 62, 85, 99, 23, 55, 89, 116, 35, 12, 113, 39, 119, 82, 131, 67, 5, 104, 42, 112, 88, 22, 38, 139, 111, 4, 124, 107, 97, 93, 25, 27, 31, 84, 140, 102, 103, 11, 129, 149, 91, 78, 132, 3, 75, 66, 36, 71, 87, 13, 144, 73, 86, 96, 142, 108, 98, 45, 109, 143, 106, 151, 90, 53, 10, 33, 122, 72, 17, 94, 121, 77, 30, 148, 74, 127, 135, 79, 70, 7, 19, 92, 47, 29, 44, 76, 147, 137, 145, 130, 105, 37, 60, 117, 83, 125, 133, 28, 128, 52, 134, 59, 49, 95, 15, 138, 118, 2, 80, 32, 54, 18, 43, 141, 46, 64, 81, 6, 40, 1, 16, 150, 65, 48, 120, 61, 20, 110, 63, 50, 146, 126, 41, 114, 68, 115, 26, 8, 57, 101, 56, 58, 136, 21])

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
    d = length(x)
    s = 0.0
    for i in 1:d-1; s += ((x[i] - 1)^2) * (1 + 10 * sin(3 * pi * x[i+1])^2); end
    f = 0.1 * s
    return f + 0.1 * (sin(3 * pi * x[1])^2 + (x[d] - 1)^2 * (1 + sin(2 * pi * x[d])^2))
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
        f_star = 1.3497838043956716e-32,
        f_star_justification = "E",
        n_minimizers_documented = 1,
        description = "DFO Benchmark Suite v2.1.0 core (Julia): CK adaptation of fifteenn_local_minima (F = f o T).",
        x0_idx = Float64[-7.857486033439637, 142.0, -448.9497944064, 40.0],
        x0_alt_idx = Vector{Float64}[Float64[5.933893239498138, 146.0, 405.2864987136003, 23.0], Float64[-2.5686818718910223, 136.0, -107.25527715840008, 30.0], Float64[1.6826056838035583, 145.0, -880.1357328384001, 36.0], Float64[-6.819969427585602, 148.0, 221.42230609920023, 43.0], Float64[3.8082494616508473, 137.0, -551.4581495808001, 28.0], Float64[-4.694325649738312, 142.0, 563.1168233472002, 32.0], Float64[-0.4430380940437324, 151.0, -209.76363233280006, 38.0], Float64[-8.945613205432892, 138.0, -722.3054082048001, 24.0], Float64[6.4906094670295715, 142.0, 131.93088491520018, 31.0]],
        suite_release = "2.1.0",
        base_function = "fifteenn_local_minima",
        m = 151,
        pattern = "CK",
        strategy = "moderate",
        scale_factors_C = Float64[1.0, 100.0],
        lb_orig = Float64[-10.0, -10.0, -10.0, -10.0],
        ub_orig = Float64[7.275073242187499, 6.176609256718996, 6.949132799999999, 6.95551853394419],
        lb_C_work = Float64[-10.0, -10.0, -1000.0000000000001, -10.0],
        ub_C_work = Float64[7.275073242187499, 6.176609256718996, 694.9132800000001, 6.95551853394419],
        ProblemData = _problem_data(),
        func_idx = _objective,
        func_F = xr -> _objective(_cell2idx(xr)),
        idx2cell = _idx2cell,
        cell2idx = _cell2idx,
        decode = xr -> _decode(_cell2idx(xr)))
end

end # module
