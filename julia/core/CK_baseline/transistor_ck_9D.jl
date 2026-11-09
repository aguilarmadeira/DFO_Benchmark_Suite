# transistor_ck_9D  DFO Benchmark Suite v2.1 core, Julia version (self-contained; Base only).
#TRANSISTOR_CK_9D  DFO Benchmark Suite v2.1.0 core: transistor, dimension 9, pattern CK, strategy baseline.
#
# Copyright (C) 2026 J. F. A. Madeira
#
# Objective F(z) = f(T(z)), where f is the original function (embedded below, translated from the MATLAB original) and T acts on the finite
# coordinates only: uniform grid of 25 levels on the interval of the domain, piecewise-affine map with |T - Id| <=
# half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
# Global optimum of F on this domain: f* = 1.60049326654736e-24 (E: exact by structure); the optimal solution is listed in
# reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
# Variables: x1:C x2:K x3:C x4:K x5:C x6:K x7:C x8:K x9:C
#
# Julia version of the MATLAB wrapper transistor_ck_9D.m: same instance, same values (level indices 1..m as in MATLAB).
# Usage:  include("transistor_ck_9D.jl");  PB = transistor_ck_9D.problem();  f = PB.func_idx(PB.x0_idx)
# Copyright (C) 2026 J. F. A. Madeira. LGPL-3.0-or-later.
module transistor_ck_9D

const NAME = "transistor_ck_9D"
const N = 9
const TYPES = Char['C', 'K', 'C', 'K', 'C', 'K', 'C', 'K', 'C']
const CARD = Float64[Inf, 25.0, Inf, 25.0, Inf, 25.0, Inf, 25.0, Inf]
const LB = Float64[-10.0, -10.0, -10.0, -10.0, -7.350085376681921, -6.992661321382304, -10.0, -10.0, -10.0]
const UB = Float64[6.8211669921875, 7.03312503175329, 6.7027328, 6.180008329862558, 10.0, 10.0, 6.7876079069934505, 6.403313357018439, 7.316462562669516]
const LBW = Float64[-10.0, -10.0, -10.0, -10.0, -7.350085376681921, -6.992661321382304, -10.0, -10.0, -10.0]
const UBW = Float64[6.8211669921875, 7.03312503175329, 6.7027328, 6.180008329862558, 10.0, 10.0, 6.7876079069934505, 6.403313357018439, 7.316462562669516]
const TB = Float64[NaN, 0.6457031448458057, NaN, 2.1350062473969214, NaN, 7.8759173348272125, NaN, 0.9355422380122942, NaN]
const TI = Float64[NaN, 0.4499874719815851, NaN, 2.000068541623683, NaN, 7.999692684218543, NaN, 0.9999877234567898, NaN]
const CW = Float64[0.8999999526168558, NaN, 1.0000064824653876, NaN, 7.999971440507853, NaN, 5.000031275929672, NaN, 2.00005248348619]
const CX = Float64[0.8999999526168565, NaN, 1.000006482465387, NaN, 7.999971440507853, NaN, 5.000031275929672, NaN, 2.00005248348619]
const M = 25
const PERM = Dict{Int,Vector{Int}}(2 => Int[2, 5, 22, 21, 12, 11, 25, 8, 17, 13, 10, 24, 3, 20, 1, 9, 15, 16, 4, 19, 23, 14, 7, 18, 6], 4 => Int[2, 23, 17, 9, 25, 5, 22, 12, 15, 21, 4, 3, 19, 20, 11, 6, 14, 1, 10, 13, 24, 18, 16, 7, 8], 6 => Int[17, 11, 5, 16, 13, 24, 3, 6, 8, 1, 12, 9, 18, 19, 22, 23, 14, 4, 15, 21, 7, 20, 25, 10, 2], 8 => Int[18, 9, 4, 6, 21, 17, 10, 16, 24, 5, 25, 15, 23, 1, 3, 14, 13, 7, 20, 22, 8, 2, 11, 19, 12])

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

const _G = [0.485 0.369 5.2095 23.3037 28.5132; 0.752 1.254 10.0677 101.779 111.8467;
            0.869 0.703 22.9274 111.461 134.3884; 0.982 1.455 20.2153 191.267 211.4823]
function _orig(x::Vector{Float64})
    s = 0.0
    for k in 1:4
        g = _G[k, :]
        alfa = (1 - x[1] * x[2]) * x[3] * (exp(x[5] * (g[1] - g[3] * x[7] * 0.001 - g[5] * x[8] * 0.001)) - 1) - g[5] + g[4] * x[2]
        beta = (1 - x[1] * x[2]) * x[4] * (exp(x[6] * (g[1] - g[2] - g[3] * x[7] * 0.001 + g[4] * x[9] * 0.001)) - 1) - g[5] * x[1] + g[4]
        s += alfa^2 + beta^2
    end
    return (x[1] * x[3] - x[2] * x[4])^2 + s
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
        cardinality = Float64[Inf, 25.0, Inf, 25.0, Inf, 25.0, Inf, 25.0, Inf],
        f_star = 1.6004932665473596e-24,
        f_star_justification = "E",
        n_minimizers_documented = 1,
        description = "DFO Benchmark Suite v2.1.0 core (Julia): CK adaptation of transistor (F = f o T).",
        x0_idx = Float64[-3.4210181832313538, 2.0, -1.6299265392639999, 20.0, 3.36499364633276, 18.0, -1.4507235938738816, 24.0, -6.985594829642967],
        x0_alt_idx = Vector{Float64}[Float64[0.9156889319419861, 4.0, -9.887757635584, 21.0, -2.0707581208170245, 15.0, -6.155900896526026, 6.0, 5.878998340733439], Float64[-7.560602247714996, 2.0, 1.7373443932160004, 19.0, -5.942264775117948, 10.0, 4.884351728651227, 14.0, 1.459863282207193], Float64[3.084042489528656, 2.0, -6.648763691008, 20.0, 5.828679699069712, 7.0, 0.17917442599908284, 22.0, -2.9592717763190537], Float64[-5.490810215473175, 3.0, 4.976338337791999, 21.0, 1.957173044768787, 4.0, -4.526002876653062, 6.0, -7.378406834845299], Float64[-1.1541031002998352, 1.0, -3.2814927585280005, 18.0, -3.6207828235490744, 23.0, -9.231180179305206, 13.0, 5.486186335531107], Float64[-9.630394279956818, 3.0, -7.691014217728, 19.0, 9.714406763487442, 20.0, 2.7931618163613834, 21.0, 1.0670512770048592], Float64[6.533695876598358, 3.0, 0.08577817395200071, 20.0, 4.278654996337661, 17.0, -2.8961048567800978, 4.0, -3.3520837815213858], Float64[-2.090437591075897, 1.0, -4.323743285248001, 21.0, 0.4071483420367352, 13.0, -7.601282159432241, 12.0, -7.771218840047632], Float64[2.2462695240974426, 2.0, 3.453049106432001, 19.0, -5.028603425113047, 10.0, 4.423059836234348, 19.0, 4.341907885593876]],
        suite_release = "2.1.0",
        base_function = "transistor",
        m = 25,
        pattern = "CK",
        strategy = "baseline",
        scale_factors_C = Float64[1.0, 1.0, 1.0, 1.0, 1.0],
        lb_orig = Float64[-10.0, -10.0, -10.0, -10.0, -7.350085376681921, -6.992661321382304, -10.0, -10.0, -10.0],
        ub_orig = Float64[6.8211669921875, 7.03312503175329, 6.7027328, 6.180008329862558, 10.0, 10.0, 6.7876079069934505, 6.403313357018439, 7.316462562669516],
        lb_C_work = Float64[-10.0, -10.0, -10.0, -10.0, -7.350085376681921, -6.992661321382304, -10.0, -10.0, -10.0],
        ub_C_work = Float64[6.8211669921875, 7.03312503175329, 6.7027328, 6.180008329862558, 10.0, 10.0, 6.7876079069934505, 6.403313357018439, 7.316462562669516],
        ProblemData = _problem_data(),
        func_idx = _objective,
        func_F = xr -> _objective(_cell2idx(xr)),
        idx2cell = _idx2cell,
        cell2idx = _cell2idx,
        decode = xr -> _decode(_cell2idx(xr)))
end

end # module
