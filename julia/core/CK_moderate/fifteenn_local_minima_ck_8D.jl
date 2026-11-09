# fifteenn_local_minima_ck_8D  DFO Benchmark Suite v2.1 core, Julia version (self-contained; Base only).
#FIFTEENN_LOCAL_MINIMA_CK_8D  DFO Benchmark Suite v2.1.0 core: fifteenn_local_minima, dimension 8, pattern CK, strategy moderate.
#
# Copyright (C) 2026 J. F. A. Madeira
#
# Objective F(z) = f(T(z)), where f is the original function (embedded below, translated from the MATLAB original) and T acts on the finite
# coordinates only: uniform grid of 25 levels on the interval of the domain, piecewise-affine map with |T - Id| <=
# half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
# Global optimum of F on this domain: f* = 1.349783804395672e-32 (E: exact by structure); the optimal solution is listed in
# reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
# Variables: x1:C x2:K x3:C x4:K x5:C x6:K x7:C x8:K
#
# Julia version of the MATLAB wrapper fifteenn_local_minima_ck_8D.m: same instance, same values (level indices 1..m as in MATLAB).
# Usage:  include("fifteenn_local_minima_ck_8D.jl");  PB = fifteenn_local_minima_ck_8D.problem();  f = PB.func_idx(PB.x0_idx)
# Copyright (C) 2026 J. F. A. Madeira. LGPL-3.0-or-later.
module fifteenn_local_minima_ck_8D

const NAME = "fifteenn_local_minima_ck_8D"
const N = 8
const TYPES = Char['C', 'K', 'C', 'K', 'C', 'K', 'C', 'K']
const CARD = Float64[Inf, 25.0, Inf, 25.0, Inf, 25.0, Inf, 25.0]
const LB = Float64[-10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0]
const UB = Float64[6.9250732421875, 6.487720367830107, 6.3891328, 6.55551853394419, 7.396366368417458, 6.164258954518399, 6.225670190730475, 6.853174852863314]
const LBW = Float64[-10.0, -10.0, -999.9999999999999, -10.0, -0.09999999999999999, -10.0, -10.0, -10.0]
const UBW = Float64[6.9250732421875005, 6.487720367830107, 638.9132799999999, 6.55551853394419, 0.07396366368417458, 6.164258954518399, 6.225670190730476, 6.853174852863314]
const TB = Float64[NaN, 0.9918135785534048, NaN, 1.0370123559627942, NaN, 0.7761726363455992, NaN, 1.2354499019088774]
const TI = Float64[NaN, 1.0, NaN, 1.0, NaN, 1.0, NaN, 1.0]
const CW = Float64[1.0, NaN, 99.99999999999989, NaN, 0.010000000000000009, NaN, 1.0, NaN]
const CX = Float64[1.0, NaN, 1.0, NaN, 1.0, NaN, 1.0, NaN]
const M = 25
const PERM = Dict{Int,Vector{Int}}(2 => Int[12, 9, 23, 19, 4, 25, 11, 10, 16, 21, 22, 15, 1, 14, 18, 7, 6, 5, 24, 13, 8, 20, 17, 3, 2], 4 => Int[25, 17, 20, 11, 19, 23, 24, 16, 15, 2, 21, 12, 14, 18, 1, 22, 5, 8, 9, 7, 6, 3, 4, 13, 10], 6 => Int[3, 21, 17, 20, 13, 8, 9, 10, 5, 16, 18, 25, 14, 2, 23, 19, 24, 11, 15, 1, 7, 4, 12, 22, 6], 8 => Int[20, 18, 6, 17, 4, 2, 13, 19, 25, 10, 5, 12, 22, 16, 23, 11, 24, 14, 1, 15, 7, 9, 21, 3, 8])

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
        cardinality = Float64[Inf, 25.0, Inf, 25.0, Inf, 25.0, Inf, 25.0],
        f_star = 1.3497838043956716e-32,
        f_star_justification = "E",
        n_minimizers_documented = 1,
        description = "DFO Benchmark Suite v2.1.0 core (Julia): CK adaptation of fifteenn_local_minima (F = f o T).",
        x0_idx = Float64[-3.669625926017761, 11.0, 188.40879759360007, 14.0, 0.0016857478183980729, 8.0, 4.775829113235934, 16.0],
        x0_alt_idx = Vector{Float64}[Float64[0.3219807028770454, 9.0, -558.9356580863998, 11.0, -0.037132590359062384, 5.0, 0.22814992483050744, 24.0], Float64[-8.008328783512116, 10.0, 518.8137148415998, 12.0, -0.09163510557792097, 2.0, -4.31952926357492, 6.0], Float64[5.6771796584129355, 12.0, -228.5307408384, 13.0, 0.042072503699998307, 22.0, -9.818357040405012, 14.0], Float64[-2.653129827976226, 10.0, -736.7249707008, 15.0, -0.012430011518860265, 19.0, 1.8034897744088578, 22.0], Float64[1.412854564189912, 10.0, 89.28732241920022, 12.0, -0.06835836101420961, 16.0, -2.74418941399657, 6.0], Float64[-6.91745492219925, 11.0, -406.32005345280004, 13.0, 0.06534924826370968, 11.0, -7.291868602401997, 13.0], Float64[3.594602286815645, 9.0, 419.69223966720017, 14.0, 0.010846733044851078, 8.0, 3.37882962398721, 20.0], Float64[-4.7357071995735165, 11.0, -75.91513620479998, 12.0, -0.027971605132609365, 3.0, -1.1688495644182186, 4.0], Float64[-0.7193079829215989, 11.0, -823.2595918847999, 13.0, -0.08247412035146795, 24.0, -5.716528752823646, 12.0]],
        suite_release = "2.1.0",
        base_function = "fifteenn_local_minima",
        m = 25,
        pattern = "CK",
        strategy = "moderate",
        scale_factors_C = Float64[1.0, 100.0, 0.01, 1.0],
        lb_orig = Float64[-10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0],
        ub_orig = Float64[6.9250732421875, 6.487720367830107, 6.3891328, 6.55551853394419, 7.396366368417458, 6.164258954518399, 6.225670190730475, 6.853174852863314],
        lb_C_work = Float64[-10.0, -10.0, -999.9999999999999, -10.0, -0.09999999999999999, -10.0, -10.0, -10.0],
        ub_C_work = Float64[6.9250732421875005, 6.487720367830107, 638.9132799999999, 6.55551853394419, 0.07396366368417458, 6.164258954518399, 6.225670190730476, 6.853174852863314],
        ProblemData = _problem_data(),
        func_idx = _objective,
        func_F = xr -> _objective(_cell2idx(xr)),
        idx2cell = _idx2cell,
        cell2idx = _cell2idx,
        decode = xr -> _decode(_cell2idx(xr)))
end

end # module
