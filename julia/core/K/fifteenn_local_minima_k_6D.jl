# fifteenn_local_minima_k_6D  DFO Benchmark Suite v2.1 core, Julia version (self-contained; Base only).
#FIFTEENN_LOCAL_MINIMA_K_6D  DFO Benchmark Suite v2.1.0 core: fifteenn_local_minima, dimension 6, pattern K, strategy baseline.
#
# Copyright (C) 2026 J. F. A. Madeira
#
# Objective F(z) = f(T(z)), where f is the original function (embedded below, translated from the MATLAB original) and T acts on the finite
# coordinates only: uniform grid of 17 levels on the interval of the domain, piecewise-affine map with |T - Id| <=
# half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
# Global optimum of F on this domain: f* = 1.349783804395672e-32 (E: exact by structure); the optimal solution is listed in
# reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
# Variables: x1:K x2:K x3:K x4:K x5:K x6:K
#
# Julia version of the MATLAB wrapper fifteenn_local_minima_k_6D.m: same instance, same values (level indices 1..m as in MATLAB).
# Usage:  include("fifteenn_local_minima_k_6D.jl");  PB = fifteenn_local_minima_k_6D.problem();  f = PB.func_idx(PB.x0_idx)
# Copyright (C) 2026 J. F. A. Madeira. LGPL-3.0-or-later.
module fifteenn_local_minima_k_6D

const NAME = "fifteenn_local_minima_k_6D"
const N = 6
const TYPES = Char['K', 'K', 'K', 'K', 'K', 'K']
const CARD = Float64[17.0, 17.0, 17.0, 17.0, 17.0, 17.0]
const LB = Float64[-10.0, -10.0, -10.0, -10.0, -10.0, -10.0]
const UB = Float64[6.5750732421875, 6.954387034496774, 6.6691328, 6.75551853394419, 6.008988457072604, 6.271951262210706]
const LBW = Float64[-10.0, -10.0, -10.0, -10.0, -10.0, -10.0]
const UBW = Float64[6.5750732421875, 6.954387034496774, 6.6691328, 6.75551853394419, 6.008988457072604, 6.271951262210706]
const TB = Float64[1.3953628540039062, 0.5964918965604848, 1.4600288, 1.5194189920866314, 1.0061795642374154, 1.1869664927698604]
const TI = Float64[1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
const CW = Float64[NaN, NaN, NaN, NaN, NaN, NaN]
const CX = Float64[NaN, NaN, NaN, NaN, NaN, NaN]
const M = 17
const PERM = Dict{Int,Vector{Int}}(1 => Int[14, 15, 12, 4, 17, 1, 7, 11, 10, 13, 9, 5, 6, 8, 16, 2, 3], 2 => Int[17, 15, 12, 8, 4, 13, 10, 9, 14, 2, 5, 3, 7, 1, 6, 11, 16], 3 => Int[17, 3, 11, 2, 13, 10, 14, 4, 5, 16, 12, 1, 15, 9, 7, 8, 6], 4 => Int[14, 9, 5, 10, 2, 16, 15, 13, 1, 3, 8, 17, 7, 4, 12, 6, 11], 5 => Int[10, 17, 6, 4, 16, 8, 3, 2, 11, 14, 5, 9, 7, 12, 13, 15, 1], 6 => Int[6, 2, 17, 11, 13, 4, 12, 10, 9, 14, 1, 8, 16, 5, 3, 7, 15])

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
        cardinality = Float64[17.0, 17.0, 17.0, 17.0, 17.0, 17.0],
        f_star = 1.3497838043956716e-32,
        f_star_justification = "E",
        n_minimizers_documented = 1,
        description = "DFO Benchmark Suite v2.1.0 core (Julia): K adaptation of fifteenn_local_minima (F = f o T).",
        x0_idx = Float64[11.0, 2.0, 9.0, 7.0, 9.0, 5.0],
        x0_alt_idx = Vector{Float64}[Float64[3.0, 1.0, 2.0, 6.0, 5.0, 3.0], Float64[16.0, 2.0, 13.0, 6.0, 16.0, 17.0], Float64[8.0, 3.0, 5.0, 7.0, 13.0, 13.0], Float64[12.0, 1.0, 16.0, 8.0, 8.0, 11.0], Float64[4.0, 2.0, 8.0, 6.0, 2.0, 9.0], Float64[14.0, 2.0, 4.0, 7.0, 15.0, 6.0], Float64[6.0, 1.0, 12.0, 7.0, 10.0, 4.0], Float64[10.0, 2.0, 7.0, 6.0, 6.0, 1.0], Float64[2.0, 2.0, 15.0, 6.0, 1.0, 15.0]],
        suite_release = "2.1.0",
        base_function = "fifteenn_local_minima",
        m = 17,
        pattern = "K",
        strategy = "baseline",
        scale_factors_C = Float64[],
        lb_orig = Float64[-10.0, -10.0, -10.0, -10.0, -10.0, -10.0],
        ub_orig = Float64[6.5750732421875, 6.954387034496774, 6.6691328, 6.75551853394419, 6.008988457072604, 6.271951262210706],
        lb_C_work = Float64[-10.0, -10.0, -10.0, -10.0, -10.0, -10.0],
        ub_C_work = Float64[6.5750732421875, 6.954387034496774, 6.6691328, 6.75551853394419, 6.008988457072604, 6.271951262210706],
        ProblemData = _problem_data(),
        func_idx = _objective,
        func_F = xr -> _objective(_cell2idx(xr)),
        idx2cell = _idx2cell,
        cell2idx = _cell2idx,
        decode = xr -> _decode(_cell2idx(xr)))
end

end # module
