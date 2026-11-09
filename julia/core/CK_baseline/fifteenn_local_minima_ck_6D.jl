# fifteenn_local_minima_ck_6D  DFO Benchmark Suite v2.1 core, Julia version (self-contained; Base only).
#FIFTEENN_LOCAL_MINIMA_CK_6D  DFO Benchmark Suite v2.1.0 core: fifteenn_local_minima, dimension 6, pattern CK, strategy baseline.
#
# Copyright (C) 2026 J. F. A. Madeira
#
# Objective F(z) = f(T(z)), where f is the original function (embedded below, translated from the MATLAB original) and T acts on the finite
# coordinates only: uniform grid of 41 levels on the interval of the domain, piecewise-affine map with |T - Id| <=
# half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
# Global optimum of F on this domain: f* = 1.349783804395672e-32 (E: exact by structure); the optimal solution is listed in
# reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
# Variables: x1:C x2:K x3:C x4:K x5:C x6:K
#
# Julia version of the MATLAB wrapper fifteenn_local_minima_ck_6D.m: same instance, same values (level indices 1..m as in MATLAB).
# Usage:  include("fifteenn_local_minima_ck_6D.jl");  PB = fifteenn_local_minima_ck_6D.problem();  f = PB.func_idx(PB.x0_idx)
# Copyright (C) 2026 J. F. A. Madeira. LGPL-3.0-or-later.
module fifteenn_local_minima_ck_6D

const NAME = "fifteenn_local_minima_ck_6D"
const N = 6
const TYPES = Char['C', 'K', 'C', 'K', 'C', 'K']
const CARD = Float64[Inf, 41.0, Inf, 41.0, Inf, 41.0]
const LB = Float64[-10.0, -10.0, -10.0, -10.0, -10.0, -10.0]
const UB = Float64[6.5750732421875, 6.954387034496774, 6.6691328, 6.75551853394419, 6.008988457072604, 6.271951262210706]
const LBW = Float64[-10.0, -10.0, -10.0, -10.0, -10.0, -10.0]
const UBW = Float64[6.5750732421875, 6.954387034496774, 6.6691328, 6.75551853394419, 6.008988457072604, 6.271951262210706]
const TB = Float64[NaN, 1.0203515724229053, NaN, 0.8910870470637242, NaN, 0.983567101992227]
const TI = Float64[NaN, 1.0, NaN, 1.0, NaN, 1.0]
const CW = Float64[1.0, NaN, 1.0, NaN, 1.0, NaN]
const CX = Float64[1.0, NaN, 1.0, NaN, 1.0, NaN]
const M = 41
const PERM = Dict{Int,Vector{Int}}(2 => Int[35, 20, 17, 4, 13, 38, 11, 19, 31, 30, 6, 16, 25, 40, 12, 1, 27, 37, 8, 21, 24, 29, 14, 9, 32, 23, 26, 34, 22, 10, 36, 33, 3, 5, 41, 7, 18, 2, 15, 28, 39], 4 => Int[35, 30, 16, 34, 24, 38, 7, 2, 36, 17, 14, 3, 22, 26, 20, 5, 33, 4, 40, 23, 15, 6, 21, 37, 12, 31, 11, 19, 8, 28, 39, 41, 1, 10, 25, 32, 18, 9, 29, 13, 27], 6 => Int[38, 2, 28, 7, 15, 41, 26, 14, 31, 27, 10, 23, 30, 11, 17, 8, 4, 9, 40, 35, 25, 13, 39, 3, 21, 32, 29, 5, 24, 34, 20, 37, 33, 16, 1, 22, 19, 12, 6, 18, 36])

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
        cardinality = Float64[Inf, 41.0, Inf, 41.0, Inf, 41.0],
        f_star = 1.3497838043956716e-32,
        f_star_justification = "E",
        n_minimizers_documented = 1,
        description = "DFO Benchmark Suite v2.1.0 core (Julia): CK adaptation of fifteenn_local_minima (F = f o T).",
        x0_idx = Float64[0.34323418140411377, 5.0, -1.2467049840639994, 17.0, -2.097741986253417, 10.0],
        x0_alt_idx = Vector{Float64}[Float64[-8.179008066654205, 2.0, -8.847829540864, 13.0, -5.669995608905983, 5.0], Float64[5.417731702327728, 3.0, 2.113792188416001, 14.0, 5.191098738754846, 40.0], Float64[-2.9345512986183167, 5.0, -5.487332368384, 16.0, 1.6188451161022819, 32.0], Float64[1.3387097716331482, 2.0, 5.346270420992001, 18.0, -3.3967433035816237, 27.0], Float64[-7.110692799091339, 3.0, -2.2548541358080003, 14.0, -8.543543977541104, 22.0], Float64[3.3782207369804382, 4.0, -7.295599894528, 15.0, 3.760885167151063, 14.0], Float64[-4.974062263965607, 1.0, 1.105643036672003, 17.0, -1.25470325253284, 9.0], Float64[-0.7008011937141418, 4.0, -3.9351027220480006, 13.0, -4.826956875185405, 1.0], Float64[-9.19876354932785, 5.0, 4.466140209152, 15.0, -9.842545294869309, 36.0]],
        suite_release = "2.1.0",
        base_function = "fifteenn_local_minima",
        m = 41,
        pattern = "CK",
        strategy = "baseline",
        scale_factors_C = Float64[1.0, 1.0, 1.0],
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
