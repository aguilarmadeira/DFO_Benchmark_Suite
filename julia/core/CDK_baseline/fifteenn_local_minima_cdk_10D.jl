# fifteenn_local_minima_cdk_10D  DFO Benchmark Suite v2.1 core, Julia version (self-contained; Base only).
#FIFTEENN_LOCAL_MINIMA_CDK_10D  DFO Benchmark Suite v2.1.0 core: fifteenn_local_minima, dimension 10, pattern CDK, strategy baseline.
#
# Copyright (C) 2026 J. F. A. Madeira
#
# Objective F(z) = f(T(z)), where f is the original function (embedded below, translated from the MATLAB original) and T acts on the finite
# coordinates only: uniform grid of 17 levels on the interval of the domain, piecewise-affine map with |T - Id| <=
# half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
# Global optimum of F on this domain: f* = 1.349783804395672e-32 (E: exact by structure); the optimal solution is listed in
# reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
# Variables: x1:C x2:K x3:D x4:C x5:D x6:K x7:C x8:K x9:D x10:C
#
# Julia version of the MATLAB wrapper fifteenn_local_minima_cdk_10D.m: same instance, same values (level indices 1..m as in MATLAB).
# Usage:  include("fifteenn_local_minima_cdk_10D.jl");  PB = fifteenn_local_minima_cdk_10D.problem();  f = PB.func_idx(PB.x0_idx)
# Copyright (C) 2026 J. F. A. Madeira. LGPL-3.0-or-later.
module fifteenn_local_minima_cdk_10D

const NAME = "fifteenn_local_minima_cdk_10D"
const N = 10
const TYPES = Char['C', 'K', 'D', 'C', 'D', 'K', 'C', 'K', 'D', 'C']
const CARD = Float64[Inf, 17.0, 17.0, Inf, 17.0, 17.0, Inf, 17.0, 17.0, Inf]
const LB = Float64[-10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0]
const UB = Float64[6.2250732421875, 6.02105370116344, 6.1091328, 6.35551853394419, 7.26909364114473, 6.0565666468260915, 6.143317249554005, 6.779490642336999, 7.1391468726884195, 7.334618065521342]
const LBW = Float64[-10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0]
const UBW = Float64[6.2250732421875, 6.02105370116344, 6.1091328, 6.35551853394419, 7.26909364114473, 6.0565666468260915, 6.143317249554005, 6.779490642336999, 7.1391468726884195, 7.334618065521342]
const TB = Float64[NaN, 1.0144744195498632, 1.0750288000000001, NaN, 0.793183525715456, 1.0388895696929374, NaN, 0.48718165146062375, 0.7119667954302642, NaN]
const TI = Float64[NaN, 1.0, 1.0, NaN, 1.0, 1.0, NaN, 1.0, 1.0, NaN]
const CW = Float64[1.0, NaN, NaN, 1.0, NaN, NaN, 1.0, NaN, NaN, 1.0]
const CX = Float64[1.0, NaN, NaN, 1.0, NaN, NaN, 1.0, NaN, NaN, 1.0]
const M = 17
const PERM = Dict{Int,Vector{Int}}(2 => Int[8, 12, 1, 9, 13, 7, 14, 6, 5, 4, 16, 2, 15, 11, 3, 17, 10], 6 => Int[8, 9, 14, 11, 3, 16, 10, 15, 4, 12, 6, 7, 2, 1, 17, 5, 13], 8 => Int[11, 17, 6, 16, 10, 12, 5, 1, 15, 8, 9, 13, 2, 3, 4, 7, 14])

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
        cardinality = Float64[Inf, 17.0, 17.0, Inf, 17.0, 17.0, Inf, 17.0, 17.0, Inf],
        f_star = 1.3497838043956716e-32,
        f_star_justification = "E",
        n_minimizers_documented = 1,
        description = "DFO Benchmark Suite v2.1.0 core (Julia): CDK adaptation of fifteenn_local_minima (F = f o T).",
        x0_idx = Float64[4.181094288825987, 13.0, 16.0, 1.3146673406419431, 12.0, 7.0, 5.650441697460966, 12.0, 7.0, -7.484635969745376],
        x0_alt_idx = Vector{Float64}[Float64[-4.161191904544831, 12.0, 9.0, -0.5858697984544463, 8.0, 5.0, 1.1258441292468682, 17.0, 2.0, -4.4546848928944405], Float64[0.02183479070663452, 12.0, 4.0, 0.12938611410871026, 3.0, 3.0, -3.3987534389672307, 5.0, 14.0, -1.4247338160435046], Float64[-8.15408102273941, 13.0, 12.0, 0.8446420266718668, 16.0, 16.0, -8.869672067199964, 11.0, 10.0, 1.6052172608074322], Float64[5.155549371242522, 12.0, 7.0, 1.5598979392350234, 11.0, 14.0, 2.693188384902731, 16.0, 5.0, 4.635168337658367], Float64[-3.115435230731965, 13.0, 15.0, -0.013665068403920344, 5.0, 12.0, -1.8314091833113668, 5.0, 1.0, -9.648886738924617], Float64[1.0675914645195004, 13.0, 10.0, 0.37461671270179053, 2.0, 9.0, -6.356006751525465, 10.0, 13.0, -6.61893566207368], Float64[-7.108324348926544, 12.0, 3.0, 1.0898726252649489, 13.0, 7.0, 4.260532640558596, 14.0, 9.0, -3.588984585222744], Float64[3.159104812145232, 13.0, 13.0, -0.4836903823739984, 9.0, 4.0, -0.26406492765550027, 4.0, 5.0, -0.5590335083718063], Float64[-5.1594141840934755, 13.0, 6.0, 0.23156553018915993, 4.0, 2.0, -4.788662495869601, 9.0, 17.0, 2.4709175684791287]],
        suite_release = "2.1.0",
        base_function = "fifteenn_local_minima",
        m = 17,
        pattern = "CDK",
        strategy = "baseline",
        scale_factors_C = Float64[1.0, 1.0, 1.0, 1.0],
        lb_orig = Float64[-10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0],
        ub_orig = Float64[6.2250732421875, 6.02105370116344, 6.1091328, 6.35551853394419, 7.26909364114473, 6.0565666468260915, 6.143317249554005, 6.779490642336999, 7.1391468726884195, 7.334618065521342],
        lb_C_work = Float64[-10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0],
        ub_C_work = Float64[6.2250732421875, 6.02105370116344, 6.1091328, 6.35551853394419, 7.26909364114473, 6.0565666468260915, 6.143317249554005, 6.779490642336999, 7.1391468726884195, 7.334618065521342],
        ProblemData = _problem_data(),
        func_idx = _objective,
        func_F = xr -> _objective(_cell2idx(xr)),
        idx2cell = _idx2cell,
        cell2idx = _cell2idx,
        decode = xr -> _decode(_cell2idx(xr)))
end

end # module
