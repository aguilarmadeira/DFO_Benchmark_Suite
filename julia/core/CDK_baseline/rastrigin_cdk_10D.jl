# rastrigin_cdk_10D  DFO Benchmark Suite v2.1 core, Julia version (self-contained; Base only).
#RASTRIGIN_CDK_10D  DFO Benchmark Suite v2.1.0 core: rastrigin, dimension 10, pattern CDK, strategy baseline.
#
# Copyright (C) 2026 J. F. A. Madeira
#
# Objective F(z) = f(T(z)), where f is the original function (embedded below, translated from the MATLAB original) and T acts on the finite
# coordinates only: uniform grid of 17 levels on the interval of the domain, piecewise-affine map with |T - Id| <=
# half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
# Global optimum of F on this domain: f* = 0 (E: exact by structure); the optimal solution is listed in
# reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
# Variables: x1:C x2:K x3:D x4:C x5:D x6:K x7:C x8:K x9:D x10:C
#
# Julia version of the MATLAB wrapper rastrigin_cdk_10D.m: same instance, same values (level indices 1..m as in MATLAB).
# Usage:  include("rastrigin_cdk_10D.jl");  PB = rastrigin_cdk_10D.problem();  f = PB.func_idx(PB.x0_idx)
# Copyright (C) 2026 J. F. A. Madeira. LGPL-3.0-or-later.
module rastrigin_cdk_10D

const NAME = "rastrigin_cdk_10D"
const N = 10
const TYPES = Char['C', 'K', 'D', 'C', 'D', 'K', 'C', 'K', 'D', 'C']
const CARD = Float64[Inf, 17.0, 17.0, Inf, 17.0, 17.0, Inf, 17.0, 17.0, Inf]
const LB = Float64[-5.12, -5.12, -3.6955815936, -3.210139775093711, -3.775091646745441, -3.3681715906305802, -5.12, -5.12, -5.12, -5.12]
const UB = Float64[3.7836375, 3.2332190011685213, 5.12, 5.12, 5.12, 5.12, 3.6885590546090206, 3.393660981729729, 3.7162186241472837, 3.23626238058141]
const LBW = Float64[-5.12, -5.12, -3.6955815936, -3.210139775093711, -3.775091646745441, -3.3681715906305802, -5.12, -5.12, -5.12, -5.12]
const UBW = Float64[3.7836375, 3.2332190011685213, 5.12, 5.12, 5.12, 5.12, 3.6885590546090206, 3.393660981729729, 3.7162186241472837, 3.23626238058141]
const TB = Float64[NaN, 0.1007618757303268, 0.1612353536000004, NaN, 0.1165109487056899, -0.1851072441441124, NaN, 0.20103811358108015, -0.14962702391715332, NaN]
const TI = Float64[NaN, 0.0, 0.0, NaN, 0.0, 0.0, NaN, 0.0, 0.0, NaN]
const CW = Float64[0.0, NaN, NaN, 0.0, NaN, NaN, 0.0, NaN, NaN, 0.0]
const CX = Float64[0.0, NaN, NaN, 0.0, NaN, NaN, 0.0, NaN, NaN, 0.0]
const M = 17
const PERM = Dict{Int,Vector{Int}}(2 => Int[14, 3, 15, 13, 9, 4, 2, 11, 1, 6, 8, 10, 5, 12, 16, 17, 7], 6 => Int[11, 12, 15, 16, 14, 7, 8, 13, 10, 9, 4, 3, 5, 2, 17, 1, 6], 8 => Int[17, 10, 5, 15, 13, 7, 16, 14, 9, 1, 6, 11, 4, 3, 12, 8, 2])

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
    for v in x; s += v^2 - 10 * cos(2 * pi * v); end
    return 10 * length(x) + s
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
        f_star = 0.0,
        f_star_justification = "E",
        n_minimizers_documented = 1,
        description = "DFO Benchmark Suite v2.1.0 core (Julia): CDK adaptation of rastrigin (F = f o T).",
        x0_idx = Float64[-4.837413848876953, 14.0, 3.0, 1.7303512868527466, 11.0, 1.0, -3.255375245920338, 14.0, 5.0, 2.1436337065203936],
        x0_alt_idx = Vector{Float64}[Float64[3.301067303466797, 15.0, 11.0, 2.0946431137518884, 6.0, 15.0, 3.0538694748549817, 2.0, 16.0, -4.742085472722076], Float64[-1.1855312805175782, 14.0, 6.0, 2.45893494065103, 2.0, 13.0, 0.068676959910138, 7.0, 12.0, -3.2814783740948856], Float64[1.109937762451172, 14.0, 14.0, 1.657492921472918, 13.0, 10.0, -2.400157930828031, 12.0, 8.0, -1.8208712754676957], Float64[-3.4288305725097654, 15.0, 9.0, 2.0217847483720606, 9.0, 8.0, -4.8689928215661995, 1.0, 4.0, -0.3602641768405048], Float64[2.2055025329589846, 14.0, 2.0, 2.3860765752712023, 4.0, 5.0, 1.4402518992091196, 7.0, 16.0, 1.1003429217866865], Float64[-2.281096051025391, 14.0, 13.0, 2.5838349955878783, 16.0, 3.0, -1.0285829915290492, 12.0, 12.0, 2.560950020413877], Float64[0.014372991943359281, 15.0, 5.0, 1.7823929764097661, 11.0, 17.0, -4.013775506473893, 16.0, 8.0, -4.324769158828593], Float64[-4.550480218505859, 14.0, 16.0, 2.146684803308907, 8.0, 13.0, 2.295469214301427, 5.0, 4.0, -2.864162060201403], Float64[2.7532849182128905, 15.0, 9.0, 2.5109766302080496, 3.0, 11.0, -0.17336567643674172, 11.0, 15.0, -1.403554961574212]],
        suite_release = "2.1.0",
        base_function = "rastrigin",
        m = 17,
        pattern = "CDK",
        strategy = "baseline",
        scale_factors_C = Float64[1.0, 1.0, 1.0, 1.0],
        lb_orig = Float64[-5.12, -5.12, -3.6955815936, -3.210139775093711, -3.775091646745441, -3.3681715906305802, -5.12, -5.12, -5.12, -5.12],
        ub_orig = Float64[3.7836375, 3.2332190011685213, 5.12, 5.12, 5.12, 5.12, 3.6885590546090206, 3.393660981729729, 3.7162186241472837, 3.23626238058141],
        lb_C_work = Float64[-5.12, -5.12, -3.6955815936, -3.210139775093711, -3.775091646745441, -3.3681715906305802, -5.12, -5.12, -5.12, -5.12],
        ub_C_work = Float64[3.7836375, 3.2332190011685213, 5.12, 5.12, 5.12, 5.12, 3.6885590546090206, 3.393660981729729, 3.7162186241472837, 3.23626238058141],
        ProblemData = _problem_data(),
        func_idx = _objective,
        func_F = xr -> _objective(_cell2idx(xr)),
        idx2cell = _idx2cell,
        cell2idx = _cell2idx,
        decode = xr -> _decode(_cell2idx(xr)))
end

end # module
