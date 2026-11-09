# tenn_local_minima_cdk_6D  DFO Benchmark Suite v2.1 core, Julia version (self-contained; Base only).
#TENN_LOCAL_MINIMA_CDK_6D  DFO Benchmark Suite v2.1.0 core: tenn_local_minima, dimension 6, pattern CDK, strategy extreme.
#
# Copyright (C) 2026 J. F. A. Madeira
#
# Objective F(z) = f(T(z)), where f is the original function (embedded below, translated from the MATLAB original) and T acts on the finite
# coordinates only: uniform grid of 25 levels on the interval of the domain, piecewise-affine map with |T - Id| <=
# half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
# Global optimum of F on this domain: f* = 7.852723858933194e-32 (E: exact by structure); the optimal solution is listed in
# reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
# Variables: x1:C x2:K x3:D x4:C x5:D x6:K
#
# Julia version of the MATLAB wrapper tenn_local_minima_cdk_6D.m: same instance, same values (level indices 1..m as in MATLAB).
# Usage:  include("tenn_local_minima_cdk_6D.jl");  PB = tenn_local_minima_cdk_6D.problem();  f = PB.func_idx(PB.x0_idx)
# Copyright (C) 2026 J. F. A. Madeira. LGPL-3.0-or-later.
module tenn_local_minima_cdk_6D

const NAME = "tenn_local_minima_cdk_6D"
const N = 6
const TYPES = Char['C', 'K', 'D', 'C', 'D', 'K']
const CARD = Float64[Inf, 25.0, 25.0, Inf, 25.0, 25.0]
const LB = Float64[-10.0, -10.0, -10.0, -10.0, -10.0, -10.0]
const UB = Float64[6.2961669921875, 7.188680587308845, 6.1987328, 6.780008329862557, 6.343473806433986, 7.315738244459228]
const LBW = Float64[-10.0, -10.0, -10.0, -1000000000.0, -10.0, -10.0]
const UBW = Float64[6.2961669921875005, 7.188680587308845, 6.1987328, 678000832.9862556, 6.343473806433986, 7.315738244459228]
const TB = Float64[NaN, 0.7429253670680271, 0.7991552000000013, NaN, 0.895649204289322, 0.8223364027870161]
const TI = Float64[NaN, 1.0, 1.0, NaN, 1.0, 1.0]
const CW = Float64[1.0, NaN, NaN, 100000000.0, NaN, NaN]
const CX = Float64[1.0, NaN, NaN, 1.0, NaN, NaN]
const M = 25
const PERM = Dict{Int,Vector{Int}}(2 => Int[17, 14, 15, 23, 1, 20, 3, 25, 13, 8, 7, 18, 21, 6, 10, 9, 4, 11, 22, 12, 19, 2, 5, 24, 16], 6 => Int[23, 11, 5, 19, 4, 10, 24, 14, 8, 6, 22, 12, 3, 9, 1, 25, 2, 16, 13, 20, 18, 7, 21, 15, 17])

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
    for i in 1:d-1; s += ((x[i] - 1)^2) * (1 + 10 * sin(pi * x[i+1])^2); end
    f = pi / d * s
    return f + pi / d * (10 * sin(pi * x[1])^2 + (x[d] - 1)^2)
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
        cardinality = Float64[Inf, 25.0, 25.0, Inf, 25.0, 25.0],
        f_star = 7.852723858933194e-32,
        f_star_justification = "E",
        n_minimizers_documented = 1,
        description = "DFO Benchmark Suite v2.1.0 core (Julia): CDK adaptation of tenn_local_minima (F = f o T).",
        x0_idx = Float64[2.4847099661827095, 6.0, 22.0, -381494903.2932794, 9.0, 13.0],
        x0_alt_idx = Vector{Float64}[Float64[-5.536059725284576, 8.0, 15.0, -308112942.66705835, 2.0, 10.0], Float64[-1.5256748795509338, 9.0, 3.0, -469553256.0447447, 20.0, 5.0], Float64[-9.546444571018219, 7.0, 20.0, -396171295.41852355, 13.0, 2.0], Float64[5.826697337627412, 7.0, 8.0, -322789334.7923025, 7.0, 23.0], Float64[-2.194072353839873, 9.0, 25.0, -517775687.31340426, 23.0, 18.0], Float64[1.81631249189377, 7.0, 13.0, -444393726.68718326, 18.0, 15.0], Float64[-6.2044571995735165, 8.0, 1.0, -371011766.0609622, 11.0, 12.0], Float64[3.7737622380256663, 8.0, 18.0, -297629805.4347409, 5.0, 7.0], Float64[-4.247007453441619, 7.0, 6.0, -459070118.8124274, 21.0, 4.0]],
        suite_release = "2.1.0",
        base_function = "tenn_local_minima",
        m = 25,
        pattern = "CDK",
        strategy = "extreme",
        scale_factors_C = Float64[1.0, 100000000.0],
        lb_orig = Float64[-10.0, -10.0, -10.0, -10.0, -10.0, -10.0],
        ub_orig = Float64[6.2961669921875, 7.188680587308845, 6.1987328, 6.780008329862557, 6.343473806433986, 7.315738244459228],
        lb_C_work = Float64[-10.0, -10.0, -10.0, -1000000000.0, -10.0, -10.0],
        ub_C_work = Float64[6.2961669921875005, 7.188680587308845, 6.1987328, 678000832.9862556, 6.343473806433986, 7.315738244459228],
        ProblemData = _problem_data(),
        func_idx = _objective,
        func_F = xr -> _objective(_cell2idx(xr)),
        idx2cell = _idx2cell,
        cell2idx = _cell2idx,
        decode = xr -> _decode(_cell2idx(xr)))
end

end # module
