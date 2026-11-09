# neumaier2_cdk_4D  DFO Benchmark Suite v2.1 core, Julia version (self-contained; Base only).
#NEUMAIER2_CDK_4D  DFO Benchmark Suite v2.1.0 core: neumaier2, dimension 4, pattern CDK, strategy baseline.
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
# Julia version of the MATLAB wrapper neumaier2_cdk_4D.m: same instance, same values (level indices 1..m as in MATLAB).
# Usage:  include("neumaier2_cdk_4D.jl");  PB = neumaier2_cdk_4D.problem();  f = PB.func_idx(PB.x0_idx)
# Copyright (C) 2026 J. F. A. Madeira. LGPL-3.0-or-later.
module neumaier2_cdk_4D

const NAME = "neumaier2_cdk_4D"
const N = 4
const TYPES = Char['C', 'K', 'D', 'C']
const CARD = Float64[Inf, 151.0, 151.0, Inf]
const LB = Float64[0.0, 0.0, 0.75801344, 0.0]
const UB = Float64[3.3762646484375, 3.242235431590713, 4.0, 3.2196751353602666]
const LBW = Float64[0.0, 0.0, 0.75801344, 0.0]
const UBW = Float64[3.3762646484375, 3.242235431590713, 4.0, 3.2196751353602666]
const TB = Float64[NaN, 2.010185967586242, 1.9899683328000002, NaN]
const TI = Float64[NaN, 2.0, 2.0, NaN]
const CW = Float64[1.0, NaN, NaN, 3.0]
const CX = Float64[1.0, NaN, NaN, 3.0]
const M = 151
const PERM = Dict{Int,Vector{Int}}(2 => Int[91, 124, 59, 110, 103, 13, 106, 38, 151, 65, 3, 105, 131, 41, 85, 140, 16, 18, 96, 121, 28, 45, 118, 119, 56, 40, 44, 114, 9, 87, 43, 83, 148, 48, 104, 123, 7, 17, 19, 60, 29, 90, 70, 2, 8, 94, 149, 75, 92, 81, 112, 84, 126, 79, 20, 31, 39, 73, 5, 23, 50, 77, 113, 14, 35, 117, 147, 80, 21, 143, 61, 25, 76, 150, 30, 55, 74, 24, 53, 98, 47, 138, 127, 141, 108, 22, 51, 145, 10, 36, 97, 57, 135, 66, 89, 109, 37, 15, 52, 82, 72, 34, 68, 33, 46, 86, 139, 116, 120, 107, 134, 11, 93, 137, 115, 12, 125, 102, 88, 128, 63, 78, 130, 95, 99, 49, 100, 62, 32, 42, 144, 133, 26, 54, 146, 136, 6, 69, 142, 122, 58, 132, 111, 129, 67, 4, 64, 101, 27, 71, 1])

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

const _b = [8, 18, 44, 114]
function _orig(x::Vector{Float64})
    s = 0.0
    for k in 1:4
        t = 0.0
        for v in x; t += v^k; end
        s += (_b[k] - t)^2
    end
    return s
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
        cardinality = Float64[Inf, 151.0, 151.0, Inf],
        f_star = 0.0,
        f_star_justification = "E",
        n_minimizers_documented = 1,
        description = "DFO Benchmark Suite v2.1.0 core (Julia): CDK adaptation of neumaier2 (F = f o T).",
        x0_idx = Float64[1.0567312693595887, 138.0, 129.0, 2.3037908715322524],
        x0_alt_idx = Vector{Float64}[Float64[1.9271744990348816, 143.0, 89.0, 2.4445929911544213], Float64[0.22585364103317263, 147.0, 15.0, 2.5853951107765902], Float64[2.9954457354545596, 135.0, 120.0, 2.7261972303987596], Float64[1.274342076778412, 144.0, 44.0, 2.416432567229988], Float64[2.144785306453705, 148.0, 149.0, 2.5572346868521567], Float64[0.44346444845199584, 138.0, 75.0, 2.6336701232184767], Float64[2.5800069212913512, 140.0, 6.0, 2.3239054600497044], Float64[0.8490118622779846, 149.0, 105.0, 2.464707579671874], Float64[1.7194550919532776, 139.0, 36.0, 2.6055096992940427]],
        suite_release = "2.1.0",
        base_function = "neumaier2",
        m = 151,
        pattern = "CDK",
        strategy = "baseline",
        scale_factors_C = Float64[1.0, 1.0],
        lb_orig = Float64[0.0, 0.0, 0.75801344, 0.0],
        ub_orig = Float64[3.3762646484375, 3.242235431590713, 4.0, 3.2196751353602666],
        lb_C_work = Float64[0.0, 0.0, 0.75801344, 0.0],
        ub_C_work = Float64[3.3762646484375, 3.242235431590713, 4.0, 3.2196751353602666],
        ProblemData = _problem_data(),
        func_idx = _objective,
        func_F = xr -> _objective(_cell2idx(xr)),
        idx2cell = _idx2cell,
        cell2idx = _cell2idx,
        decode = xr -> _decode(_cell2idx(xr)))
end

end # module
