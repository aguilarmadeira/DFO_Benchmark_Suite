# transistor_cdk_9D  DFO Benchmark Suite v2.1 core, Julia version (self-contained; Base only).
#TRANSISTOR_CDK_9D  DFO Benchmark Suite v2.1.0 core: transistor, dimension 9, pattern CDK, strategy extreme.
#
# Copyright (C) 2026 J. F. A. Madeira
#
# Objective F(z) = f(T(z)), where f is the original function (embedded below, translated from the MATLAB original) and T acts on the finite
# coordinates only: uniform grid of 17 levels on the interval of the domain, piecewise-affine map with |T - Id| <=
# half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
# Global optimum of F on this domain: f* = 1.60049326654736e-24 (E: exact by structure); the optimal solution is listed in
# reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
# Variables: x1:C x2:K x3:D x4:C x5:D x6:K x7:C x8:K x9:D
#
# Julia version of the MATLAB wrapper transistor_cdk_9D.m: same instance, same values (level indices 1..m as in MATLAB).
# Usage:  include("transistor_cdk_9D.jl");  PB = transistor_cdk_9D.problem();  f = PB.func_idx(PB.x0_idx)
# Copyright (C) 2026 J. F. A. Madeira. LGPL-3.0-or-later.
module transistor_cdk_9D

const NAME = "transistor_cdk_9D"
const N = 9
const TYPES = Char['C', 'K', 'D', 'C', 'D', 'K', 'C', 'K', 'D']
const CARD = Float64[Inf, 17.0, 17.0, Inf, 17.0, 17.0, Inf, 17.0, 17.0]
const LB = Float64[-10.0, -10.0, -10.0, -10.0, -7.350085376681921, -6.992661321382304, -10.0, -10.0, -10.0]
const UB = Float64[6.8211669921875, 7.03312503175329, 6.7027328, 6.180008329862558, 10.0, 10.0, 6.7876079069934505, 6.403313357018439, 7.316462562669516]
const LBW = Float64[-10.0, -10.0, -10.0, -1000000000.0, -7.350085376681921, -6.992661321382304, -1000000000.0, -10.0, -10.0]
const UBW = Float64[6.8211669921875, 7.03312503175329, 6.7027328, 618000832.9862559, 10.0, 10.0, 678760790.6993451, 6.403313357018439, 7.316462562669516]
const TB = Float64[NaN, 0.6457031448458057, 1.4831287999999994, NaN, 7.83123932791476, 7.8759173348272125, NaN, 1.2772779329501773, 1.90506801183529]
const TI = Float64[NaN, 0.4499874719815851, 1.000006482465387, NaN, 7.999971440507853, 7.999692684218543, NaN, 0.9999877234567898, 2.00005248348619]
const CW = Float64[0.8999999526168558, NaN, NaN, 200006854.1623683, NaN, NaN, 500003127.5929673, NaN, NaN]
const CX = Float64[0.8999999526168565, NaN, NaN, 2.000068541623683, NaN, NaN, 5.000031275929672, NaN, NaN]
const M = 17
const PERM = Dict{Int,Vector{Int}}(2 => Int[6, 10, 13, 11, 7, 2, 1, 15, 8, 14, 3, 16, 17, 9, 5, 12, 4], 6 => Int[12, 11, 16, 6, 10, 3, 4, 1, 8, 15, 9, 14, 5, 13, 17, 7, 2], 8 => Int[8, 3, 16, 6, 10, 11, 14, 17, 15, 4, 12, 1, 5, 2, 7, 13, 9])

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
        cardinality = Float64[Inf, 17.0, 17.0, Inf, 17.0, 17.0, Inf, 17.0, 17.0],
        f_star = 1.6004932665473596e-24,
        f_star_justification = "E",
        n_minimizers_documented = 1,
        description = "DFO Benchmark Suite v2.1.0 core (Julia): CDK adaptation of transistor (F = f o T).",
        x0_idx = Float64[-3.4210181832313538, 2.0, 9.0, 289818240.0398557, 11.0, 12.0, -145072359.38738823, 16.0, 4.0],
        x0_alt_idx = Vector{Float64}[Float64[0.9156889319419861, 3.0, 1.0, 360576293.96053743, 6.0, 10.0, -615590089.6526027, 4.0, 16.0], Float64[-7.560602247714996, 2.0, 12.0, 204908575.3350377, 2.0, 7.0, 488435172.8651228, 10.0, 12.0], Float64[3.084042489528656, 2.0, 4.0, 275666629.2557192, 13.0, 5.0, 17917442.599908352, 15.0, 8.0], Float64[-5.490810215473175, 2.0, 15.0, 346424683.17640114, 10.0, 3.0, -452600287.6653062, 4.0, 3.0], Float64[-1.1541031002998352, 1.0, 7.0, 158410425.61573243, 4.0, 16.0, -923118017.9305205, 9.0, 15.0], Float64[-9.630394279956818, 2.0, 3.0, 229168479.53641415, 17.0, 14.0, 279316181.6361382, 14.0, 11.0], Float64[6.533695876598358, 2.0, 11.0, 299926533.4570956, 12.0, 12.0, -289610485.67800975, 3.0, 7.0], Float64[-2.090437591075897, 1.0, 6.0, 370684587.3777776, 8.0, 9.0, -760128215.9432241, 8.0, 3.0], Float64[2.2462695240974426, 2.0, 14.0, 215016868.75227785, 3.0, 7.0, 442305983.6234348, 13.0, 14.0]],
        suite_release = "2.1.0",
        base_function = "transistor",
        m = 17,
        pattern = "CDK",
        strategy = "extreme",
        scale_factors_C = Float64[1.0, 100000000.0, 100000000.0],
        lb_orig = Float64[-10.0, -10.0, -10.0, -10.0, -7.350085376681921, -6.992661321382304, -10.0, -10.0, -10.0],
        ub_orig = Float64[6.8211669921875, 7.03312503175329, 6.7027328, 6.180008329862558, 10.0, 10.0, 6.7876079069934505, 6.403313357018439, 7.316462562669516],
        lb_C_work = Float64[-10.0, -10.0, -10.0, -1000000000.0, -7.350085376681921, -6.992661321382304, -1000000000.0, -10.0, -10.0],
        ub_C_work = Float64[6.8211669921875, 7.03312503175329, 6.7027328, 618000832.9862559, 10.0, 10.0, 678760790.6993451, 6.403313357018439, 7.316462562669516],
        ProblemData = _problem_data(),
        func_idx = _objective,
        func_F = xr -> _objective(_cell2idx(xr)),
        idx2cell = _idx2cell,
        cell2idx = _cell2idx,
        decode = xr -> _decode(_cell2idx(xr)))
end

end # module
