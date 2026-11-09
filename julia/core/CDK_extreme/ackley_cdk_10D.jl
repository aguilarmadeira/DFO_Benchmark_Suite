# ackley_cdk_10D  DFO Benchmark Suite v2.1 core, Julia version (self-contained; Base only).
#ACKLEY_CDK_10D  DFO Benchmark Suite v2.1.0 core: ackley, dimension 10, pattern CDK, strategy extreme.
#
# Copyright (C) 2026 J. F. A. Madeira
#
# Objective F(z) = f(T(z)), where f is the original function (embedded below, translated from the MATLAB original) and T acts on the finite
# coordinates only: uniform grid of 17 levels on the interval of the domain, piecewise-affine map with |T - Id| <=
# half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
# Global optimum of F on this domain: f* = 4.440892098500626e-16 (E: exact by structure); the optimal solution is listed in
# reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
# Variables: x1:C x2:K x3:D x4:C x5:D x6:K x7:C x8:K x9:D x10:C
#
# Julia version of the MATLAB wrapper ackley_cdk_10D.m: same instance, same values (level indices 1..m as in MATLAB).
# Usage:  include("ackley_cdk_10D.jl");  PB = ackley_cdk_10D.problem();  f = PB.func_idx(PB.x0_idx)
# Copyright (C) 2026 J. F. A. Madeira. LGPL-3.0-or-later.
module ackley_cdk_10D

const NAME = "ackley_cdk_10D"
const N = 10
const TYPES = Char['C', 'K', 'D', 'C', 'D', 'K', 'C', 'K', 'D', 'C']
const CARD = Float64[Inf, 17.0, 17.0, Inf, 17.0, 17.0, Inf, 17.0, 17.0, Inf]
const LB = Float64[-19.659594726562503, -30.0, -21.3513984, -21.050229071220325, -19.970767024110373, -30.0, -18.69154344416374, -20.12905364446252, -20.512336648311006, -20.415744802984953]
const UB = Float64[30.0, 18.374272214601433, 30.0, 30.0, 30.0, 19.809936626868804, 30.0, 30.0, 30.0, 30.0]
const LBW = Float64[-19.659594726562503, -30.0, -21.3513984, -21.05022907122033, -19.970767024110373, -30.0, -1869154344.4163742, -20.12905364446252, -20.512336648311006, -2041574480.2984953]
const UBW = Float64[30.0, 18.374272214601433, 30.0, 30.0, 30.0, 19.809936626868804, 3000000000.0, 30.0, 30.0, 3000000000.0]
const TB = Float64[NaN, 0.23392013412589563, 1.1148384, NaN, -1.2317293900689847, 1.1312103917930045, NaN, -1.3306585277890761, -1.5702104051943806, NaN]
const TI = Float64[NaN, 0.0, 0.0, NaN, 0.0, 0.0, NaN, 0.0, 0.0, NaN]
const CW = Float64[0.0, NaN, NaN, -3.552713678800501e-15, NaN, NaN, 0.0, NaN, NaN, 0.0]
const CX = Float64[0.0, NaN, NaN, 0.0, NaN, NaN, 0.0, NaN, NaN, 0.0]
const M = 17
const PERM = Dict{Int,Vector{Int}}(2 => Int[4, 11, 9, 7, 8, 3, 1, 12, 13, 5, 2, 16, 14, 6, 15, 10, 17], 6 => Int[2, 13, 12, 10, 17, 7, 16, 1, 5, 8, 14, 3, 6, 9, 11, 4, 15], 8 => Int[4, 8, 13, 17, 3, 14, 15, 7, 5, 10, 1, 11, 16, 12, 2, 6, 9])

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
    n = length(x)
    s1 = 0.0
    for v in x; s1 += v^2; end
    s2 = 0.0
    for v in x; s2 += cos(2 * pi * v); end
    return -20 * exp(-0.02 * sqrt(1 / n * s1)) - exp(1 / n * s2) + 20 + exp(1)
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
        f_star = 4.440892098500626e-16,
        f_star_justification = "E",
        n_minimizers_documented = 1,
        description = "DFO Benchmark Suite v2.1.0 core (Julia): CDK adaptation of ackley (F = f o T).",
        x0_idx = Float64[9.777293944358824, 12.0, 4.0, -16.98917378067344, 1.0, 2.0, 2548069533.6751747, 13.0, 10.0, 596316861.8266063],
        x0_alt_idx = Vector{Float64}[Float64[-15.755691039562228, 12.0, 15.0, -15.777235814384577, 14.0, 15.0, 1183358454.4442883, 2.0, 5.0, 1477543316.76582], Float64[24.98069525957108, 13.0, 8.0, -20.688773888292072, 9.0, 13.0, -181352624.78659916, 6.0, 1.0, 2358769771.7050343], Float64[-0.04308489561081075, 12.0, 3.0, -18.456256581970482, 5.0, 10.0, -1831493472.2226388, 11.0, 13.0, -1795583515.8655455], Float64[12.759779369831087, 13.0, 11.0, -16.223739275648896, 16.0, 8.0, 1656101508.0340724, 17.0, 9.0, -914357060.926332], Float64[-12.554974973201755, 13.0, 6.0, -13.991221969327306, 12.0, 6.0, 291390428.803185, 6.0, 5.0, -33130605.987118244], Float64[18.87023731470108, 12.0, 14.0, -19.92333938326753, 7.0, 2.0, -1073320650.4277016, 11.0, 17.0, 848095848.9520955], Float64[-6.153542840480807, 13.0, 10.0, -17.69082207694594, 4.0, 16.0, 2414274329.829009, 16.0, 13.0, 1729322303.8913102], Float64[6.649321424961087, 13.0, 2.0, -15.45830477062435, 15.0, 14.0, 764133482.3929701, 5.0, 9.0, 2610548758.8305235], Float64[-18.810920011997226, 12.0, 13.0, -20.369842844531846, 11.0, 11.0, -600577596.8379173, 10.0, 4.0, -1543804528.740056]],
        suite_release = "2.1.0",
        base_function = "ackley",
        m = 17,
        pattern = "CDK",
        strategy = "extreme",
        scale_factors_C = Float64[1.0, 1.0, 100000000.0, 100000000.0],
        lb_orig = Float64[-19.659594726562503, -30.0, -21.3513984, -21.050229071220325, -19.970767024110373, -30.0, -18.69154344416374, -20.12905364446252, -20.512336648311006, -20.415744802984953],
        ub_orig = Float64[30.0, 18.374272214601433, 30.0, 30.0, 30.0, 19.809936626868804, 30.0, 30.0, 30.0, 30.0],
        lb_C_work = Float64[-19.659594726562503, -30.0, -21.3513984, -21.05022907122033, -19.970767024110373, -30.0, -1869154344.4163742, -20.12905364446252, -20.512336648311006, -2041574480.2984953],
        ub_C_work = Float64[30.0, 18.374272214601433, 30.0, 30.0, 30.0, 19.809936626868804, 3000000000.0, 30.0, 30.0, 3000000000.0],
        ProblemData = _problem_data(),
        func_idx = _objective,
        func_F = xr -> _objective(_cell2idx(xr)),
        idx2cell = _idx2cell,
        cell2idx = _cell2idx,
        decode = xr -> _decode(_cell2idx(xr)))
end

end # module
