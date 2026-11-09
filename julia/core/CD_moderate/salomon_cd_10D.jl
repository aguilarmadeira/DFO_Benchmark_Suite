# salomon_cd_10D  DFO Benchmark Suite v2.1 core, Julia version (self-contained; Base only).
#SALOMON_CD_10D  DFO Benchmark Suite v2.1.0 core: salomon, dimension 10, pattern CD, strategy moderate.
#
# Copyright (C) 2026 J. F. A. Madeira
#
# Objective F(z) = f(T(z)), where f is the original function (embedded below, translated from the MATLAB original) and T acts on the finite
# coordinates only: uniform grid of 17 levels on the interval of the domain, piecewise-affine map with |T - Id| <=
# half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
# Global optimum of F on this domain: f* = 0 (E: exact by structure); the optimal solution is listed in
# reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
# Variables: x1:C x2:D x3:C x4:D x5:C x6:D x7:C x8:D x9:C x10:D
#
# Julia version of the MATLAB wrapper salomon_cd_10D.m: same instance, same values (level indices 1..m as in MATLAB).
# Usage:  include("salomon_cd_10D.jl");  PB = salomon_cd_10D.problem();  f = PB.func_idx(PB.x0_idx)
# Copyright (C) 2026 J. F. A. Madeira. LGPL-3.0-or-later.
module salomon_cd_10D

const NAME = "salomon_cd_10D"
const N = 10
const TYPES = Char['C', 'D', 'C', 'D', 'C', 'D', 'C', 'D', 'C', 'D']
const CARD = Float64[Inf, 17.0, Inf, 17.0, Inf, 17.0, Inf, 17.0, Inf, 17.0]
const LB = Float64[-100.0, -100.0, -63.779328, -70.41232819658475, -69.91407690731508, -62.55383214873429, -100.0, -100.0, -100.0, -100.0]
const UB = Float64[63.399169921875, 61.59325306101712, 100.0, 100.0, 100.0, 100.0, 69.5715808000383, 64.0719147336193, 70.75630804635489, 61.75997375866169]
const LBW = Float64[-100.0, -100.0, -6377.9328, -70.41232819658475, -0.6991407690731507, -62.55383214873429, -100.0, -100.0, -10000.0, -100.0]
const UBW = Float64[63.399169921875, 61.59325306101712, 10000.0, 100.0, 0.9999999999999999, 100.0, 69.5715808000383, 64.0719147336193, 7075.63080463549, 61.75997375866169]
const TB = Float64[NaN, 0.9957831631356981, NaN, 4.143065389421082, NaN, -1.5961450929589347, NaN, 2.544946708512043, NaN, 1.0999835991635507]
const TI = Float64[NaN, 0.0, NaN, 0.0, NaN, 0.0, NaN, 0.0, NaN, 0.0]
const CW = Float64[0.0, NaN, 9.094947017729282e-13, NaN, 0.0, NaN, 0.0, NaN, 0.0, NaN]
const CX = Float64[0.0, NaN, 0.0, NaN, 0.0, NaN, 0.0, NaN, 0.0, NaN]
const M = 17
const PERM = Dict{Int,Vector{Int}}()

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
    for v in x; s += v^2; end
    r = sqrt(s)
    return 1 - cos(2 * pi * r) + 0.1 * r
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
        cardinality = Float64[Inf, 17.0, Inf, 17.0, Inf, 17.0, Inf, 17.0, Inf, 17.0],
        f_star = 0.0,
        f_star_justification = "E",
        n_minimizers_documented = 1,
        description = "DFO Benchmark Suite v2.1.0 core (Julia): CD adaptation of salomon (F = f o T).",
        x0_idx = Float64[27.73538625240326, 16.0, 5629.057294336001, 2.0, 0.7855329457518486, 5.0, -34.18013340409669, 16.0, -3349.107061463994, 17.0],
        x0_alt_idx = Vector{Float64}[Float64[-52.68764269351959, 17.0, -1839.2800624639995, 2.0, 0.2531950789573303, 3.0, -81.70711625808664, 4.0, -8447.797512128967, 3.0], Float64[-12.476128220558167, 16.0, 8930.848546816, 3.0, -0.12595203480279415, 17.0, 29.81044481761532, 10.0, 4237.879059178682, 6.0], Float64[-92.89915716648102, 16.0, 1462.511190015999, 2.0, -0.6582899015973124, 14.0, -17.716538036374658, 15.0, -119.79609890409847, 9.0], Float64[53.58564555644989, 17.0, -3615.9582126079995, 2.0, 0.6476612680208945, 12.0, -65.24352089036464, 4.0, -4477.471256986877, 12.0], Float64[-26.83738338947296, 16.0, 4638.519918592002, 3.0, 0.10139696913234031, 8.0, 56.214324180943095, 9.0, -8835.14641506966, 15.0], Float64[13.374131083488464, 16.0, -314.16696012800094, 1.0, -0.27775014462778425, 6.0, 8.687341326953117, 14.0, 3850.53015623799, 1.0], Float64[-67.04889786243439, 16.0, 7940.311171072001, 2.0, 0.8750102719560288, 4.0, -48.779925522642614, 3.0, -507.1450018447904, 4.0], Float64[33.001179814338684, 15.0, 2987.6242923519994, 2.0, 0.4958631581959042, 1.0, -96.30690837663259, 8.0, -4864.82015992757, 7.0], Float64[-47.42184913158417, 16.0, -4480.713064448, 3.0, -0.036474708598614125, 15.0, 25.150936694675124, 13.0, -9963.510610592544, 10.0]],
        suite_release = "2.1.0",
        base_function = "salomon",
        m = 17,
        pattern = "CD",
        strategy = "moderate",
        scale_factors_C = Float64[1.0, 100.0, 0.01, 1.0, 100.0],
        lb_orig = Float64[-100.0, -100.0, -63.779328, -70.41232819658475, -69.91407690731508, -62.55383214873429, -100.0, -100.0, -100.0, -100.0],
        ub_orig = Float64[63.399169921875, 61.59325306101712, 100.0, 100.0, 100.0, 100.0, 69.5715808000383, 64.0719147336193, 70.75630804635489, 61.75997375866169],
        lb_C_work = Float64[-100.0, -100.0, -6377.9328, -70.41232819658475, -0.6991407690731507, -62.55383214873429, -100.0, -100.0, -10000.0, -100.0],
        ub_C_work = Float64[63.399169921875, 61.59325306101712, 10000.0, 100.0, 0.9999999999999999, 100.0, 69.5715808000383, 64.0719147336193, 7075.63080463549, 61.75997375866169],
        ProblemData = _problem_data(),
        func_idx = _objective,
        func_F = xr -> _objective(_cell2idx(xr)),
        idx2cell = _idx2cell,
        cell2idx = _cell2idx,
        decode = xr -> _decode(_cell2idx(xr)))
end

end # module
