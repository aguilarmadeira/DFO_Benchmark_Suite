# salomon_cd_5D  DFO Benchmark Suite v2.1 core, Julia version (self-contained; Base only).
#SALOMON_CD_5D  DFO Benchmark Suite v2.1.0 core: salomon, dimension 5, pattern CD, strategy extreme.
#
# Copyright (C) 2026 J. F. A. Madeira
#
# Objective F(z) = f(T(z)), where f is the original function (embedded below, translated from the MATLAB original) and T acts on the finite
# coordinates only: uniform grid of 151 levels on the interval of the domain, piecewise-affine map with |T - Id| <=
# half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
# Global optimum of F on this domain: f* = 0 (E: exact by structure); the optimal solution is listed in
# reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
# Variables: x1:C x2:D x3:C x4:D x5:C
#
# Julia version of the MATLAB wrapper salomon_cd_5D.m: same instance, same values (level indices 1..m as in MATLAB).
# Usage:  include("salomon_cd_5D.jl");  PB = salomon_cd_5D.problem();  f = PB.func_idx(PB.x0_idx)
# Copyright (C) 2026 J. F. A. Madeira. LGPL-3.0-or-later.
module salomon_cd_5D

const NAME = "salomon_cd_5D"
const N = 5
const TYPES = Char['C', 'D', 'C', 'D', 'C']
const CARD = Float64[Inf, 151.0, Inf, 151.0, Inf]
const LB = Float64[-70.399169921875, -100.0, -66.579328, -100.0, -71.18680418004234]
const UB = Float64[100.0, 66.2599197276838, 100.0, 72.41232819658475, 100.0]
const LBW = Float64[-70.399169921875, -100.0, -6657932800.0, -100.0, -7118680418.004234]
const UBW = Float64[100.0, 66.2599197276838, 10000000000.0, 72.41232819658475, 10000000000.0]
const TB = Float64[NaN, -0.2440481633897349, NaN, -0.0008496459808355894, NaN]
const TI = Float64[NaN, 0.0, NaN, 0.0, NaN]
const CW = Float64[0.0, NaN, 0.0, NaN, 0.0]
const CX = Float64[0.0, NaN, 0.0, NaN, 0.0]
const M = 151
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
        cardinality = Float64[Inf, 151.0, Inf, 151.0, Inf],
        f_star = 0.0,
        f_star_justification = "E",
        n_minimizers_documented = 1,
        description = "DFO Benchmark Suite v2.1.0 core (Julia): CD adaptation of salomon (F = f o T).",
        x0_idx = Float64[-22.39120066165924, 94.0, 2222744334.3360004, 133.0, 6283021306.684282],
        x0_alt_idx = Vector{Float64}[Float64[21.539835333824158, 98.0, -5373273022.464, 139.0, 919768313.2148838], Float64[-64.32537138462067, 88.0, 5580983586.816, 146.0, -2900102523.860441], Float64[91.4301198720932, 90.0, -2015033769.9839993, 132.0, 8713848203.004942], Float64[4.56648051738739, 100.0, 8811289915.392, 138.0, 4893977365.929617], Float64[48.49751651287079, 89.0, 1215272558.591999, 145.0, -609583096.3028774], Float64[-37.367690205574036, 93.0, -3822086320.128, 148.0, -4429453933.378202], Float64[70.46303451061249, 96.0, 4573511811.072002, 134.0, 7184496793.487183], Float64[-16.899821162223816, 85.0, -463847067.6480007, 140.0, 3364625956.4118557], Float64[27.031214833259583, 94.0, 7931751063.552002, 147.0, -1998627037.057541]],
        suite_release = "2.1.0",
        base_function = "salomon",
        m = 151,
        pattern = "CD",
        strategy = "extreme",
        scale_factors_C = Float64[1.0, 100000000.0, 100000000.0],
        lb_orig = Float64[-70.399169921875, -100.0, -66.579328, -100.0, -71.18680418004234],
        ub_orig = Float64[100.0, 66.2599197276838, 100.0, 72.41232819658475, 100.0],
        lb_C_work = Float64[-70.399169921875, -100.0, -6657932800.0, -100.0, -7118680418.004234],
        ub_C_work = Float64[100.0, 66.2599197276838, 10000000000.0, 72.41232819658475, 10000000000.0],
        ProblemData = _problem_data(),
        func_idx = _objective,
        func_F = xr -> _objective(_cell2idx(xr)),
        idx2cell = _idx2cell,
        cell2idx = _cell2idx,
        decode = xr -> _decode(_cell2idx(xr)))
end

end # module
