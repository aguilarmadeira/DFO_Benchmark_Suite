# griewank_c_5D  DFO Benchmark Suite v2.1 core, Julia version (self-contained; Base only).
#GRIEWANK_C_5D  DFO Benchmark Suite v2.1.0 core: griewank, dimension 5, pattern C, strategy baseline.
#
# Copyright (C) 2026 J. F. A. Madeira
#
# Objective F(z) = f(T(z)), where f is the original function (embedded below, translated from the MATLAB original) and T acts on the finite
# coordinates only: uniform grid of 0 levels on the interval of the domain, piecewise-affine map with |T - Id| <=
# half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
# Global optimum of F on this domain: f* = 0 (E: exact by structure); the optimal solution is listed in
# reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
# Variables: x1:C x2:C x3:C x4:C x5:C
#
# Julia version of the MATLAB wrapper griewank_c_5D.m: same instance, same values (level indices 1..m as in MATLAB).
# Usage:  include("griewank_c_5D.jl");  PB = griewank_c_5D.problem();  f = PB.func_idx(PB.x0_idx)
# Copyright (C) 2026 J. F. A. Madeira. LGPL-3.0-or-later.
module griewank_c_5D

const NAME = "griewank_c_5D"
const N = 5
const TYPES = Char['C', 'C', 'C', 'C', 'C']
const CARD = Float64[Inf, Inf, Inf, Inf, Inf]
const LB = Float64[-600.0, -600.0, -379.987968, -600.0, -600.0]
const UB = Float64[363.00439453125, 431.78174058832496, 600.0, 415.6168263223657, 405.6001639232293]
const LBW = Float64[-600.0, -600.0, -379.987968, -600.0, -600.0]
const UBW = Float64[363.00439453125, 431.78174058832496, 600.0, 415.6168263223657, 405.6001639232293]
const TB = Float64[NaN, NaN, NaN, NaN, NaN]
const TI = Float64[NaN, NaN, NaN, NaN, NaN]
const CW = Float64[0.0, 0.0, 0.0, 0.0, 0.0]
const CX = Float64[0.0, 0.0, 0.0, 0.0, 0.0]
const M = 0
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
    n = length(x)
    s = 0.0
    for v in x; s += v^2; end
    p = 1.0
    for i in 1:n; p *= cos(x[i] / sqrt(i)); end
    return 1 + 1 / 4000 * s - p
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
        cardinality = Float64[Inf, Inf, Inf, Inf, Inf],
        f_star = 0.0,
        f_star_justification = "E",
        n_minimizers_documented = 1,
        description = "DFO Benchmark Suite v2.1.0 core (Julia): C adaptation of griewank (F = f o T).",
        x0_idx = Float64[362.0639605522156, -279.6617275448227, 369.8204260761602, -311.93875105142405, -552.4020959224918],
        x0_alt_idx = Vector{Float64}[Float64[-133.07452940940857, -351.8439892320718, -77.05408733183998, -429.955033660312, 220.49720362180847], Float64[115.20004105567932, -322.1218814784986, 567.38600042496, -385.54030364621434, -94.55559003407933], Float64[-370.0638921260834, -258.4316505779847, 120.51148701695996, -341.1255736321167, -409.6083836899671], Float64[239.33732628822327, -341.9366199808808, -371.52087195648005, -438.8379796631315, 363.2909158543332], Float64[-251.56921076774597, -312.21451222730764, 310.55075377152, -394.4232496490339, 39.99606088632265], Float64[-3.294640302658081, -282.49240447373444, -173.95529760768, -370.31239621280946, -184.39405833621612], Float64[-488.5585734844208, -369.7716097818461, 508.1163281203201, -325.8976661987119, -499.4468519921038], Float64[301.40596890449524, -306.0813788813322, 23.610276741119947, -423.6100722297266, 273.4524475521964], Float64[-192.32187008857727, -276.3592711277591, -235.10654681088002, -379.19534221562895, -41.60034610369155]],
        suite_release = "2.1.0",
        base_function = "griewank",
        m = 0,
        pattern = "C",
        strategy = "baseline",
        scale_factors_C = Float64[1.0, 1.0, 1.0, 1.0, 1.0],
        lb_orig = Float64[-600.0, -600.0, -379.987968, -600.0, -600.0],
        ub_orig = Float64[363.00439453125, 431.78174058832496, 600.0, 415.6168263223657, 405.6001639232293],
        lb_C_work = Float64[-600.0, -600.0, -379.987968, -600.0, -600.0],
        ub_C_work = Float64[363.00439453125, 431.78174058832496, 600.0, 415.6168263223657, 405.6001639232293],
        ProblemData = _problem_data(),
        func_idx = _objective,
        func_F = xr -> _objective(_cell2idx(xr)),
        idx2cell = _idx2cell,
        cell2idx = _cell2idx,
        decode = xr -> _decode(_cell2idx(xr)))
end

end # module
