# sinusoidal_cdk_10D  DFO Benchmark Suite v2.1 core, Julia version (self-contained; Base only).
#SINUSOIDAL_CDK_10D  DFO Benchmark Suite v2.1.0 core: sinusoidal, dimension 10, pattern CDK, strategy moderate.
#
# Copyright (C) 2026 J. F. A. Madeira
#
# Objective F(z) = f(T(z)), where f is the original function (embedded below, translated from the MATLAB original) and T acts on the finite
# coordinates only: uniform grid of 17 levels on the interval of the domain, piecewise-affine map with |T - Id| <=
# half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
# Global optimum of F on this domain: f* = -3.5 (E: exact by structure); the optimal solution is listed in
# reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
# Variables: x1:C x2:K x3:D x4:C x5:D x6:K x7:C x8:K x9:D x10:C
#
# Julia version of the MATLAB wrapper sinusoidal_cdk_10D.m: same instance, same values (level indices 1..m as in MATLAB).
# Usage:  include("sinusoidal_cdk_10D.jl");  PB = sinusoidal_cdk_10D.problem();  f = PB.func_idx(PB.x0_idx)
# Copyright (C) 2026 J. F. A. Madeira. LGPL-3.0-or-later.
module sinusoidal_cdk_10D

const NAME = "sinusoidal_cdk_10D"
const N = 10
const TYPES = Char['C', 'K', 'D', 'C', 'D', 'K', 'C', 'K', 'D', 'C']
const CARD = Float64[Inf, 17.0, 17.0, Inf, 17.0, 17.0, Inf, 17.0, 17.0, Inf]
const LB = Float64[32.1532470703125, 26.235208047553726, 0.0, 0.0, 0.0, 0.0, 34.79734198584787, 26.401758734202474, 0.0, 0.0]
const UB = Float64[180.0, 180.0, 148.8125952, 147.71395251978342, 153.96399153063317, 149.13158502853543, 180.0, 180.0, 148.20241637215418, 153.82416663249825]
const LBW = Float64[32.15324707031249, 26.235208047553726, 0.0, 0.0, 0.0, 0.0, 0.34797341985847885, 26.401758734202474, 0.0, 0.0]
const UBW = Float64[180.0, 180.0, 148.8125952, 14771.395251978342, 153.96399153063317, 149.13158502853543, 1.8000000000000003, 180.0, 148.20241637215418, 153.82416663249825]
const TB = Float64[NaN, 122.33820301783265, 120.9102336, NaN, 115.47299364797487, 121.16941283568504, NaN, 122.40065952532592, 120.41446330237527, NaN]
const TI = Float64[NaN, 120.0, 120.0, NaN, 120.0, 120.0, NaN, 120.0, 120.0, NaN]
const CW = Float64[119.99999999999999, NaN, NaN, 12000.0, NaN, NaN, 1.2000000000000002, NaN, NaN, 120.0]
const CX = Float64[120.0, NaN, NaN, 120.0, NaN, NaN, 120.0, NaN, NaN, 120.0]
const M = 17
const PERM = Dict{Int,Vector{Int}}(2 => Int[3, 9, 1, 2, 7, 5, 6, 4, 8, 15, 16, 10, 17, 13, 11, 12, 14], 6 => Int[10, 4, 1, 15, 11, 16, 2, 9, 8, 12, 5, 14, 13, 3, 7, 6, 17], 8 => Int[5, 1, 2, 10, 17, 3, 9, 16, 7, 4, 8, 11, 15, 6, 13, 14, 12])

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

function _sind(v::Float64)
    # sine of an angle in degrees, exact at multiples of 90 (as MATLAB sind)
    r = rem(v, 360.0)
    n = r >= 0 ? floor(r / 90 + 0.5) : -floor(-r / 90 + 0.5)
    t = (r - 90 * n) * pi / 180
    q = mod(Int(n), 4)
    return q == 0 ? sin(t) : q == 1 ? cos(t) : q == 2 ? -sin(t) : -cos(t)
end
function _orig(x::Vector{Float64})
    A = 2.5; B = 5; z = 30
    p1 = 1.0
    for v in x; p1 *= _sind(v - z); end
    p2 = 1.0
    for v in x; p2 *= _sind(B * (v - z)); end
    return -(A * p1 + p2)
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
        f_star = -3.5,
        f_star_justification = "E",
        n_minimizers_documented = 1,
        description = "DFO Benchmark Suite v2.1.0 core (Julia): CDK adaptation of sinusoidal (F = f o T).",
        x0_idx = Float64[138.49029153585434, 1.0, 11.0, 7253.425656844009, 14.0, 1.0, 1.7657164495631155, 8.0, 14.0, 49.02518542106724],
        x0_alt_idx = Vector{Float64}[Float64[65.72196782827376, 2.0, 7.0, 7899.405041041311, 9.0, 15.0, 1.3587470620494018, 13.0, 9.0, 75.91240598584398], Float64[157.54866203069685, 3.0, 14.0, 6478.250395807244, 5.0, 13.0, 0.8666598941406628, 2.0, 5.0, 102.79962655062073], Float64[84.7803383231163, 1.0, 10.0, 7124.229780004547, 16.0, 10.0, 0.45969050662694927, 7.0, 1.0, 129.68684711539748], Float64[120.2982106089592, 2.0, 2.0, 7770.209164201853, 12.0, 8.0, 1.4997233858286618, 12.0, 13.0, 2.93280731002139], Float64[47.52988690137862, 3.0, 13.0, 8416.188548399155, 7.0, 5.0, 1.0927539983149481, 2.0, 9.0, 29.820027874798136], Float64[167.07784727811813, 2.0, 5.0, 6699.729041817748, 4.0, 3.0, 0.6857846108012347, 7.0, 5.0, 56.707248439574876], Float64[94.30952357053755, 2.0, 16.0, 7345.708426015051, 15.0, 17.0, 1.6406997096079223, 11.0, 17.0, 83.59446900435164], Float64[129.39425107240675, 2.0, 8.0, 7991.687810212353, 11.0, 14.0, 1.2337303220942086, 16.0, 13.0, 110.48168956912838], Float64[56.6259273648262, 1.0, 4.0, 6570.5331649782875, 6.0, 12.0, 0.8267609345804949, 6.0, 8.0, 137.3689101339051]],
        suite_release = "2.1.0",
        base_function = "sinusoidal",
        m = 17,
        pattern = "CDK",
        strategy = "moderate",
        scale_factors_C = Float64[1.0, 100.0, 0.01, 1.0],
        lb_orig = Float64[32.1532470703125, 26.235208047553726, 0.0, 0.0, 0.0, 0.0, 34.79734198584787, 26.401758734202474, 0.0, 0.0],
        ub_orig = Float64[180.0, 180.0, 148.8125952, 147.71395251978342, 153.96399153063317, 149.13158502853543, 180.0, 180.0, 148.20241637215418, 153.82416663249825],
        lb_C_work = Float64[32.15324707031249, 26.235208047553726, 0.0, 0.0, 0.0, 0.0, 0.34797341985847885, 26.401758734202474, 0.0, 0.0],
        ub_C_work = Float64[180.0, 180.0, 148.8125952, 14771.395251978342, 153.96399153063317, 149.13158502853543, 1.8000000000000003, 180.0, 148.20241637215418, 153.82416663249825],
        ProblemData = _problem_data(),
        func_idx = _objective,
        func_F = xr -> _objective(_cell2idx(xr)),
        idx2cell = _idx2cell,
        cell2idx = _cell2idx,
        decode = xr -> _decode(_cell2idx(xr)))
end

end # module
