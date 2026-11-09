# paviani_c_10D  DFO Benchmark Suite v2.1 core, Julia version (self-contained; Base only).
#PAVIANI_C_10D  DFO Benchmark Suite v2.1.0 core: paviani, dimension 10, pattern C, strategy moderate.
#
# Copyright (C) 2026 J. F. A. Madeira
#
# Objective F(z) = f(T(z)), where f is the original function (embedded below, translated from the MATLAB original) and T acts on the finite
# coordinates only: uniform grid of 0 levels on the interval of the domain, piecewise-affine map with |T - Id| <=
# half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
# Global optimum of F on this domain: f* = -45.77845205382886 (A: analytic reduction, 1-D root to machine precision); the optimal solution is listed in
# reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
# Variables: x1:C x2:C x3:C x4:C x5:C x6:C x7:C x8:C x9:C x10:C
#
# Julia version of the MATLAB wrapper paviani_c_10D.m: same instance, same values (level indices 1..m as in MATLAB).
# Usage:  include("paviani_c_10D.jl");  PB = paviani_c_10D.problem();  f = PB.func_idx(PB.x0_idx)
# Copyright (C) 2026 J. F. A. Madeira. LGPL-3.0-or-later.
module paviani_c_10D

const NAME = "paviani_c_10D"
const N = 10
const TYPES = Char['C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C']
const CARD = Float64[Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf]
const LB = Float64[3.1781763354492187, 3.3917369212010366, 3.20312627328, 3.172785281132861, 3.4030945768731646, 3.197009175449039, 3.545224663737264, 3.2314999192762484, 3.5588647308292924, 3.395080028701464]
const UB = Float64[9.999, 9.999, 9.999, 9.999, 9.999, 9.999, 9.999, 9.999, 9.999, 9.999]
const LBW = Float64[3.178176335449219, 339.1736921201037, 0.032031262732800006, 3.1727852811328616, 340.3094576873164, 0.031970091754490386, 3.5452246637372635, 323.14999192762485, 0.03558864730829291, 3.3950800287014644]
const UBW = Float64[9.999, 999.9000000000001, 0.09999000000000001, 9.999, 999.9000000000001, 0.09999, 9.999, 999.9000000000001, 0.09999, 9.999]
const TB = Float64[NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN]
const TI = Float64[NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN]
const CW = Float64[9.351, 935.1000000000001, 0.09351000000000001, 9.351, 935.1000000000001, 0.09351, 9.351, 935.1000000000001, 0.09350999999999998, 9.351]
const CX = Float64[9.351, 9.351, 9.351, 9.351, 9.351, 9.351, 9.351, 9.351, 9.351, 9.351]
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
    s = 0.0
    for v in x; s += (log(v - 2))^2 + (log(10 - v))^2; end
    p = 1.0
    for v in x; p *= v; end
    return s - p^0.2
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
        cardinality = Float64[Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf],
        f_star = -45.778452053828865,
        f_star_justification = "A",
        n_minimizers_documented = 1,
        description = "DFO Benchmark Suite v2.1.0 core (Julia): C adaptation of paviani (F = f o T).",
        x0_idx = Float64[4.4604112528574476, 601.4098307756822, 0.05195132780056168, 3.320625000158722, 496.90646097116917, 0.07964904427979069, 9.824289615362723, 709.5259988690772, 0.04107231636299896, 8.22461451589162],
        x0_alt_idx = Vector{Float64}[Float64[8.776713728080988, 620.4430989039097, 0.09925060893853288, 3.6191475097301713, 944.3972646588877, 0.07119686633212383, 8.015446009005348, 893.3407335299366, 0.08612203385778969, 9.378926330137618], Float64[5.339658053365946, 574.2194477353572, 0.0656518092336292, 3.9176700193016205, 797.215903977545, 0.06274468838445696, 5.8282821712529636, 437.79986937041537, 0.06968690604421607, 3.937170634406482], Float64[7.098151654382945, 586.0019470528314, 0.03466262503978599, 3.2609204982444324, 590.5673268592964, 0.04909117016130277, 4.019438564895588, 657.1343595212959, 0.05325177823064242, 5.09148244865248], Float64[3.6211302160084253, 626.7875216133189, 0.07883036756448461, 3.5594430078158816, 443.3859661779537, 0.0406389922136359, 8.642038892253328, 876.4688496721767, 0.03681665041706879, 6.245794262898479], Float64[9.216337128335237, 580.5638704447664, 0.047841183370641416, 3.8579655173873304, 885.4706592804989, 0.09980423784730402, 6.833195285895951, 420.9279855126553, 0.08466113360769426, 7.4001060771444775], Float64[5.779281453620196, 599.5971385729938, 0.09253084899755214, 4.020020594011831, 738.2892985991564, 0.08615071962414983, 5.024351679538575, 640.262475663536, 0.06822600579412064, 8.554417891390475], Float64[7.537775054637194, 608.9627149535502, 0.061541664803708936, 3.3632710729546433, 531.6407214809075, 0.07769854167648296, 9.268631775501309, 824.0772103243953, 0.05179087798054699, 9.708729705636475], Float64[4.040770734432936, 562.7390637849977, 0.04099094265410766, 3.6617935825260925, 384.4593607995649, 0.06404502345332877, 7.459788169143931, 368.5363461648741, 0.09963536117117247, 4.266974009905338], Float64[8.357073209656477, 603.5246383454852, 0.07524214623677647, 3.9603160920975418, 831.9501644872835, 0.0555928455056619, 5.650944562786554, 587.8708363157548, 0.0804054676617641, 5.421285824151337]],
        suite_release = "2.1.0",
        base_function = "paviani",
        m = 0,
        pattern = "C",
        strategy = "moderate",
        scale_factors_C = Float64[1.0, 100.0, 0.01, 1.0, 100.0, 0.01, 1.0, 100.0, 0.01, 1.0],
        lb_orig = Float64[3.1781763354492187, 3.3917369212010366, 3.20312627328, 3.172785281132861, 3.4030945768731646, 3.197009175449039, 3.545224663737264, 3.2314999192762484, 3.5588647308292924, 3.395080028701464],
        ub_orig = Float64[9.999, 9.999, 9.999, 9.999, 9.999, 9.999, 9.999, 9.999, 9.999, 9.999],
        lb_C_work = Float64[3.178176335449219, 339.1736921201037, 0.032031262732800006, 3.1727852811328616, 340.3094576873164, 0.031970091754490386, 3.5452246637372635, 323.14999192762485, 0.03558864730829291, 3.3950800287014644],
        ub_C_work = Float64[9.999, 999.9000000000001, 0.09999000000000001, 9.999, 999.9000000000001, 0.09999, 9.999, 999.9000000000001, 0.09999, 9.999],
        ProblemData = _problem_data(),
        func_idx = _objective,
        func_F = xr -> _objective(_cell2idx(xr)),
        idx2cell = _idx2cell,
        cell2idx = _cell2idx,
        decode = xr -> _decode(_cell2idx(xr)))
end

end # module
