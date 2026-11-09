# paviani_c_10D  DFO Benchmark Suite v2.1 core, Julia version (self-contained; Base only).
#PAVIANI_C_10D  DFO Benchmark Suite v2.1.0 core: paviani, dimension 10, pattern C, strategy extreme.
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
const LBW = Float64[3.178176335449219, 3.391736921201037, 3.2031262732800005, 3.1727852811328616, 3.403094576873164, 319700917.5449038, 354522466.3737263, 323149991.92762476, 355886473.08292925, 339508002.87014645]
const UBW = Float64[9.999, 9.999, 9.999, 9.999, 9.999, 999900000.0, 999900000.0, 999900000.0, 999900000.0000001, 999900000.0]
const TB = Float64[NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN]
const TI = Float64[NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN]
const CW = Float64[9.351, 9.351, 9.351, 9.351, 9.351, 935100000.0, 935100000.0, 935100000.0, 935100000.0000001, 935100000.0]
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
        x0_idx = Float64[4.4604112528574476, 6.0140983077568215, 5.195132780056167, 3.320625000158722, 4.969064609711692, 796490442.7979071, 982428961.5362723, 709525998.8690771, 410723163.6299898, 822461451.5891619],
        x0_alt_idx = Vector{Float64}[Float64[8.776713728080988, 6.204430989039096, 9.925060893853288, 3.6191475097301713, 9.443972646588877, 711968663.3212383, 801544600.9005346, 893340733.5299366, 861220338.5778972, 937892633.0137618], Float64[5.339658053365946, 5.742194477353572, 6.565180923362919, 3.9176700193016205, 7.97215903977545, 627446883.8445696, 582828217.1252964, 437799869.3704153, 696869060.4421608, 393717063.4406482], Float64[7.098151654382945, 5.8600194705283135, 3.466262503978599, 3.2609204982444324, 5.905673268592963, 490911701.61302775, 401943856.48955876, 657134359.5212958, 532517782.3064244, 509148244.86524796], Float64[3.6211302160084253, 6.267875216133188, 7.883036756448462, 3.5594430078158816, 4.433859661779537, 406389922.136359, 864203889.2253326, 876468849.6721766, 368166504.170688, 624579426.2898479], Float64[9.216337128335237, 5.805638704447663, 4.784118337064141, 3.8579655173873304, 8.85470659280499, 998042378.4730403, 683319528.5895951, 420927985.5126552, 846611336.0769429, 740010607.7144477], Float64[5.779281453620196, 5.995971385729939, 9.253084899755214, 4.020020594011831, 7.382892985991562, 861507196.2414985, 502435167.9538574, 640262475.6635358, 682260057.9412065, 855441789.1390475], Float64[7.537775054637194, 6.089627149535501, 6.154166480370893, 3.3632710729546433, 5.316407214809075, 776985416.7648298, 926863177.5501306, 824077210.3243953, 517908779.8054701, 970872970.5636473], Float64[4.040770734432936, 5.627390637849977, 4.099094265410765, 3.6617935825260925, 3.844593607995649, 640450234.5332878, 745978816.9143931, 368536346.164874, 996353611.711725, 426697400.9905338], Float64[8.357073209656477, 6.035246383454852, 7.524214623677646, 3.9603160920975418, 8.319501644872833, 555928455.056619, 565094456.2786554, 587870836.3157547, 804054676.6176411, 542128582.4151337]],
        suite_release = "2.1.0",
        base_function = "paviani",
        m = 0,
        pattern = "C",
        strategy = "extreme",
        scale_factors_C = Float64[1.0, 1.0, 1.0, 1.0, 1.0, 100000000.0, 100000000.0, 100000000.0, 100000000.0, 100000000.0],
        lb_orig = Float64[3.1781763354492187, 3.3917369212010366, 3.20312627328, 3.172785281132861, 3.4030945768731646, 3.197009175449039, 3.545224663737264, 3.2314999192762484, 3.5588647308292924, 3.395080028701464],
        ub_orig = Float64[9.999, 9.999, 9.999, 9.999, 9.999, 9.999, 9.999, 9.999, 9.999, 9.999],
        lb_C_work = Float64[3.178176335449219, 3.391736921201037, 3.2031262732800005, 3.1727852811328616, 3.403094576873164, 319700917.5449038, 354522466.3737263, 323149991.92762476, 355886473.08292925, 339508002.87014645],
        ub_C_work = Float64[9.999, 9.999, 9.999, 9.999, 9.999, 999900000.0, 999900000.0, 999900000.0, 999900000.0000001, 999900000.0],
        ProblemData = _problem_data(),
        func_idx = _objective,
        func_F = xr -> _objective(_cell2idx(xr)),
        idx2cell = _idx2cell,
        cell2idx = _cell2idx,
        decode = xr -> _decode(_cell2idx(xr)))
end

end # module
