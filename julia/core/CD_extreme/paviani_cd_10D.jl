# paviani_cd_10D  DFO Benchmark Suite v2.1 core, Julia version (self-contained; Base only).
#PAVIANI_CD_10D  DFO Benchmark Suite v2.1.0 core: paviani, dimension 10, pattern CD, strategy extreme.
#
# Copyright (C) 2026 J. F. A. Madeira
#
# Objective F(z) = f(T(z)), where f is the original function (embedded below, translated from the MATLAB original) and T acts on the finite
# coordinates only: uniform grid of 17 levels on the interval of the domain, piecewise-affine map with |T - Id| <=
# half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
# Global optimum of F on this domain: f* = -45.77845205382886 (A: analytic reduction, 1-D root to machine precision); the optimal solution is listed in
# reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
# Variables: x1:C x2:D x3:C x4:D x5:C x6:D x7:C x8:D x9:C x10:D
#
# Julia version of the MATLAB wrapper paviani_cd_10D.m: same instance, same values (level indices 1..m as in MATLAB).
# Usage:  include("paviani_cd_10D.jl");  PB = paviani_cd_10D.problem();  f = PB.func_idx(PB.x0_idx)
# Copyright (C) 2026 J. F. A. Madeira. LGPL-3.0-or-later.
module paviani_cd_10D

const NAME = "paviani_cd_10D"
const N = 10
const TYPES = Char['C', 'D', 'C', 'D', 'C', 'D', 'C', 'D', 'C', 'D']
const CARD = Float64[Inf, 17.0, Inf, 17.0, Inf, 17.0, Inf, 17.0, Inf, 17.0]
const LB = Float64[3.1781763354492187, 3.3917369212010366, 3.20312627328, 3.172785281132861, 3.4030945768731646, 3.197009175449039, 3.545224663737264, 3.2314999192762484, 3.5588647308292924, 3.395080028701464]
const UB = Float64[9.999, 9.999, 9.999, 9.999, 9.999, 9.999, 9.999, 9.999, 9.999, 9.999]
const LBW = Float64[3.178176335449219, 3.3917369212010366, 3.2031262732800005, 3.172785281132861, 340309457.6873164, 3.197009175449039, 354522466.3737263, 3.2314999192762484, 355886473.08292925, 3.395080028701464]
const UBW = Float64[9.999, 9.999, 9.999, 9.999, 999900000.0, 9.999, 999900000.0, 9.999, 999900000.0000001, 9.999]
const TB = Float64[NaN, 9.17309211515013, NaN, 9.145723160141607, NaN, 9.14875114693113, NaN, 9.153062489909532, NaN, 9.173510003587683]
const TI = Float64[NaN, 9.351, NaN, 9.351, NaN, 9.351, NaN, 9.351, NaN, 9.351]
const CW = Float64[9.351, NaN, 9.351, NaN, 935100000.0, NaN, 935100000.0, NaN, 935100000.0000001, NaN]
const CX = Float64[9.351, NaN, 9.351, NaN, 9.351, NaN, 9.351, NaN, 9.351, NaN]
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
        cardinality = Float64[Inf, 17.0, Inf, 17.0, Inf, 17.0, Inf, 17.0, Inf, 17.0],
        f_star = -45.778452053828865,
        f_star_justification = "A",
        n_minimizers_documented = 1,
        description = "DFO Benchmark Suite v2.1.0 core (Julia): CD adaptation of paviani (F = f o T).",
        x0_idx = Float64[4.4604112528574476, 7.0, 5.195132780056167, 1.0, 496906460.9711692, 12.0, 982428961.5362723, 10.0, 410723163.6299898, 13.0],
        x0_alt_idx = Vector{Float64}[Float64[8.776713728080988, 8.0, 9.925060893853288, 2.0, 944397264.6588876, 10.0, 801544600.9005346, 14.0, 861220338.5778972, 15.0], Float64[5.339658053365946, 7.0, 6.565180923362919, 3.0, 797215903.977545, 8.0, 582828217.1252964, 4.0, 696869060.4421608, 2.0], Float64[7.098151654382945, 7.0, 3.466262503978599, 1.0, 590567326.8592963, 5.0, 401943856.48955876, 9.0, 532517782.3064244, 5.0], Float64[3.6211302160084253, 8.0, 7.883036756448462, 2.0, 443385966.17795366, 3.0, 864203889.2253326, 14.0, 368166504.170688, 8.0], Float64[9.216337128335237, 7.0, 4.784118337064141, 3.0, 885470659.2804989, 17.0, 683319528.5895951, 3.0, 846611336.0769429, 11.0], Float64[5.779281453620196, 7.0, 9.253084899755214, 3.0, 738289298.5991561, 14.0, 502435167.9538574, 8.0, 682260057.9412065, 14.0], Float64[7.537775054637194, 8.0, 6.154166480370893, 1.0, 531640721.48090744, 12.0, 926863177.5501306, 13.0, 517908779.8054701, 16.0], Float64[4.040770734432936, 6.0, 4.099094265410765, 2.0, 384459360.7995649, 9.0, 745978816.9143931, 2.0, 996353611.711725, 3.0], Float64[8.357073209656477, 7.0, 7.524214623677646, 3.0, 831950164.4872832, 7.0, 565094456.2786554, 7.0, 804054676.6176411, 6.0]],
        suite_release = "2.1.0",
        base_function = "paviani",
        m = 17,
        pattern = "CD",
        strategy = "extreme",
        scale_factors_C = Float64[1.0, 1.0, 100000000.0, 100000000.0, 100000000.0],
        lb_orig = Float64[3.1781763354492187, 3.3917369212010366, 3.20312627328, 3.172785281132861, 3.4030945768731646, 3.197009175449039, 3.545224663737264, 3.2314999192762484, 3.5588647308292924, 3.395080028701464],
        ub_orig = Float64[9.999, 9.999, 9.999, 9.999, 9.999, 9.999, 9.999, 9.999, 9.999, 9.999],
        lb_C_work = Float64[3.178176335449219, 3.3917369212010366, 3.2031262732800005, 3.172785281132861, 340309457.6873164, 3.197009175449039, 354522466.3737263, 3.2314999192762484, 355886473.08292925, 3.395080028701464],
        ub_C_work = Float64[9.999, 9.999, 9.999, 9.999, 999900000.0, 9.999, 999900000.0, 9.999, 999900000.0000001, 9.999],
        ProblemData = _problem_data(),
        func_idx = _objective,
        func_F = xr -> _objective(_cell2idx(xr)),
        idx2cell = _idx2cell,
        cell2idx = _cell2idx,
        decode = xr -> _decode(_cell2idx(xr)))
end

end # module
