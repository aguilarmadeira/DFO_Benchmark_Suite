# paviani_cdk_10D  DFO Benchmark Suite v2.1 core, Julia version (self-contained; Base only).
#PAVIANI_CDK_10D  DFO Benchmark Suite v2.1.0 core: paviani, dimension 10, pattern CDK, strategy moderate.
#
# Copyright (C) 2026 J. F. A. Madeira
#
# Objective F(z) = f(T(z)), where f is the original function (embedded below, translated from the MATLAB original) and T acts on the finite
# coordinates only: uniform grid of 17 levels on the interval of the domain, piecewise-affine map with |T - Id| <=
# half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
# Global optimum of F on this domain: f* = -45.77845205382886 (A: analytic reduction, 1-D root to machine precision); the optimal solution is listed in
# reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
# Variables: x1:C x2:K x3:D x4:C x5:D x6:K x7:C x8:K x9:D x10:C
#
# Julia version of the MATLAB wrapper paviani_cdk_10D.m: same instance, same values (level indices 1..m as in MATLAB).
# Usage:  include("paviani_cdk_10D.jl");  PB = paviani_cdk_10D.problem();  f = PB.func_idx(PB.x0_idx)
# Copyright (C) 2026 J. F. A. Madeira. LGPL-3.0-or-later.
module paviani_cdk_10D

const NAME = "paviani_cdk_10D"
const N = 10
const TYPES = Char['C', 'K', 'D', 'C', 'D', 'K', 'C', 'K', 'D', 'C']
const CARD = Float64[Inf, 17.0, 17.0, Inf, 17.0, 17.0, Inf, 17.0, 17.0, Inf]
const LB = Float64[3.1781763354492187, 3.3917369212010366, 3.20312627328, 3.172785281132861, 3.4030945768731646, 3.197009175449039, 3.545224663737264, 3.2314999192762484, 3.5588647308292924, 3.395080028701464]
const UB = Float64[9.999, 9.999, 9.999, 9.999, 9.999, 9.999, 9.999, 9.999, 9.999, 9.999]
const LBW = Float64[3.178176335449219, 3.3917369212010366, 3.20312627328, 317.2785281132862, 3.4030945768731646, 3.197009175449039, 0.035452246637372634, 3.2314999192762484, 3.5588647308292924, 3.3950800287014644]
const UBW = Float64[9.999, 9.999, 9.999, 999.9000000000001, 9.999, 9.999, 0.09999000000000001, 9.999, 9.999, 9.999]
const TB = Float64[NaN, 9.17309211515013, 9.14951578416, NaN, 9.174511822109146, 9.14875114693113, NaN, 9.153062489909532, 9.193983091353662, NaN]
const TI = Float64[NaN, 9.351, 9.351, NaN, 9.351, 9.351, NaN, 9.351, 9.351, NaN]
const CW = Float64[9.351, NaN, NaN, 935.1000000000001, NaN, NaN, 0.09351000000000001, NaN, NaN, 9.351]
const CX = Float64[9.351, NaN, NaN, 9.351, NaN, NaN, 9.351, NaN, NaN, 9.351]
const M = 17
const PERM = Dict{Int,Vector{Int}}(2 => Int[3, 9, 7, 15, 4, 13, 8, 17, 11, 2, 10, 12, 16, 6, 14, 1, 5], 6 => Int[6, 3, 14, 7, 4, 2, 17, 9, 13, 11, 12, 5, 1, 15, 10, 16, 8], 8 => Int[14, 16, 5, 17, 13, 8, 6, 10, 12, 7, 9, 1, 15, 11, 2, 4, 3])

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
        cardinality = Float64[Inf, 17.0, 17.0, Inf, 17.0, 17.0, Inf, 17.0, 17.0, Inf],
        f_star = -45.778452053828865,
        f_star_justification = "A",
        n_minimizers_documented = 1,
        description = "DFO Benchmark Suite v2.1.0 core (Julia): CDK adaptation of paviani (F = f o T).",
        x0_idx = Float64[4.4604112528574476, 7.0, 6.0, 332.0625000158723, 5.0, 12.0, 0.09824289615362723, 10.0, 2.0, 8.22461451589162],
        x0_alt_idx = Vector{Float64}[Float64[8.776713728080988, 8.0, 17.0, 361.9147509730172, 16.0, 10.0, 0.08015446009005348, 14.0, 14.0, 9.378926330137618], Float64[5.339658053365946, 7.0, 9.0, 391.7670019301621, 12.0, 8.0, 0.05828282171252964, 4.0, 9.0, 3.937170634406482], Float64[7.098151654382945, 7.0, 2.0, 326.0920498244433, 7.0, 5.0, 0.04019438564895588, 9.0, 5.0, 5.09148244865248], Float64[3.6211302160084253, 8.0, 12.0, 355.9443007815882, 4.0, 3.0, 0.08642038892253327, 14.0, 1.0, 6.245794262898479], Float64[9.216337128335237, 7.0, 5.0, 385.7965517387331, 14.0, 17.0, 0.06833195285895952, 3.0, 13.0, 7.4001060771444775], Float64[5.779281453620196, 7.0, 15.0, 402.0020594011832, 11.0, 14.0, 0.050243516795385745, 8.0, 9.0, 8.554417891390475], Float64[7.537775054637194, 8.0, 8.0, 336.3271072954644, 6.0, 12.0, 0.09268631775501307, 13.0, 5.0, 9.708729705636475], Float64[4.040770734432936, 6.0, 3.0, 366.1793582526093, 2.0, 9.0, 0.07459788169143933, 2.0, 17.0, 4.266974009905338], Float64[8.357073209656477, 7.0, 11.0, 396.0316092097542, 13.0, 7.0, 0.05650944562786554, 7.0, 12.0, 5.421285824151337]],
        suite_release = "2.1.0",
        base_function = "paviani",
        m = 17,
        pattern = "CDK",
        strategy = "moderate",
        scale_factors_C = Float64[1.0, 100.0, 0.01, 1.0],
        lb_orig = Float64[3.1781763354492187, 3.3917369212010366, 3.20312627328, 3.172785281132861, 3.4030945768731646, 3.197009175449039, 3.545224663737264, 3.2314999192762484, 3.5588647308292924, 3.395080028701464],
        ub_orig = Float64[9.999, 9.999, 9.999, 9.999, 9.999, 9.999, 9.999, 9.999, 9.999, 9.999],
        lb_C_work = Float64[3.178176335449219, 3.3917369212010366, 3.20312627328, 317.2785281132862, 3.4030945768731646, 3.197009175449039, 0.035452246637372634, 3.2314999192762484, 3.5588647308292924, 3.3950800287014644],
        ub_C_work = Float64[9.999, 9.999, 9.999, 999.9000000000001, 9.999, 9.999, 0.09999000000000001, 9.999, 9.999, 9.999],
        ProblemData = _problem_data(),
        func_idx = _objective,
        func_F = xr -> _objective(_cell2idx(xr)),
        idx2cell = _idx2cell,
        cell2idx = _cell2idx,
        decode = xr -> _decode(_cell2idx(xr)))
end

end # module
