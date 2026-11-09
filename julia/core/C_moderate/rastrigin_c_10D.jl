# rastrigin_c_10D  DFO Benchmark Suite v2.1 core, Julia version (self-contained; Base only).
#RASTRIGIN_C_10D  DFO Benchmark Suite v2.1.0 core: rastrigin, dimension 10, pattern C, strategy moderate.
#
# Copyright (C) 2026 J. F. A. Madeira
#
# Objective F(z) = f(T(z)), where f is the original function (embedded below, translated from the MATLAB original) and T acts on the finite
# coordinates only: uniform grid of 0 levels on the interval of the domain, piecewise-affine map with |T - Id| <=
# half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
# Global optimum of F on this domain: f* = 0 (E: exact by structure); the optimal solution is listed in
# reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
# Variables: x1:C x2:C x3:C x4:C x5:C x6:C x7:C x8:C x9:C x10:C
#
# Julia version of the MATLAB wrapper rastrigin_c_10D.m: same instance, same values (level indices 1..m as in MATLAB).
# Usage:  include("rastrigin_c_10D.jl");  PB = rastrigin_c_10D.problem();  f = PB.func_idx(PB.x0_idx)
# Copyright (C) 2026 J. F. A. Madeira. LGPL-3.0-or-later.
module rastrigin_c_10D

const NAME = "rastrigin_c_10D"
const N = 10
const TYPES = Char['C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C']
const CARD = Float64[Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf]
const LB = Float64[-5.12, -5.12, -3.6955815936, -3.210139775093711, -3.775091646745441, -3.3681715906305802, -5.12, -5.12, -5.12, -5.12]
const UB = Float64[3.7836375, 3.2332190011685213, 5.12, 5.12, 5.12, 5.12, 3.6885590546090206, 3.393660981729729, 3.7162186241472837, 3.23626238058141]
const LBW = Float64[-5.12, -512.0, -0.036955815936, -3.210139775093711, -377.50916467454414, -0.0336817159063058, -5.12, -512.0, -0.0512, -5.119999999999999]
const UBW = Float64[3.7836375, 323.3219001168522, 0.05120000000000001, 5.12, 512.0, 0.0512, 3.68855905460902, 339.3660981729729, 0.037162186241472835, 3.2362623805814095]
const TB = Float64[NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN]
const TI = Float64[NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN]
const CW = Float64[0.0, 0.0, 0.0, 0.0, -1.1368683772161603e-13, 0.0, 0.0, 0.0, 0.0, 8.881784197001252e-16]
const CX = Float64[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
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
    for v in x; s += v^2 - 10 * cos(2 * pi * v); end
    return 10 * length(x) + s
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
        f_star = 0.0,
        f_star_justification = "E",
        n_minimizers_documented = 1,
        description = "DFO Benchmark Suite v2.1.0 core (Julia): C adaptation of rastrigin (F = f o T).",
        x0_idx = Float64[-4.837413848876953, 190.78568642661537, -0.02522051371859968, 1.7303512868527466, 157.13195609152302, -0.03244538686260904, -3.255375245920338, 153.30431348697107, -0.028308768715942935, 2.1436337065203945],
        x0_alt_idx = Vector{Float64}[Float64[3.301067303466797, 214.84845721187446, 0.019210017513144335, 2.0946431137518884, -121.54972810778958, 0.04138663821565696, 3.0538694748549817, -464.4605313310616, 0.03350191321889518, -4.742085472722075], Float64[-1.1855312805175782, 156.41029959053094, -0.007448301225902076, 2.45893494065103, -320.03524419219195, 0.030839206061618947, 0.068676959910138, -188.53301474868533, 0.010952016919086413, -3.2814783740948847], Float64[1.109937762451172, 171.30630055283427, 0.03698223000584193, 1.657492921472918, 283.4409208725063, 0.013801046428172967, -2.400157930828031, 87.39450183369092, -0.011597879380722334, -1.8208712754676948], Float64[-3.4288305725097654, 222.86938080696098, 0.009646874600407043, 2.0217847483720606, 84.95540478810398, 0.0032536142741349597, -4.8689928215661995, -485.68572491432127, -0.034147775680531095, -0.36026417684050394], Float64[2.2055025329589846, 164.43122318561734, -0.03055217746640896, 2.3860765752712023, -201.01684014892768, -0.013784545359311038, 1.4402518992091196, -209.75820833194507, 0.03149747799224551, 1.1003429217866874], Float64[-2.281096051025391, 188.49399397087643, 0.027419087093104658, 2.5838349955878783, 482.6554930306807, -0.02433197751334904, -1.0285829915290492, 66.16930825043119, 0.008947581692436757, 2.5609500204138778], Float64[0.014372991943359281, 200.33440499219444, -0.01277996497371136, 1.7823929764097661, 203.97380883136822, 0.04950004756491696, -4.013775506473893, 297.4122067627869, -0.013602314607372004, -4.324769158828592], Float64[-4.550480218505859, 141.89624737085092, 0.045191299585802244, 2.146684803308907, 5.488292746965726, 0.03246188793147097, 2.295469214301427, -275.6680199852252, -0.03615221090718076, -2.864162060201402], Float64[2.7532849182128905, 193.4593276249775, 0.004992247518986234, 2.5109766302080496, -273.19339145234676, 0.021914455777432958, -0.17336567643674172, 0.2594965971510419, 0.025658471027657344, -1.4035549615742111]],
        suite_release = "2.1.0",
        base_function = "rastrigin",
        m = 0,
        pattern = "C",
        strategy = "moderate",
        scale_factors_C = Float64[1.0, 100.0, 0.01, 1.0, 100.0, 0.01, 1.0, 100.0, 0.01, 1.0],
        lb_orig = Float64[-5.12, -5.12, -3.6955815936, -3.210139775093711, -3.775091646745441, -3.3681715906305802, -5.12, -5.12, -5.12, -5.12],
        ub_orig = Float64[3.7836375, 3.2332190011685213, 5.12, 5.12, 5.12, 5.12, 3.6885590546090206, 3.393660981729729, 3.7162186241472837, 3.23626238058141],
        lb_C_work = Float64[-5.12, -512.0, -0.036955815936, -3.210139775093711, -377.50916467454414, -0.0336817159063058, -5.12, -512.0, -0.0512, -5.119999999999999],
        ub_C_work = Float64[3.7836375, 323.3219001168522, 0.05120000000000001, 5.12, 512.0, 0.0512, 3.68855905460902, 339.3660981729729, 0.037162186241472835, 3.2362623805814095],
        ProblemData = _problem_data(),
        func_idx = _objective,
        func_F = xr -> _objective(_cell2idx(xr)),
        idx2cell = _idx2cell,
        cell2idx = _cell2idx,
        decode = xr -> _decode(_cell2idx(xr)))
end

end # module
