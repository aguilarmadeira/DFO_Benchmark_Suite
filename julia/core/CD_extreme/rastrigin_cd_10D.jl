# rastrigin_cd_10D  DFO Benchmark Suite v2.1 core, Julia version (self-contained; Base only).
#RASTRIGIN_CD_10D  DFO Benchmark Suite v2.1.0 core: rastrigin, dimension 10, pattern CD, strategy extreme.
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
# Julia version of the MATLAB wrapper rastrigin_cd_10D.m: same instance, same values (level indices 1..m as in MATLAB).
# Usage:  include("rastrigin_cd_10D.jl");  PB = rastrigin_cd_10D.problem();  f = PB.func_idx(PB.x0_idx)
# Copyright (C) 2026 J. F. A. Madeira. LGPL-3.0-or-later.
module rastrigin_cd_10D

const NAME = "rastrigin_cd_10D"
const N = 10
const TYPES = Char['C', 'D', 'C', 'D', 'C', 'D', 'C', 'D', 'C', 'D']
const CARD = Float64[Inf, 17.0, Inf, 17.0, Inf, 17.0, Inf, 17.0, Inf, 17.0]
const LB = Float64[-5.12, -5.12, -3.6955815936, -3.210139775093711, -3.775091646745441, -3.3681715906305802, -5.12, -5.12, -5.12, -5.12]
const UB = Float64[3.7836375, 3.2332190011685213, 5.12, 5.12, 5.12, 5.12, 3.6885590546090206, 3.393660981729729, 3.7162186241472837, 3.23626238058141]
const LBW = Float64[-5.12, -5.12, -3.6955815936000005, -3.210139775093711, -377509164.67454416, -3.3681715906305802, -512000000.0, -5.12, -512000000.0, -5.12]
const UBW = Float64[3.7836375, 3.2332190011685213, 5.120000000000001, 5.12, 512000000.00000006, 5.12, 368855905.460902, 3.393660981729729, 371621862.4147284, 3.23626238058141]
const TB = Float64[NaN, 0.1007618757303268, NaN, -0.08633735943356946, NaN, -0.1851072441441124, NaN, 0.20103811358108015, NaN, 0.10266398786338105]
const TI = Float64[NaN, 0.0, NaN, 0.0, NaN, 0.0, NaN, 0.0, NaN, 0.0]
const CW = Float64[0.0, NaN, -4.440892098500626e-16, NaN, -5.960464477539063e-08, NaN, 0.0, NaN, 5.960464477539063e-08, NaN]
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
        cardinality = Float64[Inf, 17.0, Inf, 17.0, Inf, 17.0, Inf, 17.0, Inf, 17.0],
        f_star = 0.0,
        f_star_justification = "E",
        n_minimizers_documented = 1,
        description = "DFO Benchmark Suite v2.1.0 core (Julia): CD adaptation of rastrigin (F = f o T).",
        x0_idx = Float64[-4.837413848876953, 14.0, -2.5220513718599684, 10.0, 157131956.09152305, 1.0, -325537524.59203374, 14.0, -283087687.1594293, 15.0],
        x0_alt_idx = Vector{Float64}[Float64[3.301067303466797, 15.0, 1.9210017513144337, 11.0, -121549728.10778958, 15.0, 305386947.4854981, 2.0, 335019132.18895185, 2.0], Float64[-1.1855312805175782, 14.0, -0.7448301225902081, 12.0, -320035244.19219196, 13.0, 6867695.991013765, 7.0, 109520169.1908642, 5.0], Float64[1.109937762451172, 14.0, 3.6982230005841923, 10.0, 283440920.8725063, 10.0, -240015793.08280313, 12.0, -115978793.80722332, 7.0], Float64[-3.4288305725097654, 15.0, 0.9646874600407043, 11.0, 84955404.788104, 8.0, -486899282.1566199, 1.0, -341477756.80531085, 10.0], Float64[2.2055025329589846, 14.0, -3.0552177466408965, 12.0, -201016840.14892766, 5.0, 144025189.9209119, 7.0, 314974779.9224551, 13.0], Float64[-2.281096051025391, 14.0, 2.7419087093104655, 12.0, 482655493.0306807, 3.0, -102858299.15290493, 12.0, 89475816.92436767, 16.0], Float64[0.014372991943359281, 15.0, -1.2779964973711362, 11.0, 203973808.83136827, 17.0, -401377550.6473893, 16.0, -136023146.07372004, 3.0], Float64[-4.550480218505859, 14.0, 4.5191299585802245, 11.0, 5488292.746965766, 13.0, 229546921.43014264, 5.0, -361522109.07180756, 5.0], Float64[2.7532849182128905, 15.0, 0.4992247518986237, 12.0, -273193391.45234674, 11.0, -17336567.643674135, 11.0, 256584710.27657354, 8.0]],
        suite_release = "2.1.0",
        base_function = "rastrigin",
        m = 17,
        pattern = "CD",
        strategy = "extreme",
        scale_factors_C = Float64[1.0, 1.0, 100000000.0, 100000000.0, 100000000.0],
        lb_orig = Float64[-5.12, -5.12, -3.6955815936, -3.210139775093711, -3.775091646745441, -3.3681715906305802, -5.12, -5.12, -5.12, -5.12],
        ub_orig = Float64[3.7836375, 3.2332190011685213, 5.12, 5.12, 5.12, 5.12, 3.6885590546090206, 3.393660981729729, 3.7162186241472837, 3.23626238058141],
        lb_C_work = Float64[-5.12, -5.12, -3.6955815936000005, -3.210139775093711, -377509164.67454416, -3.3681715906305802, -512000000.0, -5.12, -512000000.0, -5.12],
        ub_C_work = Float64[3.7836375, 3.2332190011685213, 5.120000000000001, 5.12, 512000000.00000006, 5.12, 368855905.460902, 3.393660981729729, 371621862.4147284, 3.23626238058141],
        ProblemData = _problem_data(),
        func_idx = _objective,
        func_F = xr -> _objective(_cell2idx(xr)),
        idx2cell = _idx2cell,
        cell2idx = _cell2idx,
        decode = xr -> _decode(_cell2idx(xr)))
end

end # module
