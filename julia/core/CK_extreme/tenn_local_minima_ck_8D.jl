# tenn_local_minima_ck_8D  DFO Benchmark Suite v2.1 core, Julia version (self-contained; Base only).
#TENN_LOCAL_MINIMA_CK_8D  DFO Benchmark Suite v2.1.0 core: tenn_local_minima, dimension 8, pattern CK, strategy extreme.
#
# Copyright (C) 2026 J. F. A. Madeira
#
# Objective F(z) = f(T(z)), where f is the original function (embedded below, translated from the MATLAB original) and T acts on the finite
# coordinates only: uniform grid of 25 levels on the interval of the domain, piecewise-affine map with |T - Id| <=
# half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
# Global optimum of F on this domain: f* = 5.889542894199896e-32 (E: exact by structure); the optimal solution is listed in
# reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
# Variables: x1:C x2:K x3:C x4:K x5:C x6:K x7:C x8:K
#
# Julia version of the MATLAB wrapper tenn_local_minima_ck_8D.m: same instance, same values (level indices 1..m as in MATLAB).
# Usage:  include("tenn_local_minima_ck_8D.jl");  PB = tenn_local_minima_ck_8D.problem();  f = PB.func_idx(PB.x0_idx)
# Copyright (C) 2026 J. F. A. Madeira. LGPL-3.0-or-later.
module tenn_local_minima_ck_8D

const NAME = "tenn_local_minima_ck_8D"
const N = 8
const TYPES = Char['C', 'K', 'C', 'K', 'C', 'K', 'C', 'K']
const CARD = Float64[Inf, 25.0, Inf, 25.0, Inf, 25.0, Inf, 25.0]
const LB = Float64[-10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0]
const UB = Float64[7.1711669921875, 6.722013920642178, 7.2627328, 6.580008329862557, 6.216201079161259, 7.20804593676692, 6.952313789346391, 6.55068177807107]
const LBW = Float64[-10.0, -10.0, -10.0, -10.0, -1000000000.0, -10.0, -1000000000.0, -10.0]
const UBW = Float64[7.1711669921875005, 6.722013920642178, 7.262732800000001, 6.580008329862557, 621620107.9161258, 7.20804593676692, 695231378.9346392, 6.55068177807107]
const TB = Float64[NaN, 1.1480092804281181, NaN, 1.0533388865750393, NaN, 0.7550287104793281, NaN, 1.0337878520473804]
const TI = Float64[NaN, 1.0, NaN, 1.0, NaN, 1.0, NaN, 1.0]
const CW = Float64[1.0000000000000018, NaN, 1.0, NaN, 100000000.0, NaN, 100000000.0, NaN]
const CX = Float64[1.0, NaN, 1.0, NaN, 1.0, NaN, 1.0, NaN]
const M = 25
const PERM = Dict{Int,Vector{Int}}(2 => Int[4, 11, 19, 14, 3, 5, 10, 25, 7, 23, 17, 2, 6, 22, 18, 21, 24, 12, 16, 1, 15, 20, 13, 8, 9], 4 => Int[18, 6, 4, 23, 5, 8, 12, 1, 16, 14, 19, 9, 13, 24, 17, 11, 7, 20, 25, 21, 15, 3, 22, 10, 2], 6 => Int[17, 12, 24, 13, 1, 2, 11, 14, 25, 5, 19, 15, 21, 22, 3, 23, 6, 18, 8, 9, 16, 4, 7, 20, 10], 8 => Int[10, 17, 11, 13, 16, 4, 25, 7, 8, 23, 2, 22, 24, 12, 19, 18, 5, 20, 21, 3, 14, 9, 15, 1, 6])

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
    d = length(x)
    s = 0.0
    for i in 1:d-1; s += ((x[i] - 1)^2) * (1 + 10 * sin(pi * x[i+1])^2); end
    f = pi / d * s
    return f + pi / d * (10 * sin(pi * x[1])^2 + (x[d] - 1)^2)
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
        cardinality = Float64[Inf, 25.0, Inf, 25.0, Inf, 25.0, Inf, 25.0],
        f_star = 5.889542894199896e-32,
        f_star_justification = "E",
        n_minimizers_documented = 1,
        description = "DFO Benchmark Suite v2.1.0 core (Julia): CK adaptation of tenn_local_minima (F = f o T).",
        x0_idx = Float64[-7.576920282840728, 14.0, -8.254392459264, 13.0, -293358630.6601405, 14.0, -336123514.5389283, 21.0],
        x0_alt_idx = Vector{Float64}[Float64[3.2892088294029254, 16.0, 3.7604695695360046, 14.0, -801409408.2717291, 11.0, -811257568.8424084, 4.0], Float64[-5.363449537754058, 17.0, -4.774225526783999, 12.0, 444959765.581161, 7.0, 303599460.53634596, 11.0], Float64[-0.9365080475807179, 15.0, 7.240636502016002, 13.0, -63091012.030427694, 4.0, -171534593.76713443, 19.0], Float64[-9.68977872133255, 15.0, -1.4266363822079988, 14.0, -424940127.0199765, 25.0, -646668648.0706146, 3.0], Float64[6.810639560222626, 17.0, -9.961331478528, 11.0, -946281964.8699323, 20.0, 567562824.0382793, 11.0], Float64[-1.8420188069343553, 15.0, 2.0535305502720025, 12.0, 300087208.9829576, 17.0, 92428769.73479915, 19.0], Float64[2.584922683238984, 16.0, -6.481164546047999, 13.0, -207963568.62863088, 14.0, -482079727.2988208, 1.0], Float64[-6.218654143810272, 16.0, 5.533697482752002, 14.0, -569812683.6181796, 9.0, -957213781.602301, 9.0], Float64[4.6474749684333805, 15.0, -3.000997613567999, 12.0, 530354827.61267066, 6.0, 257017690.50659275, 17.0]],
        suite_release = "2.1.0",
        base_function = "tenn_local_minima",
        m = 25,
        pattern = "CK",
        strategy = "extreme",
        scale_factors_C = Float64[1.0, 1.0, 100000000.0, 100000000.0],
        lb_orig = Float64[-10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0],
        ub_orig = Float64[7.1711669921875, 6.722013920642178, 7.2627328, 6.580008329862557, 6.216201079161259, 7.20804593676692, 6.952313789346391, 6.55068177807107],
        lb_C_work = Float64[-10.0, -10.0, -10.0, -10.0, -1000000000.0, -10.0, -1000000000.0, -10.0],
        ub_C_work = Float64[7.1711669921875005, 6.722013920642178, 7.262732800000001, 6.580008329862557, 621620107.9161258, 7.20804593676692, 695231378.9346392, 6.55068177807107],
        ProblemData = _problem_data(),
        func_idx = _objective,
        func_F = xr -> _objective(_cell2idx(xr)),
        idx2cell = _idx2cell,
        cell2idx = _cell2idx,
        decode = xr -> _decode(_cell2idx(xr)))
end

end # module
