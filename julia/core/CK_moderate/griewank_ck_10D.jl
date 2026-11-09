# griewank_ck_10D  DFO Benchmark Suite v2.1 core, Julia version (self-contained; Base only).
#GRIEWANK_CK_10D  DFO Benchmark Suite v2.1.0 core: griewank, dimension 10, pattern CK, strategy moderate.
#
# Copyright (C) 2026 J. F. A. Madeira
#
# Objective F(z) = f(T(z)), where f is the original function (embedded below, translated from the MATLAB original) and T acts on the finite
# coordinates only: uniform grid of 17 levels on the interval of the domain, piecewise-affine map with |T - Id| <=
# half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
# Global optimum of F on this domain: f* = 0 (E: exact by structure); the optimal solution is listed in
# reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
# Variables: x1:C x2:K x3:C x4:K x5:C x6:K x7:C x8:K x9:C x10:K
#
# Julia version of the MATLAB wrapper griewank_ck_10D.m: same instance, same values (level indices 1..m as in MATLAB).
# Usage:  include("griewank_ck_10D.jl");  PB = griewank_ck_10D.problem();  f = PB.func_idx(PB.x0_idx)
# Copyright (C) 2026 J. F. A. Madeira. LGPL-3.0-or-later.
module griewank_ck_10D

const NAME = "griewank_ck_10D"
const N = 10
const TYPES = Char['C', 'K', 'C', 'K', 'C', 'K', 'C', 'K', 'C', 'K']
const CARD = Float64[Inf, 17.0, Inf, 17.0, Inf, 17.0, Inf, 17.0, Inf, 17.0]
const LB = Float64[-270.0029296875, -256.7433826144388, -400.0, -400.0, -400.0, -400.0, -400.0, -400.0, -278.2615270814498, -287.5916191725778]
const UB = Float64[400.0, 400.0, 264.525312, 285.07788421491045, 275.49101837306193, 285.0082280032212, 291.6565654146861, 262.337520430322, 400.0, 400.0]
const LBW = Float64[-270.0029296875, -256.7433826144388, -40000.0, -400.0, -4.0, -400.0, -400.0, -400.0, -27826.15270814498, -287.5916191725778]
const UBW = Float64[400.0, 400.0, 26452.531199999998, 285.07788421491045, 2.75491018373062, 285.0082280032212, 291.65656541468616, 262.337520430322, 40000.0, 400.0]
const TB = Float64[NaN, -10.46461413402426, NaN, -14.643690129112883, NaN, -14.68287174818812, NaN, 13.960950268951251, NaN, 13.22971421542502]
const TI = Float64[NaN, 0.0, NaN, 0.0, NaN, 0.0, NaN, 0.0, NaN, 0.0]
const CW = Float64[0.0, NaN, 0.0, NaN, 8.881784197001252e-16, NaN, 5.684341886080802e-14, NaN, 0.0, NaN]
const CX = Float64[0.0, NaN, 0.0, NaN, 0.0, NaN, 0.0, NaN, 0.0, NaN]
const M = 17
const PERM = Dict{Int,Vector{Int}}(2 => Int[15, 14, 1, 9, 13, 3, 12, 10, 2, 7, 11, 16, 8, 5, 17, 6, 4], 4 => Int[15, 13, 12, 7, 17, 1, 11, 5, 2, 9, 4, 3, 16, 8, 10, 14, 6], 6 => Int[13, 5, 9, 12, 6, 4, 16, 3, 15, 14, 7, 11, 8, 10, 17, 2, 1], 8 => Int[3, 5, 6, 4, 12, 14, 7, 15, 1, 2, 10, 13, 16, 11, 9, 17, 8], 10 => Int[10, 6, 1, 11, 12, 8, 5, 7, 13, 3, 2, 17, 9, 14, 16, 15, 4])

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
        cardinality = Float64[Inf, 17.0, Inf, 17.0, Inf, 17.0, Inf, 17.0, Inf, 17.0],
        f_star = 0.0,
        f_star_justification = "E",
        n_minimizers_documented = 1,
        description = "DFO Benchmark Suite v2.1.0 core (Julia): CK adaptation of griewank (F = f o T).",
        x0_idx = Float64[64.34423542022705, 15.0, -2446.3455682559943, 3.0, 2.404730767744586, 11.0, -296.66681884502754, 14.0, 4361.585085627441, 5.0],
        x0_alt_idx = Vector{Float64}[Float64[222.3576216697693, 14.0, -32748.699795456, 1.0, 0.8974367598046946, 9.0, 198.7411709156646, 3.0, -15890.935074534525, 8.0], Float64[-107.40944528579712, 14.0, 10950.484721664005, 2.0, -1.2188649281109094, 7.0, 4.885870574524176, 8.0, 34497.870245505124, 11.0], Float64[308.7251868247986, 15.0, -19351.869505535997, 3.0, -3.335166616026513, 4.0, -229.51432918437106, 13.0, 17188.73864701822, 13.0], Float64[-21.041880130767822, 14.0, 23836.959571968, 1.0, 1.8566238557664434, 2.0, 265.89366057632105, 2.0, -120.39295146868244, 16.0], Float64[139.91585493087769, 14.0, -6465.394655232005, 2.0, -0.3150421666924079, 16.0, 72.03836023518073, 7.0, -17429.524549955582, 3.0], Float64[-189.85121202468872, 15.0, -39319.526080512, 2.0, -1.8223361746322984, 12.0, -121.81694010595976, 13.0, 32959.280770084064, 6.0], Float64[351.90896940231323, 13.0, 6931.435634688001, 3.0, -3.9386378625479015, 10.0, -356.21713986485497, 1.0, 15650.14917159716, 9.0], Float64[22.141902446746826, 14.0, -25922.695790592, 1.0, 1.2531526092450544, 7.0, 139.19084989583712, 6.0, -1658.9824268897355, 11.0], Float64[181.1367383003235, 15.0, 20328.265924607993, 2.0, -0.8631490786705496, 5.0, -54.66445044530326, 11.0, -21911.50258705171, 14.0]],
        suite_release = "2.1.0",
        base_function = "griewank",
        m = 17,
        pattern = "CK",
        strategy = "moderate",
        scale_factors_C = Float64[1.0, 100.0, 0.01, 1.0, 100.0],
        lb_orig = Float64[-270.0029296875, -256.7433826144388, -400.0, -400.0, -400.0, -400.0, -400.0, -400.0, -278.2615270814498, -287.5916191725778],
        ub_orig = Float64[400.0, 400.0, 264.525312, 285.07788421491045, 275.49101837306193, 285.0082280032212, 291.6565654146861, 262.337520430322, 400.0, 400.0],
        lb_C_work = Float64[-270.0029296875, -256.7433826144388, -40000.0, -400.0, -4.0, -400.0, -400.0, -400.0, -27826.15270814498, -287.5916191725778],
        ub_C_work = Float64[400.0, 400.0, 26452.531199999998, 285.07788421491045, 2.75491018373062, 285.0082280032212, 291.65656541468616, 262.337520430322, 40000.0, 400.0],
        ProblemData = _problem_data(),
        func_idx = _objective,
        func_F = xr -> _objective(_cell2idx(xr)),
        idx2cell = _idx2cell,
        cell2idx = _cell2idx,
        decode = xr -> _decode(_cell2idx(xr)))
end

end # module
