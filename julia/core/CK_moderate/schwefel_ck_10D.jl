# schwefel_ck_10D  DFO Benchmark Suite v2.1 core, Julia version (self-contained; Base only).
#SCHWEFEL_CK_10D  DFO Benchmark Suite v2.1.0 core: schwefel, dimension 10, pattern CK, strategy moderate.
#
# Copyright (C) 2026 J. F. A. Madeira
#
# Objective F(z) = f(T(z)), where f is the original function (embedded below, translated from the MATLAB original) and T acts on the finite
# coordinates only: uniform grid of 17 levels on the interval of the domain, piecewise-affine map with |T - Id| <=
# half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
# Global optimum of F on this domain: f* = -4189.828872721624 (A: analytic reduction, 1-D root to machine precision); the optimal solution is listed in
# reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
# Variables: x1:C x2:K x3:C x4:K x5:C x6:K x7:C x8:K x9:C x10:K
#
# Julia version of the MATLAB wrapper schwefel_ck_10D.m: same instance, same values (level indices 1..m as in MATLAB).
# Usage:  include("schwefel_ck_10D.jl");  PB = schwefel_ck_10D.problem();  f = PB.func_idx(PB.x0_idx)
# Copyright (C) 2026 J. F. A. Madeira. LGPL-3.0-or-later.
module schwefel_ck_10D

const NAME = "schwefel_ck_10D"
const N = 10
const TYPES = Char['C', 'K', 'C', 'K', 'C', 'K', 'C', 'K', 'C', 'K']
const CARD = Float64[Inf, 17.0, Inf, 17.0, Inf, 17.0, Inf, 17.0, Inf, 17.0]
const LB = Float64[-343.245849609375, -300.1884875273078, -344.09664, -322.0616409829238, -330.47947544566625, -366.20111340639335, -335.504962823721, -309.3069420891491, -344.6511054491658, -301.55848948296364]
const UB = Float64[500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0]
const LBW = Float64[-343.245849609375, -300.1884875273078, -34409.664000000004, -322.0616409829238, -3.3047947544566623, -366.20111340639335, -335.5049628237209, -309.3069420891491, -34465.11054491658, -301.55848948296364]
const UBW = Float64[500.0, 500.0, 50000.0, 500.0, 5.0, 500.0, 500.0, 500.0, 50000.0, 500.0]
const TB = Float64[NaN, 399.97643905908654, NaN, 397.24229487713455, NaN, 445.86243041210037, NaN, 398.83663223885634, NaN, 399.80518881462956]
const TI = Float64[NaN, 420.9687, NaN, 420.9687, NaN, 420.9687, NaN, 420.9687, NaN, 420.9687]
const CW = Float64[420.9687, NaN, 42096.869999999995, NaN, 4.209687000000001, NaN, 420.9687000000001, NaN, 42096.869999999995, NaN]
const CX = Float64[420.9687, NaN, 420.9687, NaN, 420.9687, NaN, 420.9687, NaN, 420.9687, NaN]
const M = 17
const PERM = Dict{Int,Vector{Int}}(2 => Int[8, 5, 16, 6, 10, 13, 12, 1, 17, 2, 11, 7, 3, 15, 14, 9, 4], 4 => Int[2, 16, 9, 14, 8, 6, 17, 7, 10, 12, 3, 13, 15, 1, 4, 11, 5], 6 => Int[1, 17, 11, 8, 7, 13, 9, 15, 3, 16, 6, 12, 14, 4, 10, 5, 2], 8 => Int[11, 12, 13, 14, 5, 10, 15, 3, 2, 16, 9, 17, 4, 6, 1, 7, 8], 10 => Int[8, 4, 10, 16, 13, 12, 6, 11, 2, 15, 1, 5, 17, 7, 9, 14, 3])

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
    for v in x; s += v * sin(sqrt(abs(v))); end
    return -s
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
        f_star = -4189.828872721624,
        f_star_justification = "A",
        n_minimizers_documented = 1,
        description = "DFO Benchmark Suite v2.1.0 core (Julia): CK adaptation of schwefel (F = f o T).",
        x0_idx = Float64[-0.2654820680618286, 13.0, -2914.730168320002, 9.0, -2.0194546196572087, 9.0, 136.24158040302484, 3.0, 9450.915534408014, 2.0],
        x0_alt_idx = Vector{Float64}[Float64[217.13383853435516, 12.0, 39627.740487680014, 9.0, 3.6148276217209774, 7.0, -97.93109083822569, 7.0, -15769.906904129173, 5.0], Float64[-207.78301537036896, 13.0, 14102.258094079996, 10.0, 1.761691602131474, 5.0, 451.5328893945125, 12.0, 46980.165769126426, 8.0], Float64[325.83349883556366, 13.0, -24388.548689920004, 8.0, -0.8401862435547982, 1.0, 217.36021815326194, 1.0, 25424.79729358061, 11.0], Float64[-104.0242487192154, 12.0, 30470.98013696001, 9.0, -2.693322263144301, 15.0, -16.812453087988615, 7.0, 3869.428818034794, 13.0], Float64[113.3750718832016, 13.0, -8019.826647040001, 10.0, 2.872892539497815, 12.0, -250.98512432923908, 12.0, -17685.939657511022, 16.0], Float64[-311.5417820215225, 13.0, 47487.968399360005, 8.0, 1.0197565199083125, 10.0, 347.4561466206234, 17.0, 45064.13301574458, 3.0], Float64[459.2376273870468, 12.0, 8997.161615359997, 9.0, -1.5821213257779598, 8.0, 64.30618466224848, 5.0, 23508.76454019876, 6.0], Float64[26.909433007240295, 12.0, -32734.976266240003, 9.0, 4.800902741696996, 5.0, -169.86648657900201, 10.0, 1953.3960646529376, 9.0], Float64[244.3087536096573, 13.0, 26014.14987776, 10.0, 2.1990248960107226, 3.0, 428.57478437086047, 16.0, -23267.426373884242, 11.0]],
        suite_release = "2.1.0",
        base_function = "schwefel",
        m = 17,
        pattern = "CK",
        strategy = "moderate",
        scale_factors_C = Float64[1.0, 100.0, 0.01, 1.0, 100.0],
        lb_orig = Float64[-343.245849609375, -300.1884875273078, -344.09664, -322.0616409829238, -330.47947544566625, -366.20111340639335, -335.504962823721, -309.3069420891491, -344.6511054491658, -301.55848948296364],
        ub_orig = Float64[500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0],
        lb_C_work = Float64[-343.245849609375, -300.1884875273078, -34409.664000000004, -322.0616409829238, -3.3047947544566623, -366.20111340639335, -335.5049628237209, -309.3069420891491, -34465.11054491658, -301.55848948296364],
        ub_C_work = Float64[500.0, 500.0, 50000.0, 500.0, 5.0, 500.0, 500.0, 500.0, 50000.0, 500.0],
        ProblemData = _problem_data(),
        func_idx = _objective,
        func_F = xr -> _objective(_cell2idx(xr)),
        idx2cell = _idx2cell,
        cell2idx = _cell2idx,
        decode = xr -> _decode(_cell2idx(xr)))
end

end # module
