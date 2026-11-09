# schwefel_c_10D  DFO Benchmark Suite v2.1 core, Julia version (self-contained; Base only).
#SCHWEFEL_C_10D  DFO Benchmark Suite v2.1.0 core: schwefel, dimension 10, pattern C, strategy extreme.
#
# Copyright (C) 2026 J. F. A. Madeira
#
# Objective F(z) = f(T(z)), where f is the original function (embedded below, translated from the MATLAB original) and T acts on the finite
# coordinates only: uniform grid of 0 levels on the interval of the domain, piecewise-affine map with |T - Id| <=
# half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
# Global optimum of F on this domain: f* = -4189.828872721624 (A: analytic reduction, 1-D root to machine precision); the optimal solution is listed in
# reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
# Variables: x1:C x2:C x3:C x4:C x5:C x6:C x7:C x8:C x9:C x10:C
#
# Julia version of the MATLAB wrapper schwefel_c_10D.m: same instance, same values (level indices 1..m as in MATLAB).
# Usage:  include("schwefel_c_10D.jl");  PB = schwefel_c_10D.problem();  f = PB.func_idx(PB.x0_idx)
# Copyright (C) 2026 J. F. A. Madeira. LGPL-3.0-or-later.
module schwefel_c_10D

const NAME = "schwefel_c_10D"
const N = 10
const TYPES = Char['C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C']
const CARD = Float64[Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf]
const LB = Float64[-343.245849609375, -300.1884875273078, -344.09664, -322.0616409829238, -330.47947544566625, -366.20111340639335, -335.504962823721, -309.3069420891491, -344.6511054491658, -301.55848948296364]
const UB = Float64[500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0]
const LBW = Float64[-343.245849609375, -300.18848752730787, -344.09664, -322.0616409829238, -330.47947544566625, -36620111340.63933, -33550496282.372093, -30930694208.91491, -34465110544.91658, -30155848948.29636]
const UBW = Float64[500.0, 500.0, 500.0, 500.0, 500.0, 50000000000.0, 49999999999.99999, 49999999999.99999, 50000000000.0, 50000000000.0]
const TB = Float64[NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN]
const TI = Float64[NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN]
const CW = Float64[420.9687, 420.9686999999999, 420.9687, 420.9687, 420.9687, 42096870000.0, 42096870000.00001, 42096870000.00001, 42096870000.0, 42096870000.0]
const CX = Float64[420.9687, 420.9687, 420.9687, 420.9687, 420.9687, 420.9687, 420.9687, 420.9687, 420.9687, 420.9687]
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
        cardinality = Float64[Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf],
        f_star = -4189.828872721624,
        f_star_justification = "A",
        n_minimizers_documented = 1,
        description = "DFO Benchmark Suite v2.1.0 core (Julia): C adaptation of schwefel (F = f o T).",
        x0_idx = Float64[-0.2654820680618286, 313.76500221700974, -29.1473016832, 64.83117630599907, -201.9454619657209, 4620050908.932243, 13624158040.302475, -22836444868.294464, 9450915534.408012, -23957401651.708267],
        x0_alt_idx = Vector{Float64}[Float64[217.13383853435516, 257.78473765748606, 396.2774048768001, 100.78139384169555, 361.48276217209775, -6143394878.957855, -9793109083.82257, -854540317.8923035, -15769906904.129173, -9946807364.669552], Float64[-207.78301537036896, 280.83543482905463, 141.02258094080003, 136.73161137739186, 176.1691602131474, -16906840666.847946, 45153288939.45124, 25375075256.742165, 46980165769.126434, 4063786922.369171], Float64[325.83349883556366, 295.10491403050185, -243.88548689919998, 57.641132798859815, -84.01862435547983, -34293945401.131943, 21736021815.326195, -29101818629.037113, 25424797293.580605, 18074381209.407894], Float64[-104.0242487192154, 239.12464947097828, 304.7098013696001, 93.5913503345563, -269.3322263144301, 41050175114.098724, -1681245308.7988663, -2872203054.402645, 3869428818.0347977, 32084975496.446613], Float64[113.3750718832016, 288.51900055291094, -80.19826647039997, 129.5415678702526, 287.2892539497815, 23663070379.814713, -25098512432.923912, 23357412520.231827, -17685939657.511024, 46095569783.48534], Float64[-311.5417820215225, 311.5696977244794, 474.8796839936001, 34.016704132545044, 101.97565199083124, 12899624591.92463, 34745614662.062325, 49587028094.866295, 45064133015.74458, -19954374712.55435], Float64[459.2376273870468, 243.88114253812728, 89.97161615360005, 69.96692166824141, -158.21213257779598, 2136178804.0345306, 6430618466.224846, -9137576815.14529, 23508764540.19876, -5943780425.515629], Float64[26.909433007240295, 266.93183970969596, -327.3497626624, 105.91713920393778, 480.0902741696997, -15250925930.24947, -16986648657.900204, 17092038759.489182, 1953396064.6529388, 8066813861.52309], Float64[244.3087536096573, 316.3261907916285, 260.1414987776001, 141.8673567396341, 219.90248960107226, -26014371718.139565, 42857478437.08604, 43321654334.123634, -23267426373.884247, 22077408148.561817]],
        suite_release = "2.1.0",
        base_function = "schwefel",
        m = 0,
        pattern = "C",
        strategy = "extreme",
        scale_factors_C = Float64[1.0, 1.0, 1.0, 1.0, 1.0, 100000000.0, 100000000.0, 100000000.0, 100000000.0, 100000000.0],
        lb_orig = Float64[-343.245849609375, -300.1884875273078, -344.09664, -322.0616409829238, -330.47947544566625, -366.20111340639335, -335.504962823721, -309.3069420891491, -344.6511054491658, -301.55848948296364],
        ub_orig = Float64[500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0],
        lb_C_work = Float64[-343.245849609375, -300.18848752730787, -344.09664, -322.0616409829238, -330.47947544566625, -36620111340.63933, -33550496282.372093, -30930694208.91491, -34465110544.91658, -30155848948.29636],
        ub_C_work = Float64[500.0, 500.0, 500.0, 500.0, 500.0, 50000000000.0, 49999999999.99999, 49999999999.99999, 50000000000.0, 50000000000.0],
        ProblemData = _problem_data(),
        func_idx = _objective,
        func_F = xr -> _objective(_cell2idx(xr)),
        idx2cell = _idx2cell,
        cell2idx = _cell2idx,
        decode = xr -> _decode(_cell2idx(xr)))
end

end # module
