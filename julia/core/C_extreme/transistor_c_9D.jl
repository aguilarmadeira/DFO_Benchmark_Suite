# transistor_c_9D  DFO Benchmark Suite v2.1 core, Julia version (self-contained; Base only).
#TRANSISTOR_C_9D  DFO Benchmark Suite v2.1.0 core: transistor, dimension 9, pattern C, strategy extreme.
#
# Copyright (C) 2026 J. F. A. Madeira
#
# Objective F(z) = f(T(z)), where f is the original function (embedded below, translated from the MATLAB original) and T acts on the finite
# coordinates only: uniform grid of 0 levels on the interval of the domain, piecewise-affine map with |T - Id| <=
# half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
# Global optimum of F on this domain: f* = 1.60049326654736e-24 (E: exact by structure); the optimal solution is listed in
# reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
# Variables: x1:C x2:C x3:C x4:C x5:C x6:C x7:C x8:C x9:C
#
# Julia version of the MATLAB wrapper transistor_c_9D.m: same instance, same values (level indices 1..m as in MATLAB).
# Usage:  include("transistor_c_9D.jl");  PB = transistor_c_9D.problem();  f = PB.func_idx(PB.x0_idx)
# Copyright (C) 2026 J. F. A. Madeira. LGPL-3.0-or-later.
module transistor_c_9D

const NAME = "transistor_c_9D"
const N = 9
const TYPES = Char['C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C']
const CARD = Float64[Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf]
const LB = Float64[-10.0, -10.0, -10.0, -10.0, -7.350085376681921, -6.992661321382304, -10.0, -10.0, -10.0]
const UB = Float64[6.8211669921875, 7.03312503175329, 6.7027328, 6.180008329862558, 10.0, 10.0, 6.7876079069934505, 6.403313357018439, 7.316462562669516]
const LBW = Float64[-10.0, -10.0, -10.0, -10.0, -735008537.6681921, -699266132.1382303, -1000000000.0, -1000000000.0, -999999999.9999999]
const UBW = Float64[6.8211669921875, 7.03312503175329, 6.7027328, 6.1800083298625585, 1000000000.0, 1000000000.0, 678760790.6993451, 640331335.701844, 731646256.2669514]
const TB = Float64[NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN]
const TI = Float64[NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN]
const CW = Float64[0.8999999526168558, 0.4499874719815846, 1.0000064824653876, 2.0000685416236834, 799997144.0507853, 799969268.4218543, 500003127.5929673, 99998772.34567904, 200005248.3486191]
const CX = Float64[0.8999999526168565, 0.4499874719815851, 1.000006482465387, 2.000068541623683, 7.999971440507853, 7.999692684218543, 5.000031275929672, 0.9999877234567898, 2.00005248348619]
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

const _G = [0.485 0.369 5.2095 23.3037 28.5132; 0.752 1.254 10.0677 101.779 111.8467;
            0.869 0.703 22.9274 111.461 134.3884; 0.982 1.455 20.2153 191.267 211.4823]
function _orig(x::Vector{Float64})
    s = 0.0
    for k in 1:4
        g = _G[k, :]
        alfa = (1 - x[1] * x[2]) * x[3] * (exp(x[5] * (g[1] - g[3] * x[7] * 0.001 - g[5] * x[8] * 0.001)) - 1) - g[5] + g[4] * x[2]
        beta = (1 - x[1] * x[2]) * x[4] * (exp(x[6] * (g[1] - g[2] - g[3] * x[7] * 0.001 + g[4] * x[9] * 0.001)) - 1) - g[5] * x[1] + g[4]
        s += alfa^2 + beta^2
    end
    return (x[1] * x[3] - x[2] * x[4])^2 + s
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
        cardinality = Float64[Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf],
        f_star = 1.6004932665473596e-24,
        f_star_justification = "E",
        n_minimizers_documented = 1,
        description = "DFO Benchmark Suite v2.1.0 core (Julia): C adaptation of transistor (F = f o T).",
        x0_idx = Float64[-3.4210181832313538, -9.267894946051756, -1.6299265392639999, 2.8981824003985572, 336499364.63327587, 511954060.36448646, -145072359.38738823, 545388699.709187, -698559482.9642966],
        x0_alt_idx = Vector{Float64}[Float64[0.9156889319419861, -8.216467474955873, -9.887757635584, 3.605762939605375, -207075812.0817024, 300802647.495239, -615590089.6526027, -644861928.3398107, 587899834.0733439], Float64[-7.560602247714996, -9.408085275531207, 1.7373443932160004, 2.0490857533503757, -594226477.5117949, -40288096.370468736, 488435172.8651228, -113230996.82425463, 145986328.22071946], Float64[3.084042489528656, -9.104339561659064, -6.648763691008, 2.756666292557192, 582867969.9069712, -251439509.2397163, 17917442.599908352, 418399934.69130135, -295927177.6319053], Float64[-5.490810215473175, -8.613673408480985, 4.976338337791999, 3.4642468317640116, 195717304.47687876, -462590922.10896385, -452600287.6653062, -685756615.3794689, -737840683.4845297], Float64[-1.1541031002998352, -9.805291209056318, -3.2814927585280005, 1.5841042561573229, -362078282.35490745, 885529636.9793091, -923118017.9305205, -154125683.86391282, 548618633.5531107], Float64[-9.630394279956818, -8.753863737960437, -7.691014217728, 2.291684795364141, 971440676.3487444, 674378224.1100614, 279316181.6361382, 377505247.6516433, 106705127.70048606], Float64[6.533695876598358, -8.51242483719027, 0.08577817395200071, 2.999265334570957, 427865499.63376594, 463226811.240814, -289610485.67800975, -812745380.3973547, -335208378.1521385], Float64[-2.090437591075897, -9.704042637765603, -4.323743285248001, 3.706845873777775, 40714834.20367348, 122136067.37510633, -760128215.9432241, -281114448.88179874, -777121884.0047631], Float64[2.2462695240974426, -9.213376484587526, 3.453049106432001, 2.150168687522779, -502860342.51130474, -89015345.49414122, 442305983.6234348, 250516482.63375735, 434190788.5593878]],
        suite_release = "2.1.0",
        base_function = "transistor",
        m = 0,
        pattern = "C",
        strategy = "extreme",
        scale_factors_C = Float64[1.0, 1.0, 1.0, 1.0, 100000000.0, 100000000.0, 100000000.0, 100000000.0, 100000000.0],
        lb_orig = Float64[-10.0, -10.0, -10.0, -10.0, -7.350085376681921, -6.992661321382304, -10.0, -10.0, -10.0],
        ub_orig = Float64[6.8211669921875, 7.03312503175329, 6.7027328, 6.180008329862558, 10.0, 10.0, 6.7876079069934505, 6.403313357018439, 7.316462562669516],
        lb_C_work = Float64[-10.0, -10.0, -10.0, -10.0, -735008537.6681921, -699266132.1382303, -1000000000.0, -1000000000.0, -999999999.9999999],
        ub_C_work = Float64[6.8211669921875, 7.03312503175329, 6.7027328, 6.1800083298625585, 1000000000.0, 1000000000.0, 678760790.6993451, 640331335.701844, 731646256.2669514],
        ProblemData = _problem_data(),
        func_idx = _objective,
        func_F = xr -> _objective(_cell2idx(xr)),
        idx2cell = _idx2cell,
        cell2idx = _cell2idx,
        decode = xr -> _decode(_cell2idx(xr)))
end

end # module
