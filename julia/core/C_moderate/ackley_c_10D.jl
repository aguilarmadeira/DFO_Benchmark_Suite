# ackley_c_10D  DFO Benchmark Suite v2.1 core, Julia version (self-contained; Base only).
#ACKLEY_C_10D  DFO Benchmark Suite v2.1.0 core: ackley, dimension 10, pattern C, strategy moderate.
#
# Copyright (C) 2026 J. F. A. Madeira
#
# Objective F(z) = f(T(z)), where f is the original function (embedded below, translated from the MATLAB original) and T acts on the finite
# coordinates only: uniform grid of 0 levels on the interval of the domain, piecewise-affine map with |T - Id| <=
# half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
# Global optimum of F on this domain: f* = 4.440892098500626e-16 (E: exact by structure); the optimal solution is listed in
# reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
# Variables: x1:C x2:C x3:C x4:C x5:C x6:C x7:C x8:C x9:C x10:C
#
# Julia version of the MATLAB wrapper ackley_c_10D.m: same instance, same values (level indices 1..m as in MATLAB).
# Usage:  include("ackley_c_10D.jl");  PB = ackley_c_10D.problem();  f = PB.func_idx(PB.x0_idx)
# Copyright (C) 2026 J. F. A. Madeira. LGPL-3.0-or-later.
module ackley_c_10D

const NAME = "ackley_c_10D"
const N = 10
const TYPES = Char['C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C']
const CARD = Float64[Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf]
const LB = Float64[-19.659594726562503, -30.0, -21.3513984, -21.050229071220325, -19.970767024110373, -30.0, -18.69154344416374, -20.12905364446252, -20.512336648311006, -20.415744802984953]
const UB = Float64[30.0, 18.374272214601433, 30.0, 30.0, 30.0, 19.809936626868804, 30.0, 30.0, 30.0, 30.0]
const LBW = Float64[-19.659594726562503, -3000.0, -0.21351398400000002, -21.05022907122033, -1997.0767024110373, -0.3, -18.69154344416374, -2012.905364446252, -0.20512336648311008, -20.415744802984953]
const UBW = Float64[30.0, 1837.4272214601433, 0.30000000000000004, 30.0, 3000.0, 0.19809936626868804, 30.0, 3000.0, 0.3, 30.0]
const TB = Float64[NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN]
const TI = Float64[NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN]
const CW = Float64[0.0, 4.547473508864641e-13, 2.7755575615628914e-17, -3.552713678800501e-15, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
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
    n = length(x)
    s1 = 0.0
    for v in x; s1 += v^2; end
    s2 = 0.0
    for v in x; s2 += cos(2 * pi * v); end
    return -20 * exp(-0.02 * sqrt(1 / n * s1)) - exp(1 / n * s2) + 20 + exp(1)
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
        f_star = 4.440892098500626e-16,
        f_star_justification = "E",
        n_minimizers_documented = 1,
        description = "DFO Benchmark Suite v2.1.0 core (Julia): C adaptation of ackley (F = f o T).",
        x0_idx = Float64[9.777293944358824, 293.5204082094888, -0.10818199560192002, -16.98917378067344, -1843.147217212801, -0.2634984078428499, 25.48069533675175, 1612.8452570755226, 0.07755806709817509, 5.963168618266064],
        x0_alt_idx = Vector{Float64}[Float64[-15.755691039562228, 432.87016356019285, 0.24922373726208005, -15.777235814384577, 1997.5811573675833, 0.13167100290194905, 11.833584544442882, -1761.4927029256935, -0.0732690217843053, 14.7754331676582], Float64[24.98069525957108, 731.4767821688438, -0.004657576427519988, -20.688773888292072, 432.00566157088724, 0.06977699880939014, -1.8135262478659904, -399.91773660941317, -0.20217574101912547, 23.58769771705034], Float64[-0.04308489561081075, 393.0559477457059, -0.14022526820352002, -18.456256581970482, -683.0445116943858, -0.030205623186281783, -18.314934722226386, 1224.763503391173, 0.17308604164223979, -17.955835158655457], Float64[12.759779369831087, 479.32008201042845, 0.09492305534976006, -16.223739275648896, 2707.158540354576, -0.09209962727884072, 16.561015080340724, 2849.4447433917585, 0.0441793224074197, -9.14357060926332], Float64[-12.554974973201755, 618.6698373611312, -0.04064463642624, -13.991221969327306, 1592.1083670893022, -0.15399363137139957, 2.9139042880318513, -524.8932166094582, -0.08472739682740044, -0.33130605987118145], Float64[18.87023731470108, 280.24900293799374, 0.19844747452416003, -19.92333938326753, -14.423976210250885, -0.2539762533670716, -10.733206504277016, 1099.7880233911283, 0.2905343858339649, 8.480958489520958], Float64[-6.153542840480807, 578.8556215466447, 0.06287978274816, -17.69082207694594, -1129.4741494755235, 0.17928177528084044, 24.142743298290092, 2724.469263391713, 0.16162766659914468, 17.293223038913098], Float64[6.649321424961087, 647.4245487827056, -0.19100153094144, -15.45830477062435, 2260.728902573437, 0.11738777118828159, 7.6413348239297, -912.9749702938086, 0.0327209473643246, 26.105487588305238], Float64[-18.810920011997226, 309.0037143595673, 0.16640420192256003, -20.369842844531846, 1145.6787293081638, 0.017405149192609604, -6.005775968379174, 711.7062697067768, -0.09618577187049555, -15.438045287400561]],
        suite_release = "2.1.0",
        base_function = "ackley",
        m = 0,
        pattern = "C",
        strategy = "moderate",
        scale_factors_C = Float64[1.0, 100.0, 0.01, 1.0, 100.0, 0.01, 1.0, 100.0, 0.01, 1.0],
        lb_orig = Float64[-19.659594726562503, -30.0, -21.3513984, -21.050229071220325, -19.970767024110373, -30.0, -18.69154344416374, -20.12905364446252, -20.512336648311006, -20.415744802984953],
        ub_orig = Float64[30.0, 18.374272214601433, 30.0, 30.0, 30.0, 19.809936626868804, 30.0, 30.0, 30.0, 30.0],
        lb_C_work = Float64[-19.659594726562503, -3000.0, -0.21351398400000002, -21.05022907122033, -1997.0767024110373, -0.3, -18.69154344416374, -2012.905364446252, -0.20512336648311008, -20.415744802984953],
        ub_C_work = Float64[30.0, 1837.4272214601433, 0.30000000000000004, 30.0, 3000.0, 0.19809936626868804, 30.0, 3000.0, 0.3, 30.0],
        ProblemData = _problem_data(),
        func_idx = _objective,
        func_F = xr -> _objective(_cell2idx(xr)),
        idx2cell = _idx2cell,
        cell2idx = _cell2idx,
        decode = xr -> _decode(_cell2idx(xr)))
end

end # module
