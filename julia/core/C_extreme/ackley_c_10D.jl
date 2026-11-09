# ackley_c_10D  DFO Benchmark Suite v2.1 core, Julia version (self-contained; Base only).
#ACKLEY_C_10D  DFO Benchmark Suite v2.1.0 core: ackley, dimension 10, pattern C, strategy extreme.
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
const LBW = Float64[-19.659594726562503, -30.0, -21.3513984, -21.05022907122033, -19.970767024110373, -3000000000.0, -1869154344.4163742, -2012905364.4462519, -2051233664.8311005, -2041574480.2984953]
const UBW = Float64[30.0, 18.374272214601433, 30.0, 30.0, 30.0, 1980993662.6868804, 3000000000.0, 3000000000.0, 3000000000.0, 3000000000.0]
const TB = Float64[NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN]
const TI = Float64[NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN]
const CW = Float64[0.0, 3.552713678800501e-15, 0.0, -3.552713678800501e-15, 0.0, -4.76837158203125e-07, 0.0, 0.0, 2.384185791015625e-07, 0.0]
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
        x0_idx = Float64[9.777293944358824, 2.9352040820948915, -10.818199560192001, -16.98917378067344, -18.43147217212801, -2634984078.428499, 2548069533.6751747, 1612845257.0755224, 775580670.981751, 596316861.8266063],
        x0_alt_idx = Vector{Float64}[Float64[-15.755691039562228, 4.328701635601931, 24.922373726208, -15.777235814384577, 19.975811573675834, 1316710029.0194902, 1183358454.4442883, -1761492702.9256935, -732690217.8430529, 1477543316.76582], Float64[24.98069525957108, 7.314767821688434, -0.46575764275199916, -20.688773888292072, 4.320056615708872, 697769988.0939012, -181352624.78659916, -399917736.60941315, -2021757410.1912544, 2358769771.7050343], Float64[-0.04308489561081075, 3.9305594774570594, -14.022526820352002, -18.456256581970482, -6.830445116943858, -302056231.86281824, -1831493472.2226388, 1224763503.391173, 1730860416.422398, -1795583515.8655455], Float64[12.759779369831087, 4.793200820104282, 9.492305534976001, -16.223739275648896, 27.071585403545754, -920996272.7884073, 1656101508.0340724, 2849444743.391759, 441793224.0741973, -914357060.926332], Float64[-12.554974973201755, 6.186698373611314, -4.064463642624002, -13.991221969327306, 15.921083670893019, -1539936313.713996, 291390428.803185, -524893216.6094582, -847273968.274004, -33130605.987118244], Float64[18.87023731470108, 2.802490029379939, 19.844747452416, -19.92333938326753, -0.144239762102508, -2539762533.6707163, -1073320650.4277016, 1099788023.3911285, 2905343858.339649, 848095848.9520955], Float64[-6.153542840480807, 5.78855621546645, 6.287978274816002, -17.69082207694594, -11.294741494755234, 1792817752.808404, 2414274329.829009, 2724469263.391713, 1616276665.991447, 1729322303.8913102], Float64[6.649321424961087, 6.474245487827055, -19.100153094144, -15.45830477062435, 22.607289025734374, 1173877711.8828158, 764133482.3929701, -912974970.2938085, 327209473.6432462, 2610548758.8305235], Float64[-18.810920011997226, 3.090037143595673, 16.640420192256002, -20.369842844531846, 11.45678729308164, 174051491.92609596, -600577596.8379173, 711706269.7067766, -961857718.7049553, -1543804528.740056]],
        suite_release = "2.1.0",
        base_function = "ackley",
        m = 0,
        pattern = "C",
        strategy = "extreme",
        scale_factors_C = Float64[1.0, 1.0, 1.0, 1.0, 1.0, 100000000.0, 100000000.0, 100000000.0, 100000000.0, 100000000.0],
        lb_orig = Float64[-19.659594726562503, -30.0, -21.3513984, -21.050229071220325, -19.970767024110373, -30.0, -18.69154344416374, -20.12905364446252, -20.512336648311006, -20.415744802984953],
        ub_orig = Float64[30.0, 18.374272214601433, 30.0, 30.0, 30.0, 19.809936626868804, 30.0, 30.0, 30.0, 30.0],
        lb_C_work = Float64[-19.659594726562503, -30.0, -21.3513984, -21.05022907122033, -19.970767024110373, -3000000000.0, -1869154344.4163742, -2012905364.4462519, -2051233664.8311005, -2041574480.2984953],
        ub_C_work = Float64[30.0, 18.374272214601433, 30.0, 30.0, 30.0, 1980993662.6868804, 3000000000.0, 3000000000.0, 3000000000.0, 3000000000.0],
        ProblemData = _problem_data(),
        func_idx = _objective,
        func_F = xr -> _objective(_cell2idx(xr)),
        idx2cell = _idx2cell,
        cell2idx = _cell2idx,
        decode = xr -> _decode(_cell2idx(xr)))
end

end # module
