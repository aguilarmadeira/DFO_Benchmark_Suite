# salomon_c_10D  DFO Benchmark Suite v2.1 core, Julia version (self-contained; Base only).
#SALOMON_C_10D  DFO Benchmark Suite v2.1.0 core: salomon, dimension 10, pattern C, strategy moderate.
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
# Julia version of the MATLAB wrapper salomon_c_10D.m: same instance, same values (level indices 1..m as in MATLAB).
# Usage:  include("salomon_c_10D.jl");  PB = salomon_c_10D.problem();  f = PB.func_idx(PB.x0_idx)
# Copyright (C) 2026 J. F. A. Madeira. LGPL-3.0-or-later.
module salomon_c_10D

const NAME = "salomon_c_10D"
const N = 10
const TYPES = Char['C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C']
const CARD = Float64[Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf]
const LB = Float64[-100.0, -100.0, -63.779328, -70.41232819658475, -69.91407690731508, -62.55383214873429, -100.0, -100.0, -100.0, -100.0]
const UB = Float64[63.399169921875, 61.59325306101712, 100.0, 100.0, 100.0, 100.0, 69.5715808000383, 64.0719147336193, 70.75630804635489, 61.75997375866169]
const LBW = Float64[-100.0, -10000.0, -0.6377932799999999, -70.41232819658475, -6991.4076907315075, -0.6255383214873428, -100.0, -9999.999999999998, -1.0, -100.0]
const UBW = Float64[63.399169921875, 6159.325306101713, 1.0, 100.0, 10000.0, 0.9999999999999999, 69.5715808000383, 6407.191473361929, 0.7075630804635489, 61.75997375866169]
const TB = Float64[NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN]
const TI = Float64[NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN]
const CW = Float64[0.0, 1.8189894035458565e-12, 1.1102230246251565e-16, 0.0, 0.0, 1.1102230246251565e-16, 0.0, 1.8189894035458565e-12, -1.1102230246251565e-16, 0.0]
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
    for v in x; s += v^2; end
    r = sqrt(s)
    return 1 - cos(2 * pi * r) + 0.1 * r
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
        description = "DFO Benchmark Suite v2.1.0 core (Julia): C adaptation of salomon (F = f o T).",
        x0_idx = Float64[27.73538625240326, 5390.889169003138, 0.5629057294336002, -63.243796273279855, 7855.329457518487, -0.2267376135757917, -34.18013340409669, 5412.091363591033, -0.3349107061463995, 57.34273063580841],
        x0_alt_idx = Vector{Float64}[Float64[-52.68764269351959, 5856.384136668619, -0.18392800624639988, -55.79136209558664, 2531.950789573303, -0.4287275825180059, -81.70711625808664, -6493.228940086223, -0.8447797512128967, -75.95056522002102], Float64[-12.476128220558167, 4725.896358052452, 0.8930848546815999, -48.338927917893436, -1259.5203480279415, 0.9852022000774935, 29.81044481761532, -1175.6627284980077, 0.4237879059178682, -47.676229735451145], Float64[-92.89915716648102, 5014.059909464417, 0.14625111900160004, -64.7342831088185, -6582.899015973124, 0.6589107117862244, -17.716538036374658, 4141.903483090207, -0.011979609890409915, -19.401894250881284], Float64[53.58564555644989, 6011.549125890449, -0.3615958212607999, -57.28184893112528, 6476.612680208945, 0.4569207428440102, -65.24352089036464, -6902.272494823778, -0.4477471256986878, 8.872441233688605], Float64[-26.83738338947296, 4881.061347274279, 0.46385199185920034, -49.829414753432076, 1013.969691323402, 0.1306292545527411, 56.214324180943095, -1584.7062832355623, -0.8835146415069659, 37.14677671825851], Float64[13.374131083488464, 5346.55631493976, -0.03141669601280006, -69.6315969970169, -2777.5014462778427, -0.07136071438947311, 8.687341326953117, 3732.859928352653, 0.385053015623799, -96.14651913757093], Float64[-67.04889786243439, 5575.60939426722, 0.7940311171072003, -62.17916281932368, 8750.102719560287, -0.2733506833316872, -48.779925522642614, -8172.4603753246065, -0.050714500184479094, -67.87218365300106], Float64[33.001179814338684, 4445.121615651051, 0.2987624292352, -54.72672864163047, 4958.631581959041, -0.5996421716229563, -96.30690837663259, -2854.8941637363932, -0.486482015992757, -39.59784816843118], Float64[-47.42184913158417, 5442.610832077082, -0.4480713064447999, -47.27429446393726, -364.7470859861414, 0.8142876109725431, 25.150936694675124, 2462.672047851822, -0.9963510610592543, -11.323512683861296]],
        suite_release = "2.1.0",
        base_function = "salomon",
        m = 0,
        pattern = "C",
        strategy = "moderate",
        scale_factors_C = Float64[1.0, 100.0, 0.01, 1.0, 100.0, 0.01, 1.0, 100.0, 0.01, 1.0],
        lb_orig = Float64[-100.0, -100.0, -63.779328, -70.41232819658475, -69.91407690731508, -62.55383214873429, -100.0, -100.0, -100.0, -100.0],
        ub_orig = Float64[63.399169921875, 61.59325306101712, 100.0, 100.0, 100.0, 100.0, 69.5715808000383, 64.0719147336193, 70.75630804635489, 61.75997375866169],
        lb_C_work = Float64[-100.0, -10000.0, -0.6377932799999999, -70.41232819658475, -6991.4076907315075, -0.6255383214873428, -100.0, -9999.999999999998, -1.0, -100.0],
        ub_C_work = Float64[63.399169921875, 6159.325306101713, 1.0, 100.0, 10000.0, 0.9999999999999999, 69.5715808000383, 6407.191473361929, 0.7075630804635489, 61.75997375866169],
        ProblemData = _problem_data(),
        func_idx = _objective,
        func_F = xr -> _objective(_cell2idx(xr)),
        idx2cell = _idx2cell,
        cell2idx = _cell2idx,
        decode = xr -> _decode(_cell2idx(xr)))
end

end # module
