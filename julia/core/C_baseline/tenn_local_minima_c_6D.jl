# tenn_local_minima_c_6D  DFO Benchmark Suite v2.1 core, Julia version (self-contained; Base only).
#TENN_LOCAL_MINIMA_C_6D  DFO Benchmark Suite v2.1.0 core: tenn_local_minima, dimension 6, pattern C, strategy baseline.
#
# Copyright (C) 2026 J. F. A. Madeira
#
# Objective F(z) = f(T(z)), where f is the original function (embedded below, translated from the MATLAB original) and T acts on the finite
# coordinates only: uniform grid of 0 levels on the interval of the domain, piecewise-affine map with |T - Id| <=
# half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
# Global optimum of F on this domain: f* = 7.852723858933194e-32 (E: exact by structure); the optimal solution is listed in
# reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
# Variables: x1:C x2:C x3:C x4:C x5:C x6:C
#
# Julia version of the MATLAB wrapper tenn_local_minima_c_6D.m: same instance, same values (level indices 1..m as in MATLAB).
# Usage:  include("tenn_local_minima_c_6D.jl");  PB = tenn_local_minima_c_6D.problem();  f = PB.func_idx(PB.x0_idx)
# Copyright (C) 2026 J. F. A. Madeira. LGPL-3.0-or-later.
module tenn_local_minima_c_6D

const NAME = "tenn_local_minima_c_6D"
const N = 6
const TYPES = Char['C', 'C', 'C', 'C', 'C', 'C']
const CARD = Float64[Inf, Inf, Inf, Inf, Inf, Inf]
const LB = Float64[-10.0, -10.0, -10.0, -10.0, -10.0, -10.0]
const UB = Float64[6.2961669921875, 7.188680587308845, 6.1987328, 6.780008329862557, 6.343473806433986, 7.315738244459228]
const LBW = Float64[-10.0, -10.0, -10.0, -10.0, -10.0, -10.0]
const UBW = Float64[6.2961669921875, 7.188680587308845, 6.1987328, 6.780008329862557, 6.343473806433986, 7.315738244459228]
const TB = Float64[NaN, NaN, NaN, NaN, NaN, NaN]
const TI = Float64[NaN, NaN, NaN, NaN, NaN, NaN]
const CW = Float64[1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
const CX = Float64[1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
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
        cardinality = Float64[Inf, Inf, Inf, Inf, Inf, Inf],
        f_star = 7.852723858933194e-32,
        f_star_justification = "E",
        n_minimizers_documented = 1,
        description = "DFO Benchmark Suite v2.1.0 core (Julia): C adaptation of tenn_local_minima (F = f o T).",
        x0_idx = Float64[2.4847099661827095, -6.078119975735202, 3.9490527887360027, -3.8149490329327937, -4.363895960065215, -1.653451615438179],
        x0_alt_idx = Vector{Float64}[Float64[-5.536059725284576, -5.01709030985194, -0.3274126704639979, -3.081129426670583, -9.484278061705314, -3.805111397294059], Float64[-1.5256748795509338, -4.52194313243975, -8.336066166784, -4.695532560447447, 3.077234863900971, -7.280869506445866], Float64[-9.546444571018219, -5.913070916597806, 2.9382518620160027, -3.961712954185236, -2.0431472377391264, -9.432529288301746], Float64[5.826697337627412, -5.417923739185617, -5.194807902207999, -3.2278933479230254, -5.690038087108694, 5.629089184689418], Float64[-2.194072353839873, -4.356894073302354, 6.079510126592002, -5.177756873134043, 5.264029836020649, 2.153331075537608], Float64[1.81631249189377, -5.559394361303384, -1.9291433697279992, -4.443937266871832, 1.6171389866510797, 0.0016712936817295798], Float64[-6.2044571995735165, -5.315750512100562, -9.937796866048, -3.7101176606096216, -3.5032431149890186, -2.1499884881741513], Float64[3.7737622380256663, -4.820603334688373, 1.3365211627520015, -2.976298054347409, -7.150133964358584, -5.625746597325957], Float64[-4.247007453441619, -6.023103622689403, -6.672132333567999, -4.590701188124274, 3.9378877089771684, -7.777406379181838]],
        suite_release = "2.1.0",
        base_function = "tenn_local_minima",
        m = 0,
        pattern = "C",
        strategy = "baseline",
        scale_factors_C = Float64[1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
        lb_orig = Float64[-10.0, -10.0, -10.0, -10.0, -10.0, -10.0],
        ub_orig = Float64[6.2961669921875, 7.188680587308845, 6.1987328, 6.780008329862557, 6.343473806433986, 7.315738244459228],
        lb_C_work = Float64[-10.0, -10.0, -10.0, -10.0, -10.0, -10.0],
        ub_C_work = Float64[6.2961669921875, 7.188680587308845, 6.1987328, 6.780008329862557, 6.343473806433986, 7.315738244459228],
        ProblemData = _problem_data(),
        func_idx = _objective,
        func_F = xr -> _objective(_cell2idx(xr)),
        idx2cell = _idx2cell,
        cell2idx = _cell2idx,
        decode = xr -> _decode(_cell2idx(xr)))
end

end # module
