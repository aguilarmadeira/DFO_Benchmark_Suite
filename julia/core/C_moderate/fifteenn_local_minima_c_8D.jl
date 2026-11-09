# fifteenn_local_minima_c_8D  DFO Benchmark Suite v2.1 core, Julia version (self-contained; Base only).
#FIFTEENN_LOCAL_MINIMA_C_8D  DFO Benchmark Suite v2.1.0 core: fifteenn_local_minima, dimension 8, pattern C, strategy moderate.
#
# Copyright (C) 2026 J. F. A. Madeira
#
# Objective F(z) = f(T(z)), where f is the original function (embedded below, translated from the MATLAB original) and T acts on the finite
# coordinates only: uniform grid of 0 levels on the interval of the domain, piecewise-affine map with |T - Id| <=
# half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
# Global optimum of F on this domain: f* = 1.349783804395672e-32 (E: exact by structure); the optimal solution is listed in
# reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
# Variables: x1:C x2:C x3:C x4:C x5:C x6:C x7:C x8:C
#
# Julia version of the MATLAB wrapper fifteenn_local_minima_c_8D.m: same instance, same values (level indices 1..m as in MATLAB).
# Usage:  include("fifteenn_local_minima_c_8D.jl");  PB = fifteenn_local_minima_c_8D.problem();  f = PB.func_idx(PB.x0_idx)
# Copyright (C) 2026 J. F. A. Madeira. LGPL-3.0-or-later.
module fifteenn_local_minima_c_8D

const NAME = "fifteenn_local_minima_c_8D"
const N = 8
const TYPES = Char['C', 'C', 'C', 'C', 'C', 'C', 'C', 'C']
const CARD = Float64[Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf]
const LB = Float64[-10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0]
const UB = Float64[6.9250732421875, 6.487720367830107, 6.3891328, 6.55551853394419, 7.396366368417458, 6.164258954518399, 6.225670190730475, 6.853174852863314]
const LBW = Float64[-10.0, -1000.0, -0.1, -10.0, -1000.0, -0.1, -10.0, -1000.0]
const UBW = Float64[6.9250732421875005, 648.7720367830107, 0.063891328, 6.555518533944191, 739.6366368417457, 0.061642589545183986, 6.225670190730476, 685.3174852863315]
const TB = Float64[NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN]
const TI = Float64[NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN]
const CW = Float64[1.0, 100.0, 0.009999999999999995, 1.0, 100.0, 0.009999999999999995, 1.0, 100.0]
const CX = Float64[1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
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
    for i in 1:d-1; s += ((x[i] - 1)^2) * (1 + 10 * sin(3 * pi * x[i+1])^2); end
    f = 0.1 * s
    return f + 0.1 * (sin(3 * pi * x[1])^2 + (x[d] - 1)^2 * (1 + sin(2 * pi * x[d])^2))
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
        cardinality = Float64[Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf],
        f_star = 1.3497838043956716e-32,
        f_star_justification = "E",
        n_minimizers_documented = 1,
        description = "DFO Benchmark Suite v2.1.0 core (Julia): C adaptation of fifteenn_local_minima (F = f o T).",
        x0_idx = Float64[-3.669625926017761, -304.90726204209614, 0.01884087975936001, -0.9120477185595828, 16.85747818398056, -0.04982237320445359, 4.775829113235934, 46.228583226301225],
        x0_alt_idx = Vector{Float64}[Float64[0.3219807028770454, -420.2534539158047, -0.05589356580864, -2.8358251741907488, -371.32590359062385, -0.06990813876923976, 0.22814992483050744, 592.439513360652], Float64[-8.008328783512116, -372.7579631442776, 0.05188137148416, -2.1118229059424607, -916.3510557792098, -0.08999390433402593, -4.31952926357492, -630.4537836607899], Float64[5.6771796584129355, -270.98191149100535, -0.022853074083840005, -1.3878206376941709, 420.7250369999831, 0.03824598350268582, -9.818357040405012, -84.2428535264387], Float64[-2.653129827976226, -404.42162365862896, -0.07367249707008, -0.6638183694458846, -124.30011518860272, 0.01816021793789961, 1.8034897744088578, 461.9680766079125], Float64[1.412854564189912, -356.92613288710186, 0.00892873224192002, -2.256623359592118, -683.5836101420962, -0.0019255476268865285, -2.74418941399657, -672.47000905574], Float64[-6.91745492219925, -309.43064211557476, -0.04063200534528001, -1.8635935568287607, 653.4924826370968, -0.03437178430846423, -7.291868602401997, -126.25907892138878], Float64[3.594602286815645, -448.9015277145036, 0.04196922396672001, -1.1395912885804744, 108.46733044851067, -0.0544575498732504, 3.37882962398721, 331.4966398551728], Float64[-4.7357071995735165, -347.1254760612313, -0.007591513620480009, -2.732396278726706, -279.71605132609375, -0.08690378655482807, -1.1688495644182186, -802.9414458084797], Float64[-0.7193079829215989, -299.62998528970434, -0.08232595918848, -2.008394010478419, -824.7412035146797, 0.053696572398675174, -5.716528752823646, -256.7305156741286]],
        suite_release = "2.1.0",
        base_function = "fifteenn_local_minima",
        m = 0,
        pattern = "C",
        strategy = "moderate",
        scale_factors_C = Float64[1.0, 100.0, 0.01, 1.0, 100.0, 0.01, 1.0, 100.0],
        lb_orig = Float64[-10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0],
        ub_orig = Float64[6.9250732421875, 6.487720367830107, 6.3891328, 6.55551853394419, 7.396366368417458, 6.164258954518399, 6.225670190730475, 6.853174852863314],
        lb_C_work = Float64[-10.0, -1000.0, -0.1, -10.0, -1000.0, -0.1, -10.0, -1000.0],
        ub_C_work = Float64[6.9250732421875005, 648.7720367830107, 0.063891328, 6.555518533944191, 739.6366368417457, 0.061642589545183986, 6.225670190730476, 685.3174852863315],
        ProblemData = _problem_data(),
        func_idx = _objective,
        func_F = xr -> _objective(_cell2idx(xr)),
        idx2cell = _idx2cell,
        cell2idx = _cell2idx,
        decode = xr -> _decode(_cell2idx(xr)))
end

end # module
