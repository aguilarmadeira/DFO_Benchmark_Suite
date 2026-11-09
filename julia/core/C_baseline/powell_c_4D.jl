# powell_c_4D  DFO Benchmark Suite v2.1 core, Julia version (self-contained; Base only).
#POWELL_C_4D  DFO Benchmark Suite v2.1.0 core: powell, dimension 4, pattern C, strategy baseline.
#
# Copyright (C) 2026 J. F. A. Madeira
#
# Objective F(z) = f(T(z)), where f is the original function (embedded below, translated from the MATLAB original) and T acts on the finite
# coordinates only: uniform grid of 0 levels on the interval of the domain, piecewise-affine map with |T - Id| <=
# half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
# Global optimum of F on this domain: f* = 0 (E: exact by structure); the optimal solution is listed in
# reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
# Variables: x1:C x2:C x3:C x4:C
#
# Julia version of the MATLAB wrapper powell_c_4D.m: same instance, same values (level indices 1..m as in MATLAB).
# Usage:  include("powell_c_4D.jl");  PB = powell_c_4D.problem();  f = PB.func_idx(PB.x0_idx)
# Copyright (C) 2026 J. F. A. Madeira. LGPL-3.0-or-later.
module powell_c_4D

const NAME = "powell_c_4D"
const N = 4
const TYPES = Char['C', 'C', 'C', 'C']
const CARD = Float64[Inf, Inf, Inf, Inf]
const LB = Float64[-6.0063232421875, -6.781547528323934, -10.0, -10.0]
const UB = Float64[10.0, 10.0, 6.1539328, 6.469804248229904]
const LBW = Float64[-6.0063232421875, -6.781547528323934, -10.0, -10.0]
const UBW = Float64[10.0, 10.0, 6.1539328, 6.469804248229904]
const TB = Float64[NaN, NaN, NaN, NaN]
const TI = Float64[NaN, NaN, NaN, NaN]
const CW = Float64[0.0, 0.0, 0.0, 0.0]
const CX = Float64[0.0, 0.0, 0.0, 0.0]
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

_orig(x::Vector{Float64}) = (x[1] + 10 * x[2])^2 + 5 * (x[3] - x[4])^2 + (x[2] - 2 * x[3])^4 + 10 * (x[1] - x[4])^4

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
        cardinality = Float64[Inf, Inf, Inf, Inf],
        f_star = 0.0,
        f_star_justification = "E",
        n_minimizers_documented = 1,
        description = "DFO Benchmark Suite v2.1.0 core (Julia): C adaptation of powell (F = f o T).",
        x0_idx = Float64[9.00742038488388, 1.7435093093385685, 4.427400462336001, -2.5848153301388885],
        x0_alt_idx = Vector{Float64}[Float64[1.1293081641197196, 2.22692837394049, -2.938792894463999, -1.864561500041372], Float64[5.0683642745018, 1.0529106456215338, -7.823742173184, -1.1443076699438546], Float64[-2.8097479462623602, 1.3521700665655825, 0.31783995801600007, -2.7288660961583924], Float64[7.0378923296928395, 2.3880680621411337, -4.691171524608, -2.008612266060876], Float64[-0.8402198910713201, 1.2140503338221746, 3.4504106065920013, -1.2883584359633549], Float64[3.09883621931076, 1.697469398424099, -1.4345386721279993, -0.8973634996247064], Float64[-4.7792760014534, 1.9353422714821882, -8.800732028928, -2.4819219258392433], Float64[8.022656357288358, 0.7613245431632292, 1.8220941803519999, -1.761668095741726], Float64[0.14454413652419973, 1.797222538738783, -5.544099176448, -1.0414142656442102]],
        suite_release = "2.1.0",
        base_function = "powell",
        m = 0,
        pattern = "C",
        strategy = "baseline",
        scale_factors_C = Float64[1.0, 1.0, 1.0, 1.0],
        lb_orig = Float64[-6.0063232421875, -6.781547528323934, -10.0, -10.0],
        ub_orig = Float64[10.0, 10.0, 6.1539328, 6.469804248229904],
        lb_C_work = Float64[-6.0063232421875, -6.781547528323934, -10.0, -10.0],
        ub_C_work = Float64[10.0, 10.0, 6.1539328, 6.469804248229904],
        ProblemData = _problem_data(),
        func_idx = _objective,
        func_F = xr -> _objective(_cell2idx(xr)),
        idx2cell = _idx2cell,
        cell2idx = _cell2idx,
        decode = xr -> _decode(_cell2idx(xr)))
end

end # module
