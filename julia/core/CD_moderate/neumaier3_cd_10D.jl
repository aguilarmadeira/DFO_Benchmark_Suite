# neumaier3_cd_10D  DFO Benchmark Suite v2.1 core, Julia version (self-contained; Base only).
#NEUMAIER3_CD_10D  DFO Benchmark Suite v2.1.0 core: neumaier3, dimension 10, pattern CD, strategy moderate.
#
# Copyright (C) 2026 J. F. A. Madeira
#
# Objective F(z) = f(T(z)), where f is the original function (embedded below, translated from the MATLAB original) and T acts on the finite
# coordinates only: uniform grid of 17 levels on the interval of the domain, piecewise-affine map with |T - Id| <=
# half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
# Global optimum of F on this domain: f* = -210 (E: exact by structure); the optimal solution is listed in
# reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
# Variables: x1:C x2:D x3:C x4:D x5:C x6:D x7:C x8:D x9:C x10:D
#
# Julia version of the MATLAB wrapper neumaier3_cd_10D.m: same instance, same values (level indices 1..m as in MATLAB).
# Usage:  include("neumaier3_cd_10D.jl");  PB = neumaier3_cd_10D.problem();  f = PB.func_idx(PB.x0_idx)
# Copyright (C) 2026 J. F. A. Madeira. LGPL-3.0-or-later.
module neumaier3_cd_10D

const NAME = "neumaier3_cd_10D"
const N = 10
const TYPES = Char['C', 'D', 'C', 'D', 'C', 'D', 'C', 'D', 'C', 'D']
const CARD = Float64[Inf, 17.0, Inf, 17.0, Inf, 17.0, Inf, 17.0, Inf, 17.0]
const LB = Float64[-100.0, -100.0, -100.0, -72.69804248229904, -66.21159756847209, -71.1692167641189, -100.0, -100.0, -100.0, -100.0]
const UB = Float64[61.813232421875, 69.88954935731343, 72.739328, 100.0, 100.0, 100.0, 62.20825900073035, 69.96665157572455, 61.65233829210158, 65.62204272417893]
const LBW = Float64[-100.0, -100.0, -10000.0, -72.69804248229904, -0.6621159756847208, -71.1692167641189, -100.0, -100.0, -10000.0, -100.0]
const UBW = Float64[61.813232421875, 69.88954935731343, 7273.9328, 100.0, 1.0, 100.0, 62.20825900073035, 69.96665157572455, 6165.233829210158, 65.62204272417893]
const TB = Float64[NaN, 16.799065183152976, NaN, 24.44460641399418, NaN, 25.113467665697996, NaN, 27.474988681793405, NaN, 13.865154372873008]
const TI = Float64[NaN, 18.0, NaN, 28.0, NaN, 30.0, NaN, 24.0, NaN, 10.0]
const CW = Float64[10.0, NaN, 2399.999999999998, NaN, 0.30000000000000016, NaN, 27.999999999999986, NaN, 1800.0, NaN]
const CX = Float64[10.0, NaN, 24.0, NaN, 30.0, NaN, 28.0, NaN, 18.0, NaN]
const M = 17
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
    s1 = 0.0
    for v in x; s1 += (v - 1)^2; end
    s2 = 0.0
    for i in 2:d; s2 += x[i] * x[i-1]; end
    return s1 - s2
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
        f_star = -210.0,
        f_star_justification = "E",
        n_minimizers_documented = 1,
        description = "DFO Benchmark Suite v2.1.0 core (Julia): CD adaptation of neumaier3 (F = f o T).",
        x0_idx = Float64[31.55226171016693, 2.0, -8391.451377664, 15.0, -0.41860537068207576, 11.0, 48.27545108330551, 9.0, -9326.393231576432, 12.0],
        x0_alt_idx = Vector{Float64}[Float64[-48.09018862247467, 2.0, 3631.205851136001, 15.0, 0.7090360463301728, 9.0, 2.8122366228932094, 14.0, 1981.4316324334013, 15.0], Float64[-8.26896345615387, 1.0, -4909.026525184, 16.0, 0.3381506633261442, 7.0, -52.159623999173974, 3.0, -2143.911594302649, 2.0], Float64[-87.91141378879547, 2.0, 7113.630703616, 17.0, -0.1825873996795122, 4.0, -97.62283845958629, 8.0, -6269.254821038701, 5.0], Float64[41.98163020610809, 3.0, -1559.2654766079995, 15.0, -0.5534727826835407, 2.0, 18.560931828134073, 13.0, 5740.077683459584, 7.0], Float64[-37.66082012653351, 2.0, -9436.178833408, 16.0, 0.560545663419469, 16.0, -26.902282632278244, 2.0, 1614.7344567235305, 10.0], Float64[2.1604050397872925, 2.0, 1923.159375872001, 16.0, 0.18966028041544047, 13.0, -72.36549709269056, 8.0, -2510.608770012521, 13.0], Float64[-77.48204529285431, 2.0, -5953.753980928001, 15.0, -0.3310777825902159, 11.0, 34.30962703337494, 12.0, -6635.951996748572, 16.0], Float64[21.59695541858673, 1.0, 5405.584228351998, 16.0, 0.9464163144236606, 7.0, -11.153587427037365, 1.0, 5373.380507749709, 3.0], Float64[-58.04549491405487, 2.0, -2471.329128448, 16.0, 0.425678251418004, 5.0, -56.6168018874497, 6.0, 546.5296405252102, 5.0]],
        suite_release = "2.1.0",
        base_function = "neumaier3",
        m = 17,
        pattern = "CD",
        strategy = "moderate",
        scale_factors_C = Float64[1.0, 100.0, 0.01, 1.0, 100.0],
        lb_orig = Float64[-100.0, -100.0, -100.0, -72.69804248229904, -66.21159756847209, -71.1692167641189, -100.0, -100.0, -100.0, -100.0],
        ub_orig = Float64[61.813232421875, 69.88954935731343, 72.739328, 100.0, 100.0, 100.0, 62.20825900073035, 69.96665157572455, 61.65233829210158, 65.62204272417893],
        lb_C_work = Float64[-100.0, -100.0, -10000.0, -72.69804248229904, -0.6621159756847208, -71.1692167641189, -100.0, -100.0, -10000.0, -100.0],
        ub_C_work = Float64[61.813232421875, 69.88954935731343, 7273.9328, 100.0, 1.0, 100.0, 62.20825900073035, 69.96665157572455, 6165.233829210158, 65.62204272417893],
        ProblemData = _problem_data(),
        func_idx = _objective,
        func_F = xr -> _objective(_cell2idx(xr)),
        idx2cell = _idx2cell,
        cell2idx = _cell2idx,
        decode = xr -> _decode(_cell2idx(xr)))
end

end # module
