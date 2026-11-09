# neumaier3_cdk_10D  DFO Benchmark Suite v2.1 core, Julia version (self-contained; Base only).
#NEUMAIER3_CDK_10D  DFO Benchmark Suite v2.1.0 core: neumaier3, dimension 10, pattern CDK, strategy extreme.
#
# Copyright (C) 2026 J. F. A. Madeira
#
# Objective F(z) = f(T(z)), where f is the original function (embedded below, translated from the MATLAB original) and T acts on the finite
# coordinates only: uniform grid of 17 levels on the interval of the domain, piecewise-affine map with |T - Id| <=
# half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
# Global optimum of F on this domain: f* = -210 (E: exact by structure); the optimal solution is listed in
# reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
# Variables: x1:C x2:K x3:D x4:C x5:D x6:K x7:C x8:K x9:D x10:C
#
# Julia version of the MATLAB wrapper neumaier3_cdk_10D.m: same instance, same values (level indices 1..m as in MATLAB).
# Usage:  include("neumaier3_cdk_10D.jl");  PB = neumaier3_cdk_10D.problem();  f = PB.func_idx(PB.x0_idx)
# Copyright (C) 2026 J. F. A. Madeira. LGPL-3.0-or-later.
module neumaier3_cdk_10D

const NAME = "neumaier3_cdk_10D"
const N = 10
const TYPES = Char['C', 'K', 'D', 'C', 'D', 'K', 'C', 'K', 'D', 'C']
const CARD = Float64[Inf, 17.0, 17.0, Inf, 17.0, 17.0, Inf, 17.0, 17.0, Inf]
const LB = Float64[-100.0, -100.0, -100.0, -72.69804248229904, -66.21159756847209, -71.1692167641189, -100.0, -100.0, -100.0, -100.0]
const UB = Float64[61.813232421875, 69.88954935731343, 72.739328, 100.0, 100.0, 100.0, 62.20825900073035, 69.96665157572455, 61.65233829210158, 65.62204272417893]
const LBW = Float64[-100.0, -100.0, -100.0, -72.69804248229906, -66.21159756847209, -71.1692167641189, -10000000000.0, -100.0, -100.0, -10000000000.0]
const UBW = Float64[61.813232421875, 69.88954935731343, 72.739328, 100.0, 100.0, 100.0, 6220825900.073034, 69.96665157572455, 61.65233829210158, 6562204272.417892]
const TB = Float64[NaN, 16.799065183152976, 18.758287999999993, NaN, 27.282426063793466, 25.113467665697996, NaN, 27.474988681793405, 21.239253719076174, NaN]
const TI = Float64[NaN, 18.0, 24.0, NaN, 30.0, 30.0, NaN, 24.0, 18.0, NaN]
const CW = Float64[10.0, NaN, NaN, 27.999999999999986, NaN, NaN, 2799999999.999998, NaN, NaN, 1000000000.0]
const CX = Float64[10.0, NaN, NaN, 28.0, NaN, NaN, 28.0, NaN, NaN, 10.0]
const M = 17
const PERM = Dict{Int,Vector{Int}}(2 => Int[6, 5, 10, 12, 14, 3, 16, 11, 17, 1, 2, 7, 4, 15, 8, 13, 9], 6 => Int[3, 6, 1, 10, 2, 15, 14, 16, 9, 13, 8, 11, 4, 17, 5, 7, 12], 8 => Int[16, 15, 4, 6, 12, 2, 9, 10, 8, 5, 1, 7, 3, 13, 11, 17, 14])

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
        cardinality = Float64[Inf, 17.0, 17.0, Inf, 17.0, 17.0, Inf, 17.0, 17.0, Inf],
        f_star = -210.0,
        f_star_justification = "E",
        n_minimizers_documented = 1,
        description = "DFO Benchmark Suite v2.1.0 core (Julia): CDK adaptation of neumaier3 (F = f o T).",
        x0_idx = Float64[31.55226171016693, 2.0, 2.0, 75.54463371762529, 3.0, 11.0, 4827545108.330551, 9.0, 2.0, 1541049719.5351238],
        x0_alt_idx = Vector{Float64}[Float64[-48.09018862247467, 2.0, 14.0, 83.09702624600573, 14.0, 9.0, 281223662.2893219, 14.0, 13.0, 4435989110.790092], Float64[-8.26896345615387, 1.0, 6.0, 90.64941877438613, 11.0, 7.0, -5215962399.9173975, 3.0, 9.0, -9211582305.126198], Float64[-87.91141378879547, 2.0, 17.0, 98.20181130276657, 6.0, 4.0, -9762283845.958628, 8.0, 5.0, -6316642913.871227], Float64[41.98163020610809, 3.0, 9.0, 81.58654774032965, 2.0, 2.0, 1856093182.813406, 13.0, 17.0, -3421703522.6162567], Float64[-37.66082012653351, 2.0, 2.0, 89.13894026871006, 13.0, 16.0, -2690228263.227825, 2.0, 12.0, -526764131.36128426], Float64[2.1604050397872925, 2.0, 12.0, 93.23881049840227, 9.0, 13.0, -7236549709.269056, 8.0, 8.0, 2368175259.8936863], Float64[-77.48204529285431, 2.0, 5.0, 76.62354693596532, 4.0, 11.0, 3430962703.337494, 12.0, 4.0, 5263114651.148655], Float64[21.59695541858673, 1.0, 15.0, 84.17593946434576, 16.0, 7.0, -1115358742.7037373, 1.0, 16.0, -8384456764.767634], Float64[-58.04549491405487, 2.0, 8.0, 91.72833199272617, 11.0, 5.0, -5661680188.744969, 6.0, 11.0, -5489517373.512664]],
        suite_release = "2.1.0",
        base_function = "neumaier3",
        m = 17,
        pattern = "CDK",
        strategy = "extreme",
        scale_factors_C = Float64[1.0, 1.0, 100000000.0, 100000000.0],
        lb_orig = Float64[-100.0, -100.0, -100.0, -72.69804248229904, -66.21159756847209, -71.1692167641189, -100.0, -100.0, -100.0, -100.0],
        ub_orig = Float64[61.813232421875, 69.88954935731343, 72.739328, 100.0, 100.0, 100.0, 62.20825900073035, 69.96665157572455, 61.65233829210158, 65.62204272417893],
        lb_C_work = Float64[-100.0, -100.0, -100.0, -72.69804248229906, -66.21159756847209, -71.1692167641189, -10000000000.0, -100.0, -100.0, -10000000000.0],
        ub_C_work = Float64[61.813232421875, 69.88954935731343, 72.739328, 100.0, 100.0, 100.0, 6220825900.073034, 69.96665157572455, 61.65233829210158, 6562204272.417892],
        ProblemData = _problem_data(),
        func_idx = _objective,
        func_F = xr -> _objective(_cell2idx(xr)),
        idx2cell = _idx2cell,
        cell2idx = _cell2idx,
        decode = xr -> _decode(_cell2idx(xr)))
end

end # module
