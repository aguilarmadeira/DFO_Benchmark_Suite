# fletcher_powell_cdk_3D  DFO Benchmark Suite v2.1 core, Julia version (self-contained; Base only).
#FLETCHER_POWELL_CDK_3D  DFO Benchmark Suite v2.1.0 core: fletcher_powell, dimension 3, pattern CDK, strategy baseline.
#
# Copyright (C) 2026 J. F. A. Madeira
#
# Objective F(z) = f(T(z)), where f is the original function (embedded below, translated from the MATLAB original) and T acts on the finite
# coordinates only: uniform grid of 151 levels on the interval of the domain, piecewise-affine map with |T - Id| <=
# half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
# Global optimum of F on this domain: f* = 0 (E: exact by structure); the optimal solution is listed in
# reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
# Variables: x1:C x2:K x3:D
#
# Julia version of the MATLAB wrapper fletcher_powell_cdk_3D.m: same instance, same values (level indices 1..m as in MATLAB).
# Usage:  include("fletcher_powell_cdk_3D.jl");  PB = fletcher_powell_cdk_3D.problem();  f = PB.func_idx(PB.x0_idx)
# Copyright (C) 2026 J. F. A. Madeira. LGPL-3.0-or-later.
module fletcher_powell_cdk_3D

const NAME = "fletcher_powell_cdk_3D"
const N = 3
const TYPES = Char['C', 'K', 'D']
const CARD = Float64[Inf, 151.0, 151.0]
const LB = Float64[-10.0, -10.0, -10.0]
const UB = Float64[7.1000732421875, 7.351917898694305, 7.173132799999999]
const LBW = Float64[-10.0, -10.0, -10.0]
const UBW = Float64[7.1000732421875, 7.351917898694305, 7.173132799999999]
const TB = Float64[NaN, -0.051567071415265175, -0.039582976000001935]
const TI = Float64[NaN, 0.0, 0.0]
const CW = Float64[1.0, NaN, NaN]
const CX = Float64[1.0, NaN, NaN]
const M = 151
const PERM = Dict{Int,Vector{Int}}(2 => Int[36, 71, 104, 99, 74, 107, 147, 131, 93, 5, 144, 41, 51, 45, 16, 135, 43, 129, 112, 150, 35, 30, 122, 138, 21, 102, 62, 125, 97, 67, 53, 19, 87, 26, 79, 80, 92, 47, 143, 64, 54, 109, 126, 60, 39, 4, 10, 142, 146, 25, 57, 96, 130, 134, 103, 66, 91, 90, 120, 94, 95, 85, 7, 15, 33, 48, 58, 117, 118, 38, 44, 2, 27, 115, 68, 111, 29, 110, 77, 32, 59, 139, 6, 56, 3, 141, 22, 34, 61, 69, 55, 72, 116, 65, 84, 151, 124, 73, 88, 14, 89, 42, 70, 23, 106, 105, 37, 8, 136, 52, 121, 28, 145, 82, 40, 13, 148, 11, 98, 137, 127, 149, 108, 132, 75, 119, 100, 31, 114, 17, 49, 9, 123, 20, 1, 50, 18, 128, 113, 24, 83, 12, 101, 63, 78, 81, 140, 86, 76, 133, 46])

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
    x[1] == 0 && return Inf
    teta = x[1] > 0 ? 1 / (2 * pi) * atan(x[2] / x[1]) : 1 / (2 * pi) * atan(x[2] / x[1]) + 0.5
    nrm = sqrt(x[1]^2 + x[2]^2 + x[3]^2)
    return 100 * ((x[3] - 10 * teta)^2 + (nrm - 1)^2) + x[3]^2
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
        cardinality = Float64[Inf, 151.0, 151.0],
        f_star = 0.0,
        f_star_justification = "E",
        n_minimizers_documented = 1,
        description = "DFO Benchmark Suite v2.1.0 core (Julia): CDK adaptation of fletcher_powell (F = f o T).",
        x0_idx = Float64[-5.741680979728699, 31.0, 26.0],
        x0_alt_idx = Vector{Float64}[Float64[4.703725087642669, 20.0, 101.0], Float64[-3.712717211246491, 25.0, 56.0], Float64[0.49550393819808924, 34.0, 132.0], Float64[-7.920938360691071, 22.0, 85.0], Float64[5.8058782458305345, 26.0, 17.0], Float64[-2.610564053058625, 30.0, 115.0], Float64[1.597657096385955, 18.0, 47.0], Float64[-6.8187852025032045, 27.0, 146.0], Float64[3.651669800281523, 31.0, 77.0]],
        suite_release = "2.1.0",
        base_function = "fletcher_powell",
        m = 151,
        pattern = "CDK",
        strategy = "baseline",
        scale_factors_C = Float64[1.0],
        lb_orig = Float64[-10.0, -10.0, -10.0],
        ub_orig = Float64[7.1000732421875, 7.351917898694305, 7.173132799999999],
        lb_C_work = Float64[-10.0, -10.0, -10.0],
        ub_C_work = Float64[7.1000732421875, 7.351917898694305, 7.173132799999999],
        ProblemData = _problem_data(),
        func_idx = _objective,
        func_F = xr -> _objective(_cell2idx(xr)),
        idx2cell = _idx2cell,
        cell2idx = _cell2idx,
        decode = xr -> _decode(_cell2idx(xr)))
end

end # module
