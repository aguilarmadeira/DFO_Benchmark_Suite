# sphere_cdk_3D  DFO Benchmark Suite v2.1 core, Julia version (self-contained; Base only).
#SPHERE_CDK_3D  DFO Benchmark Suite v2.1.0 core: sphere, dimension 3, pattern CDK, strategy baseline.
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
# Julia version of the MATLAB wrapper sphere_cdk_3D.m: same instance, same values (level indices 1..m as in MATLAB).
# Usage:  include("sphere_cdk_3D.jl");  PB = sphere_cdk_3D.problem();  f = PB.func_idx(PB.x0_idx)
# Copyright (C) 2026 J. F. A. Madeira. LGPL-3.0-or-later.
module sphere_cdk_3D

const NAME = "sphere_cdk_3D"
const N = 3
const TYPES = Char['C', 'K', 'D']
const CARD = Float64[Inf, 151.0, 151.0]
const LB = Float64[-3.1116375, -5.12, -3.7471911936]
const UB = Float64[5.12, 3.149641497739166, 5.12]
const LBW = Float64[-3.1116375, -5.12, -3.7471911936]
const UBW = Float64[5.12, 3.149641497739166, 5.12]
const TB = Float64[NaN, 0.0071777285982825845, -0.022970892287999867]
const TI = Float64[NaN, 0.0, 0.0]
const CW = Float64[0.0, NaN, NaN]
const CX = Float64[0.0, NaN, NaN]
const M = 151
const PERM = Dict{Int,Vector{Int}}(2 => Int[2, 17, 50, 49, 147, 14, 35, 78, 130, 75, 33, 134, 148, 142, 90, 83, 73, 53, 85, 19, 82, 41, 67, 145, 93, 39, 88, 151, 9, 12, 101, 84, 45, 1, 135, 139, 136, 122, 31, 36, 132, 51, 13, 150, 29, 57, 74, 146, 91, 80, 54, 87, 64, 95, 20, 124, 44, 6, 16, 38, 11, 23, 127, 70, 106, 77, 40, 63, 113, 43, 137, 129, 140, 7, 104, 123, 115, 60, 141, 138, 28, 56, 125, 107, 109, 58, 120, 27, 66, 32, 55, 133, 111, 144, 128, 37, 30, 26, 97, 34, 46, 24, 65, 47, 10, 5, 114, 119, 110, 25, 96, 98, 3, 62, 112, 126, 118, 61, 72, 8, 68, 81, 4, 100, 149, 69, 99, 143, 102, 15, 117, 116, 79, 108, 76, 59, 52, 121, 86, 131, 18, 22, 48, 21, 42, 105, 89, 92, 71, 94, 103])

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
    return s
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
        description = "DFO Benchmark Suite v2.1.0 core (Julia): CDK adaptation of sphere (F = f o T).",
        x0_idx = Float64[4.866780682373047, 102.0, 10.0],
        x0_alt_idx = Vector{Float64}[Float64[0.8152716003417968, 111.0, 115.0], Float64[2.8410261413574216, 115.0, 40.0], Float64[-1.2104829406738282, 103.0, 145.0], Float64[3.8539034118652338, 108.0, 69.0], Float64[-0.19760567016601582, 117.0, 30.0], Float64[1.8281488708496094, 106.0, 100.0], Float64[-2.2233602111816406, 109.0, 60.0], Float64[4.36034204711914, 113.0, 130.0], Float64[0.30883296508789027, 102.0, 90.0]],
        suite_release = "2.1.0",
        base_function = "sphere",
        m = 151,
        pattern = "CDK",
        strategy = "baseline",
        scale_factors_C = Float64[1.0],
        lb_orig = Float64[-3.1116375, -5.12, -3.7471911936],
        ub_orig = Float64[5.12, 3.149641497739166, 5.12],
        lb_C_work = Float64[-3.1116375, -5.12, -3.7471911936],
        ub_C_work = Float64[5.12, 3.149641497739166, 5.12],
        ProblemData = _problem_data(),
        func_idx = _objective,
        func_F = xr -> _objective(_cell2idx(xr)),
        idx2cell = _idx2cell,
        cell2idx = _cell2idx,
        decode = xr -> _decode(_cell2idx(xr)))
end

end # module
