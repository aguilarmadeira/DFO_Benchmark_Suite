# wood_cdk_4D  DFO Benchmark Suite v2.1 core, Julia version (self-contained; Base only).
#WOOD_CDK_4D  DFO Benchmark Suite v2.1.0 core: wood, dimension 4, pattern CDK, strategy extreme.
#
# Copyright (C) 2026 J. F. A. Madeira
#
# Objective F(z) = f(T(z)), where f is the original function (embedded below, translated from the MATLAB original) and T acts on the finite
# coordinates only: uniform grid of 151 levels on the interval of the domain, piecewise-affine map with |T - Id| <=
# half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
# Global optimum of F on this domain: f* = 0 (E: exact by structure); the optimal solution is listed in
# reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
# Variables: x1:C x2:K x3:D x4:C
#
# Julia version of the MATLAB wrapper wood_cdk_4D.m: same instance, same values (level indices 1..m as in MATLAB).
# Usage:  include("wood_cdk_4D.jl");  PB = wood_cdk_4D.problem();  f = PB.func_idx(PB.x0_idx)
# Copyright (C) 2026 J. F. A. Madeira. LGPL-3.0-or-later.
module wood_cdk_4D

const NAME = "wood_cdk_4D"
const N = 4
const TYPES = Char['C', 'K', 'D', 'C']
const CARD = Float64[Inf, 151.0, 151.0, Inf]
const LB = Float64[-10.0, -10.0, -10.0, -10.0]
const UB = Float64[6.1211669921875, 6.566458365086623, 6.4227327999999995, 7.3514369012911285]
const LBW = Float64[-10.0, -10.0, -10.0, -1000000000.0]
const UBW = Float64[6.121166992187501, 6.566458365086623, 6.4227327999999995, 735143690.129113]
const TB = Float64[NaN, 1.0443055767244172, 0.9484885333333324, NaN]
const TI = Float64[NaN, 1.0, 1.0, NaN]
const CW = Float64[1.0, NaN, NaN, 100000000.0]
const CX = Float64[1.0, NaN, NaN, 1.0]
const M = 151
const PERM = Dict{Int,Vector{Int}}(2 => Int[30, 135, 36, 60, 80, 141, 142, 61, 111, 151, 55, 102, 78, 62, 148, 43, 45, 70, 53, 116, 33, 103, 34, 118, 100, 42, 82, 146, 128, 17, 115, 139, 39, 73, 86, 88, 108, 68, 2, 83, 46, 145, 26, 98, 105, 147, 64, 28, 38, 87, 110, 91, 72, 127, 48, 94, 93, 112, 51, 15, 130, 3, 10, 89, 52, 150, 25, 49, 92, 96, 41, 16, 132, 84, 47, 6, 138, 131, 27, 1, 44, 74, 124, 4, 54, 125, 50, 134, 24, 90, 31, 58, 123, 8, 37, 104, 109, 114, 20, 59, 79, 63, 5, 149, 121, 101, 56, 75, 137, 57, 22, 40, 76, 35, 106, 117, 85, 12, 32, 119, 81, 65, 18, 144, 136, 29, 99, 13, 19, 122, 67, 107, 23, 71, 95, 11, 120, 140, 77, 126, 21, 129, 7, 143, 133, 66, 14, 69, 97, 9, 113])

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

_orig(x::Vector{Float64}) = 100 * (x[1]^2 - x[2])^2 + (1 - x[1])^2 + 90 * (x[3]^2 - x[4])^2 + (1 - x[3])^2 + 10.1 * ((1 - x[2])^2 + (1 - x[4])^2) + 19.8 * (1 - x[2]) * (1 - x[4])

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
        cardinality = Float64[Inf, 151.0, 151.0, Inf],
        f_star = 0.0,
        f_star_justification = "E",
        n_minimizers_documented = 1,
        description = "DFO Benchmark Suite v2.1.0 core (Julia): CDK adaptation of wood (F = f o T).",
        x0_idx = Float64[4.365786015987396, 57.0, 106.0, 631078429.2467339],
        x0_alt_idx = Vector{Float64}[Float64[-3.5688508629798887, 67.0, 32.0, 706959348.6401353], Float64[0.3984675765037533, 56.0, 136.0, 540021325.974652], Float64[-7.536169302463531, 59.0, 61.0, 615902245.3680534], Float64[2.3821267962455757, 63.0, 21.0, 691783164.761455], Float64[-5.55251008272171, 53.0, 91.0, 490156721.8018453], Float64[-1.5851916432380673, 62.0, 52.0, 566037641.1952467], Float64[-9.519828522205353, 64.0, 122.0, 641918560.5886481], Float64[5.640995514392854, 54.0, 82.0, 717799479.9820495], Float64[-2.293641364574432, 58.0, 8.0, 550861457.3165665]],
        suite_release = "2.1.0",
        base_function = "wood",
        m = 151,
        pattern = "CDK",
        strategy = "extreme",
        scale_factors_C = Float64[1.0, 100000000.0],
        lb_orig = Float64[-10.0, -10.0, -10.0, -10.0],
        ub_orig = Float64[6.1211669921875, 6.566458365086623, 6.4227327999999995, 7.3514369012911285],
        lb_C_work = Float64[-10.0, -10.0, -10.0, -1000000000.0],
        ub_C_work = Float64[6.121166992187501, 6.566458365086623, 6.4227327999999995, 735143690.129113],
        ProblemData = _problem_data(),
        func_idx = _objective,
        func_F = xr -> _objective(_cell2idx(xr)),
        idx2cell = _idx2cell,
        cell2idx = _cell2idx,
        decode = xr -> _decode(_cell2idx(xr)))
end

end # module
