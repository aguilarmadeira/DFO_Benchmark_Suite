# schaffer1_dk_2D  DFO Benchmark Suite v2.1 core, Julia version (self-contained; Base only).
#SCHAFFER1_DK_2D  DFO Benchmark Suite v2.1.0 core: schaffer1, dimension 2, pattern DK, strategy baseline.
#
# Copyright (C) 2026 J. F. A. Madeira
#
# Objective F(z) = f(T(z)), where f is the original function (embedded below, translated from the MATLAB original) and T acts on the finite
# coordinates only: uniform grid of 151 levels on the interval of the domain, piecewise-affine map with |T - Id| <=
# half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
# Global optimum of F on this domain: f* = 0 (E: exact by structure); the optimal solution is listed in
# reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
# Variables: x1:D x2:K
#
# Julia version of the MATLAB wrapper schaffer1_dk_2D.m: same instance, same values (level indices 1..m as in MATLAB).
# Usage:  include("schaffer1_dk_2D.jl");  PB = schaffer1_dk_2D.problem();  f = PB.func_idx(PB.x0_idx)
# Copyright (C) 2026 J. F. A. Madeira. LGPL-3.0-or-later.
module schaffer1_dk_2D

const NAME = "schaffer1_dk_2D"
const N = 2
const TYPES = Char['D', 'K']
const CARD = Float64[151.0, 151.0]
const LB = Float64[-72.149169921875, -100.0]
const UB = Float64[100.0, 69.3710308387949]
const LBW = Float64[-72.149169921875, -100.0]
const UBW = Float64[100.0, 69.3710308387949]
const TB = Float64[0.15348144531249375, 0.4934782976849732]
const TI = Float64[0.0, 0.0]
const CW = Float64[NaN, NaN]
const CX = Float64[NaN, NaN]
const M = 151
const PERM = Dict{Int,Vector{Int}}(2 => Int[117, 81, 128, 63, 87, 43, 129, 35, 44, 28, 116, 85, 113, 29, 102, 68, 124, 46, 33, 75, 110, 119, 57, 23, 104, 127, 12, 107, 115, 151, 98, 150, 8, 53, 139, 18, 56, 1, 3, 138, 47, 137, 90, 74, 38, 59, 78, 16, 112, 148, 82, 9, 134, 66, 14, 86, 131, 54, 62, 120, 118, 109, 10, 22, 45, 83, 11, 31, 39, 5, 106, 52, 13, 123, 67, 20, 136, 80, 25, 71, 126, 93, 65, 48, 147, 94, 60, 77, 24, 143, 76, 6, 95, 122, 51, 32, 108, 2, 55, 133, 89, 58, 96, 69, 4, 73, 17, 140, 145, 37, 146, 130, 72, 132, 149, 61, 42, 141, 97, 111, 19, 27, 30, 103, 26, 125, 100, 144, 121, 79, 36, 114, 88, 41, 64, 40, 91, 7, 84, 15, 99, 50, 135, 142, 70, 101, 92, 34, 21, 105, 49])

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
    r = sqrt(x[1]^2 + x[2]^2)
    return 0.5 + (sin(r)^2 - 0.5) / ((1 + 0.001 * (x[1]^2 + x[2]^2))^2)
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
        cardinality = Float64[151.0, 151.0],
        f_star = 0.0,
        f_star_justification = "E",
        n_minimizers_documented = 1,
        description = "DFO Benchmark Suite v2.1.0 core (Julia): DK adaptation of schaffer1 (F = f o T).",
        x0_idx = Float64[25.0, 16.0],
        x0_alt_idx = Vector{Float64}[Float64[119.0, 6.0], Float64[44.0, 10.0], Float64[83.0, 13.0], Float64[6.0, 2.0], Float64[143.0, 11.0], Float64[68.0, 16.0], Float64[106.0, 3.0], Float64[29.0, 7.0], Float64[124.0, 17.0]],
        suite_release = "2.1.0",
        base_function = "schaffer1",
        m = 151,
        pattern = "DK",
        strategy = "baseline",
        scale_factors_C = Float64[],
        lb_orig = Float64[-72.149169921875, -100.0],
        ub_orig = Float64[100.0, 69.3710308387949],
        lb_C_work = Float64[-72.149169921875, -100.0],
        ub_C_work = Float64[100.0, 69.3710308387949],
        ProblemData = _problem_data(),
        func_idx = _objective,
        func_F = xr -> _objective(_cell2idx(xr)),
        idx2cell = _idx2cell,
        cell2idx = _cell2idx,
        decode = xr -> _decode(_cell2idx(xr)))
end

end # module
