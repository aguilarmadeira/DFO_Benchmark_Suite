# hartman_4_dk_6D  DFO Benchmark Suite v2.1 core, Julia version (self-contained; Base only).
#HARTMAN_4_DK_6D  DFO Benchmark Suite v2.1.0 core: hartman_4, dimension 6, pattern DK, strategy baseline.
#
# Copyright (C) 2026 J. F. A. Madeira
#
# Objective F(z) = f(T(z)), where f is the original function (embedded below, translated from the MATLAB original) and T acts on the finite
# coordinates only: uniform grid of 17 levels on the interval of the domain, piecewise-affine map with |T - Id| <=
# half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
# Global optimum of F on this domain: f* = -3.322368011391339 (L: literature value at a refined point); the optimal solution is listed in
# reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
# Variables: x1:D x2:K x3:D x4:K x5:D x6:K
#
# Julia version of the MATLAB wrapper hartman_4_dk_6D.m: same instance, same values (level indices 1..m as in MATLAB).
# Usage:  include("hartman_4_dk_6D.jl");  PB = hartman_4_dk_6D.problem();  f = PB.func_idx(PB.x0_idx)
# Copyright (C) 2026 J. F. A. Madeira. LGPL-3.0-or-later.
module hartman_4_dk_6D

const NAME = "hartman_4_dk_6D"
const N = 6
const TYPES = Char['D', 'K', 'D', 'K', 'D', 'K']
const CARD = Float64[17.0, 17.0, 17.0, 17.0, 17.0, 17.0]
const LB = Float64[0.0, 0.14795966062084032, 0.14470336, 0.0, 0.18109077248821803, 0.16527817653443508]
const UB = Float64[0.848441162109375, 1.0, 1.0, 0.816347355268638, 1.0, 1.0]
const LBW = Float64[0.0, 0.14795966062084032, 0.14470336, 0.0, 0.18109077248821803, 0.16527817653443508]
const UBW = Float64[0.848441162109375, 1.0, 1.0, 0.816347355268638, 1.0, 1.0]
const TB = Float64[0.21211029052734376, 0.14795966062084032, 0.4654396, 0.2551085485214494, 0.3346362526466772, 0.6348092022338153]
const TI = Float64[0.20169, 0.150011, 0.476874, 0.275332, 0.311652, 0.6573]
const CW = Float64[NaN, NaN, NaN, NaN, NaN, NaN]
const CX = Float64[NaN, NaN, NaN, NaN, NaN, NaN]
const M = 17
const PERM = Dict{Int,Vector{Int}}(2 => Int[7, 13, 4, 8, 11, 6, 17, 2, 10, 1, 9, 5, 15, 16, 14, 12, 3], 4 => Int[7, 2, 5, 15, 10, 16, 11, 17, 9, 3, 14, 1, 12, 13, 6, 4, 8], 6 => Int[17, 15, 13, 14, 7, 1, 8, 16, 9, 4, 5, 10, 12, 2, 11, 6, 3])

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

const _A3 = [3 10 30 1 0.3689 0.1170 0.2673;
             0.1 10 35 1.2 0.4699 0.4387 0.7470;
             3 10 30 3 0.1091 0.8732 0.5547;
             0.1 10 35 3.2 0.03815 0.5743 0.8828]
const _A6 = [10 3 17 3.5 1.7 8 1 0.1312 0.1696 0.5569 0.0124 0.8283 0.5886;
             0.05 10 17 0.1 8 14 1.2 0.2329 0.4135 0.8307 0.3736 0.1004 0.9991;
             3 3.5 1.7 10 17 8 3 0.2348 0.1451 0.3522 0.2883 0.3047 0.6650;
             17 8 0.05 10 0.1 14 3.2 0.4047 0.8828 0.8732 0.5743 0.1091 0.0381]
function _orig(x::Vector{Float64})
    n = length(x)
    A = n == 3 ? _A3 : _A6
    f = 0.0
    for r in 1:4
        e = 0.0
        for j in 1:n; e += A[r, j] * (x[j] - A[r, n + 1 + j])^2; end
        f += A[r, n + 1] * exp(-e)
    end
    return -f
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
        cardinality = Float64[17.0, 17.0, 17.0, 17.0, 17.0, 17.0],
        f_star = -3.322368011391339,
        f_star_justification = "L",
        n_minimizers_documented = 1,
        description = "DFO Benchmark Suite v2.1.0 core (Julia): DK adaptation of hartman_4 (F = f o T).",
        x0_idx = Float64[5.0, 1.0, 4.0, 12.0, 6.0, 16.0],
        x0_alt_idx = Vector{Float64}[Float64[9.0, 2.0, 16.0, 11.0, 1.0, 14.0], Float64[1.0, 2.0, 8.0, 11.0, 13.0, 12.0], Float64[17.0, 1.0, 3.0, 12.0, 8.0, 9.0], Float64[9.0, 2.0, 11.0, 10.0, 5.0, 7.0], Float64[13.0, 2.0, 7.0, 11.0, 16.0, 5.0], Float64[5.0, 1.0, 14.0, 11.0, 12.0, 1.0], Float64[15.0, 2.0, 10.0, 12.0, 7.0, 15.0], Float64[7.0, 3.0, 2.0, 11.0, 3.0, 12.0], Float64[11.0, 1.0, 13.0, 11.0, 14.0, 10.0]],
        suite_release = "2.1.0",
        base_function = "hartman_4",
        m = 17,
        pattern = "DK",
        strategy = "baseline",
        scale_factors_C = Float64[],
        lb_orig = Float64[0.0, 0.14795966062084032, 0.14470336, 0.0, 0.18109077248821803, 0.16527817653443508],
        ub_orig = Float64[0.848441162109375, 1.0, 1.0, 0.816347355268638, 1.0, 1.0],
        lb_C_work = Float64[0.0, 0.14795966062084032, 0.14470336, 0.0, 0.18109077248821803, 0.16527817653443508],
        ub_C_work = Float64[0.848441162109375, 1.0, 1.0, 0.816347355268638, 1.0, 1.0],
        ProblemData = _problem_data(),
        func_idx = _objective,
        func_F = xr -> _objective(_cell2idx(xr)),
        idx2cell = _idx2cell,
        cell2idx = _cell2idx,
        decode = xr -> _decode(_cell2idx(xr)))
end

end # module
