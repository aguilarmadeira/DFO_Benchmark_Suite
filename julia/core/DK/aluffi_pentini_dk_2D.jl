# aluffi_pentini_dk_2D  DFO Benchmark Suite v2.1 core, Julia version (self-contained; Base only).
#ALUFFI_PENTINI_DK_2D  DFO Benchmark Suite v2.1.0 core: aluffi_pentini, dimension 2, pattern DK, strategy baseline.
#
# Copyright (C) 2026 J. F. A. Madeira
#
# Objective F(z) = f(T(z)), where f is the original function (embedded below, translated from the MATLAB original) and T acts on the finite
# coordinates only: uniform grid of 151 levels on the interval of the domain, piecewise-affine map with |T - Id| <=
# half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
# Global optimum of F on this domain: f* = -0.3523860738000364 (A: analytic reduction, 1-D root to machine precision); the optimal solution is listed in
# reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
# Variables: x1:D x2:K
#
# Julia version of the MATLAB wrapper aluffi_pentini_dk_2D.m: same instance, same values (level indices 1..m as in MATLAB).
# Usage:  include("aluffi_pentini_dk_2D.jl");  PB = aluffi_pentini_dk_2D.problem();  f = PB.func_idx(PB.x0_idx)
# Copyright (C) 2026 J. F. A. Madeira. LGPL-3.0-or-later.
module aluffi_pentini_dk_2D

const NAME = "aluffi_pentini_dk_2D"
const N = 2
const TYPES = Char['D', 'K']
const CARD = Float64[151.0, 151.0]
const LB = Float64[-6.9031982421875, -10.0]
const UB = Float64[10.0, 7.317349997459736]
const LBW = Float64[-6.9031982421875, -10.0]
const UBW = Float64[10.0, 7.317349997459736]
const TB = Float64[-1.0434228515625001, 0.044062998526646524]
const TI = Float64[-1.04668053249732, 0.0]
const CW = Float64[NaN, NaN]
const CX = Float64[NaN, NaN]
const M = 151
const PERM = Dict{Int,Vector{Int}}(2 => Int[50, 142, 72, 131, 125, 126, 111, 22, 87, 35, 71, 6, 134, 85, 135, 49, 145, 95, 96, 67, 46, 53, 146, 47, 105, 12, 15, 141, 138, 123, 51, 147, 33, 100, 88, 130, 10, 132, 18, 56, 25, 97, 75, 36, 86, 2, 7, 119, 13, 19, 70, 109, 66, 90, 14, 21, 115, 37, 42, 54, 74, 128, 122, 113, 84, 92, 60, 151, 8, 28, 32, 55, 133, 11, 139, 48, 137, 98, 112, 43, 143, 93, 108, 104, 65, 129, 58, 68, 89, 73, 103, 107, 78, 24, 38, 102, 29, 77, 116, 83, 40, 144, 101, 41, 76, 23, 149, 3, 34, 91, 16, 17, 27, 117, 26, 64, 1, 63, 121, 127, 9, 20, 106, 94, 81, 82, 99, 79, 110, 5, 61, 62, 114, 69, 57, 59, 136, 45, 148, 140, 39, 124, 118, 80, 150, 44, 30, 31, 52, 120, 4])

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

_orig(x::Vector{Float64}) = 0.25 * x[1]^4 - 0.5 * x[1]^2 + 0.1 * x[1] + 0.5 * x[2]^2

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
        f_star = -0.35238607380003645,
        f_star_justification = "A",
        n_minimizers_documented = 1,
        description = "DFO Benchmark Suite v2.1.0 core (Julia): DK adaptation of aluffi_pentini (F = f o T).",
        x0_idx = Float64[52.0, 20.0],
        x0_alt_idx = Vector{Float64}[Float64[88.0, 24.0], Float64[14.0, 33.0], Float64[135.0, 23.0], Float64[61.0, 26.0], Float64[97.0, 30.0], Float64[24.0, 19.0], Float64[117.0, 29.0], Float64[43.0, 31.0], Float64[79.0, 20.0]],
        suite_release = "2.1.0",
        base_function = "aluffi_pentini",
        m = 151,
        pattern = "DK",
        strategy = "baseline",
        scale_factors_C = Float64[],
        lb_orig = Float64[-6.9031982421875, -10.0],
        ub_orig = Float64[10.0, 7.317349997459736],
        lb_C_work = Float64[-6.9031982421875, -10.0],
        ub_C_work = Float64[10.0, 7.317349997459736],
        ProblemData = _problem_data(),
        func_idx = _objective,
        func_F = xr -> _objective(_cell2idx(xr)),
        idx2cell = _idx2cell,
        cell2idx = _cell2idx,
        decode = xr -> _decode(_cell2idx(xr)))
end

end # module
