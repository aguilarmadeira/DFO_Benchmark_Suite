# shekel_47_k_4D  DFO Benchmark Suite v2.1 core, Julia version (self-contained; Base only).
#SHEKEL_47_K_4D  DFO Benchmark Suite v2.1.0 core: shekel_47, dimension 4, pattern K, strategy baseline.
#
# Copyright (C) 2026 J. F. A. Madeira
#
# Objective F(z) = f(T(z)), where f is the original function (embedded below, translated from the MATLAB original) and T acts on the finite
# coordinates only: uniform grid of 25 levels on the interval of the domain, piecewise-affine map with |T - Id| <=
# half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
# Global optimum of F on this domain: f* = -10.40294056681866 (L: literature value at a refined point); the optimal solution is listed in
# reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
# Variables: x1:K x2:K x3:K x4:K
#
# Julia version of the MATLAB wrapper shekel_47_k_4D.m: same instance, same values (level indices 1..m as in MATLAB).
# Usage:  include("shekel_47_k_4D.jl");  PB = shekel_47_k_4D.problem();  f = PB.func_idx(PB.x0_idx)
# Copyright (C) 2026 J. F. A. Madeira. LGPL-3.0-or-later.
module shekel_47_k_4D

const NAME = "shekel_47_k_4D"
const N = 4
const TYPES = Char['K', 'K', 'K', 'K']
const CARD = Float64[25.0, 25.0, 25.0, 25.0]
const LB = Float64[1.3487915039062501, 1.5352893359752071, 1.8390336, 1.9793835901707624]
const UB = Float64[10.0, 10.0, 10.0, 10.0]
const LBW = Float64[1.3487915039062501, 1.5352893359752071, 1.8390336, 1.9793835901707624]
const UBW = Float64[10.0, 10.0, 10.0, 10.0]
const TB = Float64[3.87206064860026, 4.004163279649105, 3.8792752, 3.984537692628072]
const TI = Float64[4.000572916244606, 4.000689366380763, 3.999489711456878, 3.999606160106975]
const CW = Float64[NaN, NaN, NaN, NaN]
const CX = Float64[NaN, NaN, NaN, NaN]
const M = 25
const PERM = Dict{Int,Vector{Int}}(1 => Int[22, 12, 24, 16, 7, 2, 3, 14, 18, 5, 8, 10, 11, 4, 21, 17, 6, 15, 25, 23, 19, 1, 13, 9, 20], 2 => Int[9, 20, 17, 24, 23, 2, 4, 6, 16, 21, 11, 13, 1, 18, 3, 25, 22, 8, 7, 10, 19, 14, 15, 5, 12], 3 => Int[1, 11, 8, 22, 12, 10, 14, 20, 3, 15, 25, 4, 17, 13, 9, 2, 5, 23, 19, 7, 18, 6, 24, 21, 16], 4 => Int[22, 2, 24, 7, 21, 16, 6, 10, 20, 23, 5, 17, 8, 3, 14, 9, 13, 11, 12, 25, 18, 19, 1, 15, 4])

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

const _A = [4 4 4 4 0.1; 1 1 1 1 0.2; 8 8 8 8 0.2; 6 6 6 6 0.4; 3 7 3 7 0.4;
            2 9 2 9 0.6; 5 5 3 3 0.3; 8 1 8 1 0.7; 6 2 6 2 0.5; 7 3.6 7 3.6 0.5]
function _orig(x::Vector{Float64})
    s = 0.0
    for r in 1:7
        e = 0.0
        for j in 1:4; e += (x[j] - _A[r, j])^2; end
        s += 1 / (e + _A[r, 5])
    end
    return -s
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
        cardinality = Float64[25.0, 25.0, 25.0, 25.0],
        f_star = -10.402940566818662,
        f_star_justification = "L",
        n_minimizers_documented = 1,
        description = "DFO Benchmark Suite v2.1.0 core (Julia): K adaptation of shekel_47 (F = f o T).",
        x0_idx = Float64[3.0, 14.0, 20.0, 19.0],
        x0_alt_idx = Vector{Float64}[Float64[23.0, 12.0, 9.0, 20.0], Float64[11.0, 13.0, 24.0, 21.0], Float64[17.0, 14.0, 13.0, 19.0], Float64[5.0, 12.0, 1.0, 20.0], Float64[20.0, 13.0, 18.0, 21.0], Float64[8.0, 14.0, 6.0, 18.0], Float64[14.0, 12.0, 23.0, 19.0], Float64[2.0, 13.0, 11.0, 20.0], Float64[24.0, 14.0, 5.0, 21.0]],
        suite_release = "2.1.0",
        base_function = "shekel_47",
        m = 25,
        pattern = "K",
        strategy = "baseline",
        scale_factors_C = Float64[],
        lb_orig = Float64[1.3487915039062501, 1.5352893359752071, 1.8390336, 1.9793835901707624],
        ub_orig = Float64[10.0, 10.0, 10.0, 10.0],
        lb_C_work = Float64[1.3487915039062501, 1.5352893359752071, 1.8390336, 1.9793835901707624],
        ub_C_work = Float64[10.0, 10.0, 10.0, 10.0],
        ProblemData = _problem_data(),
        func_idx = _objective,
        func_F = xr -> _objective(_cell2idx(xr)),
        idx2cell = _idx2cell,
        cell2idx = _cell2idx,
        decode = xr -> _decode(_cell2idx(xr)))
end

end # module
