# shekel_47_cdk_4D  DFO Benchmark Suite v2.1 core, Julia version (self-contained; Base only).
#SHEKEL_47_CDK_4D  DFO Benchmark Suite v2.1.0 core: shekel_47, dimension 4, pattern CDK, strategy baseline.
#
# Copyright (C) 2026 J. F. A. Madeira
#
# Objective F(z) = f(T(z)), where f is the original function (embedded below, translated from the MATLAB original) and T acts on the finite
# coordinates only: uniform grid of 151 levels on the interval of the domain, piecewise-affine map with |T - Id| <=
# half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
# Global optimum of F on this domain: f* = -10.40294056681866 (L: literature value at a refined point); the optimal solution is listed in
# reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
# Variables: x1:C x2:K x3:D x4:C
#
# Julia version of the MATLAB wrapper shekel_47_cdk_4D.m: same instance, same values (level indices 1..m as in MATLAB).
# Usage:  include("shekel_47_cdk_4D.jl");  PB = shekel_47_cdk_4D.problem();  f = PB.func_idx(PB.x0_idx)
# Copyright (C) 2026 J. F. A. Madeira. LGPL-3.0-or-later.
module shekel_47_cdk_4D

const NAME = "shekel_47_cdk_4D"
const N = 4
const TYPES = Char['C', 'K', 'D', 'C']
const CARD = Float64[Inf, 151.0, 151.0, Inf]
const LB = Float64[1.3487915039062501, 1.5352893359752071, 1.8390336, 1.9793835901707624]
const UB = Float64[10.0, 10.0, 10.0, 10.0]
const LBW = Float64[1.3487915039062501, 1.5352893359752071, 1.8390336, 1.9793835901707624]
const UBW = Float64[10.0, 10.0, 10.0, 10.0]
const TB = Float64[NaN, 4.018271130755813, 4.015291306666667, NaN]
const TI = Float64[NaN, 4.000689366380763, 3.999489711456878, NaN]
const CW = Float64[4.000572916244606, NaN, NaN, 3.999606160106975]
const CX = Float64[4.000572916244606, NaN, NaN, 3.999606160106975]
const M = 151
const PERM = Dict{Int,Vector{Int}}(2 => Int[137, 121, 31, 29, 100, 86, 62, 6, 143, 18, 148, 117, 35, 40, 124, 49, 65, 44, 38, 126, 134, 87, 111, 129, 17, 76, 72, 58, 9, 97, 20, 59, 78, 33, 2, 94, 82, 25, 53, 54, 39, 5, 146, 15, 22, 51, 136, 123, 73, 91, 36, 89, 102, 151, 1, 107, 23, 125, 110, 13, 45, 115, 63, 116, 37, 56, 4, 24, 32, 70, 135, 145, 139, 141, 103, 85, 104, 83, 30, 128, 34, 147, 109, 61, 21, 46, 81, 138, 47, 3, 19, 41, 88, 130, 149, 119, 64, 101, 113, 50, 55, 27, 66, 42, 92, 68, 14, 16, 74, 52, 8, 120, 108, 132, 118, 12, 112, 75, 131, 142, 150, 140, 98, 57, 122, 105, 99, 96, 79, 67, 77, 43, 144, 90, 7, 84, 127, 11, 10, 80, 28, 106, 114, 60, 48, 71, 133, 95, 93, 26, 69])

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
        cardinality = Float64[Inf, 151.0, 151.0, Inf],
        f_star = -10.402940566818662,
        f_star_justification = "L",
        n_minimizers_documented = 1,
        description = "DFO Benchmark Suite v2.1.0 core (Julia): CDK adaptation of shekel_47 (F = f o T).",
        x0_idx = Float64[2.1640665233135223, 83.0, 117.0, 8.045789004685505],
        x0_alt_idx = Vector{Float64}[Float64[9.260760992765427, 72.0, 49.0, 8.396544824357337], Float64[4.901362961530685, 77.0, 147.0, 8.747300644029169], Float64[7.1317526519298555, 79.0, 79.0, 7.975637840751137], Float64[2.721663945913315, 69.0, 3.0, 8.32639366042297], Float64[8.19625682234764, 78.0, 108.0, 8.677149480094803], Float64[3.8368587911129, 82.0, 34.0, 7.745141159252504], Float64[6.0672484815120695, 70.0, 138.0, 8.095896978924337], Float64[1.6318144381046296, 74.0, 64.0, 8.446652798596169], Float64[9.539559704065322, 83.0, 24.0, 8.797408618268001]],
        suite_release = "2.1.0",
        base_function = "shekel_47",
        m = 151,
        pattern = "CDK",
        strategy = "baseline",
        scale_factors_C = Float64[1.0, 1.0],
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
