# shekel_45_cdk_4D  DFO Benchmark Suite v2.1 core, Julia version (self-contained; Base only).
#SHEKEL_45_CDK_4D  DFO Benchmark Suite v2.1.0 core: shekel_45, dimension 4, pattern CDK, strategy baseline.
#
# Copyright (C) 2026 J. F. A. Madeira
#
# Objective F(z) = f(T(z)), where f is the original function (embedded below, translated from the MATLAB original) and T acts on the finite
# coordinates only: uniform grid of 151 levels on the interval of the domain, piecewise-affine map with |T - Id| <=
# half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
# Global optimum of F on this domain: f* = -10.15319967905823 (L: literature value at a refined point); the optimal solution is listed in
# reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
# Variables: x1:C x2:K x3:D x4:C
#
# Julia version of the MATLAB wrapper shekel_45_cdk_4D.m: same instance, same values (level indices 1..m as in MATLAB).
# Usage:  include("shekel_45_cdk_4D.jl");  PB = shekel_45_cdk_4D.problem();  f = PB.func_idx(PB.x0_idx)
# Copyright (C) 2026 J. F. A. Madeira. LGPL-3.0-or-later.
module shekel_45_cdk_4D

const NAME = "shekel_45_cdk_4D"
const N = 4
const TYPES = Char['C', 'K', 'D', 'C']
const CARD = Float64[Inf, 151.0, 151.0, Inf]
const LB = Float64[1.9175415039062502, 1.3019560026418735, 1.6990336, 1.8793835901707623]
const UB = Float64[10.0, 10.0, 10.0, 10.0]
const LBW = Float64[1.9175415039062502, 1.3019560026418735, 1.6990336, 1.8793835901707623]
const UBW = Float64[10.0, 10.0, 10.0, 10.0]
const TB = Float64[NaN, 4.027343121814086, 4.0233041919999994, NaN]
const TI = Float64[NaN, 4.000133276236526, 4.000037153472134, NaN]
const CW = Float64[4.000037153719125, NaN, NaN, 4.00013327718805]
const CX = Float64[4.000037153719125, NaN, NaN, 4.00013327718805]
const M = 151
const PERM = Dict{Int,Vector{Int}}(2 => Int[88, 102, 74, 25, 151, 141, 115, 132, 143, 47, 68, 116, 94, 67, 6, 83, 63, 150, 45, 73, 39, 76, 29, 58, 104, 117, 24, 49, 5, 8, 136, 109, 65, 93, 61, 35, 44, 89, 21, 26, 33, 95, 114, 7, 15, 122, 52, 41, 60, 134, 27, 128, 139, 77, 46, 147, 85, 30, 149, 133, 32, 105, 50, 71, 23, 96, 80, 79, 106, 91, 62, 31, 66, 10, 48, 2, 140, 145, 92, 142, 12, 130, 51, 37, 81, 20, 121, 72, 19, 64, 78, 100, 135, 131, 111, 42, 28, 120, 75, 146, 13, 86, 22, 101, 119, 16, 124, 70, 18, 103, 59, 138, 144, 90, 125, 53, 36, 126, 123, 17, 54, 97, 112, 69, 4, 110, 99, 84, 57, 38, 34, 127, 1, 11, 55, 107, 148, 113, 98, 137, 14, 56, 82, 43, 87, 9, 3, 118, 40, 129, 108])

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
    for r in 1:5
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
        f_star = -10.153199679058229,
        f_star_justification = "L",
        n_minimizers_documented = 1,
        description = "DFO Benchmark Suite v2.1.0 core (Julia): CDK adaptation of shekel_45 (F = f o T).",
        x0_idx = Float64[9.246216028928757, 33.0, 87.0, 6.861336098158462],
        x0_alt_idx = Vector{Float64}[Float64[5.268130987882614, 22.0, 19.0, 7.216465095672859], Float64[7.257173508405685, 27.0, 117.0, 7.571594093187258], Float64[3.279088467359543, 29.0, 49.0, 6.790310298655582], Float64[8.251694768667221, 19.0, 146.0, 7.145439296169981], Float64[4.2736097276210785, 28.0, 78.0, 7.500568293684379], Float64[6.26265224814415, 32.0, 4.0, 6.556939814574692], Float64[2.284567207098007, 20.0, 108.0, 6.912068812089089], Float64[9.506685882806778, 24.0, 34.0, 7.267197809603488], Float64[5.528600841760635, 33.0, 138.0, 7.622326807117885]],
        suite_release = "2.1.0",
        base_function = "shekel_45",
        m = 151,
        pattern = "CDK",
        strategy = "baseline",
        scale_factors_C = Float64[1.0, 1.0],
        lb_orig = Float64[1.9175415039062502, 1.3019560026418735, 1.6990336, 1.8793835901707623],
        ub_orig = Float64[10.0, 10.0, 10.0, 10.0],
        lb_C_work = Float64[1.9175415039062502, 1.3019560026418735, 1.6990336, 1.8793835901707623],
        ub_C_work = Float64[10.0, 10.0, 10.0, 10.0],
        ProblemData = _problem_data(),
        func_idx = _objective,
        func_F = xr -> _objective(_cell2idx(xr)),
        idx2cell = _idx2cell,
        cell2idx = _cell2idx,
        decode = xr -> _decode(_cell2idx(xr)))
end

end # module
