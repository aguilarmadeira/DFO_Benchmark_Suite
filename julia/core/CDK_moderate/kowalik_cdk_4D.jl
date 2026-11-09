# kowalik_cdk_4D  DFO Benchmark Suite v2.1 core, Julia version (self-contained; Base only).
#KOWALIK_CDK_4D  DFO Benchmark Suite v2.1.0 core: kowalik, dimension 4, pattern CDK, strategy moderate.
#
# Copyright (C) 2026 J. F. A. Madeira
#
# Objective F(z) = f(T(z)), where f is the original function (embedded below, translated from the MATLAB original) and T acts on the finite
# coordinates only: uniform grid of 151 levels on the interval of the domain, piecewise-affine map with |T - Id| <=
# half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
# Global optimum of F on this domain: f* = 0.0003074859878056062 (L: literature value at a refined point); the optimal solution is listed in
# reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
# Variables: x1:C x2:K x3:D x4:C
#
# Julia version of the MATLAB wrapper kowalik_cdk_4D.m: same instance, same values (level indices 1..m as in MATLAB).
# Usage:  include("kowalik_cdk_4D.jl");  PB = kowalik_cdk_4D.problem();  f = PB.func_idx(PB.x0_idx)
# Copyright (C) 2026 J. F. A. Madeira. LGPL-3.0-or-later.
module kowalik_cdk_4D

const NAME = "kowalik_cdk_4D"
const N = 4
const TYPES = Char['C', 'K', 'D', 'C']
const CARD = Float64[Inf, 151.0, 151.0, Inf]
const LB = Float64[0.059979711914062504, 0.08174305746075294, 0.0, 0.05673411078717201]
const UB = Float64[0.42, 0.42, 0.3474645888, 0.42]
const LBW = Float64[0.059979711914062483, 0.08174305746075294, 0.0, 5.673411078717201]
const UBW = Float64[0.42, 0.42, 0.3474645888, 42.0]
const TB = Float64[NaN, 0.189985279073312, 0.122770821376, NaN]
const TI = Float64[NaN, 0.1908362411132952, 0.1231173000994744, NaN]
const CW = Float64[0.19283345290075488, NaN, NaN, 13.576599043790779]
const CX = Float64[0.1928334529007549, NaN, NaN, 0.1357659904379078]
const M = 151
const PERM = Dict{Int,Vector{Int}}(2 => Int[14, 97, 127, 96, 88, 54, 51, 125, 18, 5, 91, 126, 21, 114, 52, 131, 86, 4, 70, 150, 28, 105, 149, 46, 23, 61, 6, 72, 13, 98, 121, 24, 48, 77, 136, 44, 64, 20, 56, 117, 69, 49, 40, 3, 120, 137, 111, 109, 74, 142, 8, 65, 146, 47, 12, 95, 116, 62, 106, 103, 35, 118, 92, 68, 76, 31, 144, 60, 66, 27, 134, 84, 101, 90, 138, 16, 85, 122, 7, 45, 75, 15, 25, 102, 145, 36, 63, 78, 19, 37, 32, 50, 141, 53, 139, 130, 99, 112, 108, 38, 135, 81, 10, 73, 71, 94, 129, 11, 57, 148, 26, 82, 39, 151, 34, 29, 110, 147, 58, 9, 113, 59, 128, 41, 17, 79, 100, 123, 30, 124, 1, 132, 143, 115, 140, 2, 89, 133, 104, 22, 87, 119, 67, 42, 33, 43, 80, 83, 107, 93, 55])

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

const _a = [0.1957, 0.1947, 0.1735, 0.16, 0.0844, 0.0627, 0.0456, 0.0342, 0.0323, 0.0235, 0.0246]
const _b = [0.25, 0.5, 1, 2, 4, 6, 8, 10, 12, 14, 16]
function _orig(x::Vector{Float64})
    s = 0.0
    for k in 1:11; s += (_a[k] - x[1] * (1 + x[2] * _b[k]) / (1 + x[3] * _b[k] + x[4] * _b[k]^2))^2; end
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
        cardinality = Float64[Inf, 151.0, 151.0, Inf],
        f_star = 0.00030748598780560617,
        f_star_justification = "L",
        n_minimizers_documented = 1,
        description = "DFO Benchmark Suite v2.1.0 core (Julia): CDK adaptation of kowalik (F = f o T).",
        x0_idx = Float64[0.10515803908109664, 105.0, 93.0, 41.924351126777836],
        x0_alt_idx = Vector{Float64}[Float64[0.33298337763547897, 109.0, 53.0, 37.703144000981126], Float64[0.1515669043421745, 114.0, 123.0, 39.29177033864655], Float64[0.24438463486433026, 101.0, 84.0, 40.88039667631198], Float64[0.06085866769552229, 111.0, 8.0, 37.38541873344804], Float64[0.415253638780117, 115.0, 113.0, 38.97404507111347], Float64[0.23805615323781965, 105.0, 39.0, 39.83644222584613], Float64[0.32665489600896835, 107.0, 143.0, 41.42506856351156], Float64[0.149457410466671, 116.0, 69.0, 37.93009062064761], Float64[0.3698995204567909, 105.0, 29.0, 39.518716958313036]],
        suite_release = "2.1.0",
        base_function = "kowalik",
        m = 151,
        pattern = "CDK",
        strategy = "moderate",
        scale_factors_C = Float64[1.0, 100.0],
        lb_orig = Float64[0.059979711914062504, 0.08174305746075294, 0.0, 0.05673411078717201],
        ub_orig = Float64[0.42, 0.42, 0.3474645888, 0.42],
        lb_C_work = Float64[0.059979711914062483, 0.08174305746075294, 0.0, 5.673411078717201],
        ub_C_work = Float64[0.42, 0.42, 0.3474645888, 42.0],
        ProblemData = _problem_data(),
        func_idx = _objective,
        func_F = xr -> _objective(_cell2idx(xr)),
        idx2cell = _idx2cell,
        cell2idx = _cell2idx,
        decode = xr -> _decode(_cell2idx(xr)))
end

end # module
