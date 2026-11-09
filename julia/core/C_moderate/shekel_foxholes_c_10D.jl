# shekel_foxholes_c_10D  DFO Benchmark Suite v2.1 core, Julia version (self-contained; Base only).
#SHEKEL_FOXHOLES_C_10D  DFO Benchmark Suite v2.1.0 core: shekel_foxholes, dimension 10, pattern C, strategy moderate.
#
# Copyright (C) 2026 J. F. A. Madeira
#
# Objective F(z) = f(T(z)), where f is the original function (embedded below, translated from the MATLAB original) and T acts on the finite
# coordinates only: uniform grid of 0 levels on the interval of the domain, piecewise-affine map with |T - Id| <=
# half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
# Global optimum of F on this domain: f* = -10.20879279215383 (L: literature value at a refined point); the optimal solution is listed in
# reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
# Variables: x1:C x2:C x3:C x4:C x5:C x6:C x7:C x8:C x9:C x10:C
#
# Julia version of the MATLAB wrapper shekel_foxholes_c_10D.m: same instance, same values (level indices 1..m as in MATLAB).
# Usage:  include("shekel_foxholes_c_10D.jl");  PB = shekel_foxholes_c_10D.problem();  f = PB.func_idx(PB.x0_idx)
# Copyright (C) 2026 J. F. A. Madeira. LGPL-3.0-or-later.
module shekel_foxholes_c_10D

const NAME = "shekel_foxholes_c_10D"
const N = 10
const TYPES = Char['C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C']
const CARD = Float64[Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf]
const LB = Float64[0.0, 1.6130671137529848, 1.4526336, 1.5936693044564765, 1.3191721876920974, 1.6072196351668362, 1.8508327247039666, 0.0, 0.0, 0.0]
const UB = Float64[8.12620849609375, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 8.606919836403957, 8.29433714144818, 8.594062897207758]
const LBW = Float64[0.0, 161.30671137529845, 0.014526336000000008, 1.5936693044564763, 131.91721876920968, 0.01607219635166835, 1.8508327247039666, 0.0, 0.0, 0.0]
const UBW = Float64[8.12620849609375, 1000.0, 0.1, 10.0, 1000.0, 0.1, 10.0, 860.6919836403957, 0.0829433714144818, 8.594062897207758]
const TB = Float64[NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN]
const TI = Float64[NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN]
const CW = Float64[8.024966814917121, 915.1928383187155, 0.05113991012282247, 7.620959243851917, 456.4022051654047, 0.04711005066398696, 2.996030602029421, 612.5993024510582, 0.007340577389717023, 4.981999379605735]
const CX = Float64[8.024966814917121, 9.151928383187155, 5.113991012282246, 7.620959243851917, 4.564022051654048, 4.711005066398696, 2.996030602029421, 6.125993024510581, 0.7340577389717022, 4.981999379605735]
const M = 0
const PERM = Dict{Int,Vector{Int}}()

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

const _DATA = [0.806 9.681 0.667 4.783 9.095 3.517 9.325 6.544 0.211 5.122 2.020; 0.517 9.400 2.041 3.788 7.931 2.882 2.672 3.568 1.284 7.033 7.374;
0.100 8.025 9.152 5.114 7.621 4.564 4.711 2.996 6.126 0.734 4.982; 0.908 2.196 0.415 5.649 6.979 9.510 9.166 6.304 6.054 9.377 1.426;
0.965 8.074 8.777 3.467 1.863 6.708 6.349 4.534 0.276 7.633 1.567; 0.669 7.650 5.658 0.720 2.764 3.278 5.283 7.474 6.274 1.409 8.208;
0.524 1.256 3.605 8.623 6.905 4.584 8.133 6.071 6.888 4.187 5.448; 0.902 8.314 2.261 4.224 1.781 4.124 0.932 8.129 8.658 1.208 5.762;
0.531 0.226 8.858 1.420 0.945 1.622 4.698 6.228 9.096 0.972 7.637; 0.876 7.305 2.228 1.242 5.928 9.133 1.826 4.060 5.204 8.713 8.247;
0.462 0.652 7.027 0.508 4.876 8.807 4.632 5.808 6.937 3.291 7.016; 0.491 2.699 3.516 5.874 4.119 4.461 7.496 8.817 0.690 6.593 9.789;
0.463 8.327 3.897 2.017 9.570 9.825 1.150 1.395 3.885 6.354 0.109; 0.714 2.132 7.006 7.136 2.641 1.882 5.943 7.273 7.691 2.880 0.564;
0.352 4.707 5.579 4.080 0.581 9.698 8.542 8.077 8.515 9.231 4.670; 0.869 8.304 7.559 8.567 0.322 7.128 8.392 1.472 8.524 2.277 7.826;
0.813 8.632 4.409 4.832 5.768 7.050 6.715 1.711 4.323 4.405 4.591; 0.811 4.887 9.112 0.170 8.967 9.693 9.867 7.508 7.770 8.382 6.740;
0.828 2.440 6.686 4.299 1.007 7.008 1.427 9.398 8.480 9.950 1.675; 0.964 6.306 8.583 6.084 1.138 4.350 3.134 7.853 6.061 7.457 2.258;
0.789 0.652 2.343 1.370 0.821 1.310 1.063 0.689 8.819 8.833 9.070; 0.360 5.558 1.272 5.756 9.857 2.279 2.764 1.284 1.677 1.244 1.234;
0.369 3.352 7.549 9.817 9.437 8.687 4.167 2.570 6.540 0.228 0.027; 0.992 8.798 0.880 2.370 0.168 1.701 3.680 1.231 2.390 2.499 0.064;
0.332 1.460 8.057 1.336 7.217 7.914 3.615 9.981 9.198 5.292 1.224; 0.817 0.432 8.645 8.774 0.249 8.081 7.461 4.416 0.652 4.002 4.644;
0.632 0.679 2.800 5.523 3.049 2.968 7.225 6.730 4.199 9.614 9.229; 0.883 4.263 1.074 7.286 5.599 8.291 5.200 9.214 8.272 4.398 4.506;
0.608 9.496 4.830 3.150 8.270 5.079 1.231 5.731 9.494 1.883 9.732; 0.326 4.138 2.562 2.532 9.661 5.611 5.500 6.886 2.341 9.699 6.500]
function _orig(x::Vector{Float64})
    d = length(x)
    s = 0.0
    for r in 1:30
        e = 0.0
        for j in 1:d; e += (x[j] - _DATA[r, 1 + j])^2; end
        s += 1 / (_DATA[r, 1] + e)
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
        cardinality = Float64[Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf],
        f_star = -10.208792792153835,
        f_star_justification = "L",
        n_minimizers_documented = 1,
        description = "DFO Benchmark Suite v2.1.0 core (Julia): C adaptation of shekel_foxholes (F = f o T).",
        x0_idx = Float64[6.8604562938213345, 711.6152935318814, 0.033426272583680014, 3.3197505343202653, 660.8542101878204, 0.0883104652178473, 8.84886584794312, 312.5796371553033, 0.06115599448995777, 2.146313383365142],
        x0_alt_idx = Vector{Float64}[Float64[2.860838049650192, 652.9412774552561, 0.09291594272768003, 3.687374325670565, 388.8853748961303, 0.07788156653965225, 6.564843186008939, 546.355576143794, 0.036389555076066726, 3.64848608953582], Float64[4.8606471717357635, 677.101166427984, 0.05065776324608001, 2.8786019846999054, 195.1809526380201, 0.0674526678614572, 3.8031164379185873, 825.3056096228142, 0.01522253212530105, 5.150658795706499], Float64[0.8610289275646209, 692.0572881730063, 0.028092715950080005, 3.246225776050205, 784.1206607157087, 0.050605985381295945, 9.640063240639275, 245.9478477817722, 0.0768420878264189, 6.652831501877177], Float64[7.384215825796127, 633.3832720963811, 0.06723281616896001, 3.6138495674005044, 590.4162384575986, 0.040177086703100884, 7.356040578705094, 524.8978812607926, 0.05567506487565323, 8.155004208047856], Float64[3.3845975816249845, 685.1544627522269, 0.044667768872960006, 3.9814733587508044, 311.3324565264921, 0.023330404222939638, 5.0720179167709105, 803.8479147398126, 0.03450804192488755, 1.073332878957515], Float64[5.384406703710556, 709.3143517249548, 0.08446430683136003, 3.0046444274485795, 978.5365776377607, 0.09633269497030506, 2.7879952548367286, 224.49015289877067, 0.013341018974121878, 2.5755055851281927], Float64[1.3847884595394133, 638.3686460113886, 0.06189925953536002, 3.3722682187988795, 706.5677423460706, 0.08590379629210998, 8.147237971401248, 458.2660918872612, 0.07496057467523973, 4.077678291298872], Float64[6.360504013299941, 662.5285349841165, 0.01964108005376001, 3.7398920101491795, 512.8633200879602, 0.06905711381194873, 5.863215309467066, 737.2161253662814, 0.05379355172447406, 5.57985099746955], Float64[2.3608857691287994, 714.2997256399624, 0.07913075019776003, 2.9311196691785195, 240.89448479627018, 0.058628215133753686, 3.5791926475328832, 157.85836352523953, 0.02902711231058301, 7.082023703640228]],
        suite_release = "2.1.0",
        base_function = "shekel_foxholes",
        m = 0,
        pattern = "C",
        strategy = "moderate",
        scale_factors_C = Float64[1.0, 100.0, 0.01, 1.0, 100.0, 0.01, 1.0, 100.0, 0.01, 1.0],
        lb_orig = Float64[0.0, 1.6130671137529848, 1.4526336, 1.5936693044564765, 1.3191721876920974, 1.6072196351668362, 1.8508327247039666, 0.0, 0.0, 0.0],
        ub_orig = Float64[8.12620849609375, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 8.606919836403957, 8.29433714144818, 8.594062897207758],
        lb_C_work = Float64[0.0, 161.30671137529845, 0.014526336000000008, 1.5936693044564763, 131.91721876920968, 0.01607219635166835, 1.8508327247039666, 0.0, 0.0, 0.0],
        ub_C_work = Float64[8.12620849609375, 1000.0, 0.1, 10.0, 1000.0, 0.1, 10.0, 860.6919836403957, 0.0829433714144818, 8.594062897207758],
        ProblemData = _problem_data(),
        func_idx = _objective,
        func_F = xr -> _objective(_cell2idx(xr)),
        idx2cell = _idx2cell,
        cell2idx = _cell2idx,
        decode = xr -> _decode(_cell2idx(xr)))
end

end # module
