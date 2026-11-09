# schwefel_c_10D  DFO Benchmark Suite v2.1 core, Julia version (self-contained; Base only).
#SCHWEFEL_C_10D  DFO Benchmark Suite v2.1.0 core: schwefel, dimension 10, pattern C, strategy moderate.
#
# Copyright (C) 2026 J. F. A. Madeira
#
# Objective F(z) = f(T(z)), where f is the original function (embedded below, translated from the MATLAB original) and T acts on the finite
# coordinates only: uniform grid of 0 levels on the interval of the domain, piecewise-affine map with |T - Id| <=
# half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
# Global optimum of F on this domain: f* = -4189.828872721624 (A: analytic reduction, 1-D root to machine precision); the optimal solution is listed in
# reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
# Variables: x1:C x2:C x3:C x4:C x5:C x6:C x7:C x8:C x9:C x10:C
#
# Julia version of the MATLAB wrapper schwefel_c_10D.m: same instance, same values (level indices 1..m as in MATLAB).
# Usage:  include("schwefel_c_10D.jl");  PB = schwefel_c_10D.problem();  f = PB.func_idx(PB.x0_idx)
# Copyright (C) 2026 J. F. A. Madeira. LGPL-3.0-or-later.
module schwefel_c_10D

const NAME = "schwefel_c_10D"
const N = 10
const TYPES = Char['C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C']
const CARD = Float64[Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf]
const LB = Float64[-343.245849609375, -300.1884875273078, -344.09664, -322.0616409829238, -330.47947544566625, -366.20111340639335, -335.504962823721, -309.3069420891491, -344.6511054491658, -301.55848948296364]
const UB = Float64[500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0]
const LBW = Float64[-343.245849609375, -30018.848752730784, -3.4409664, -322.0616409829238, -33047.94754456663, -3.662011134063933, -335.5049628237209, -30930.69420891491, -3.4465110544916584, -301.55848948296364]
const UBW = Float64[500.0, 50000.0, 5.0, 500.0, 50000.0, 5.0, 500.0, 50000.0, 5.0, 500.0]
const TB = Float64[NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN]
const TI = Float64[NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN]
const CW = Float64[420.9687, 42096.869999999995, 4.209687000000001, 420.9687, 42096.87000000001, 4.209687000000001, 420.9687000000001, 42096.87000000001, 4.209687, 420.9687]
const CX = Float64[420.9687, 420.9687, 420.9687, 420.9687, 420.9687, 420.9687, 420.9687, 420.9687, 420.9687, 420.9687]
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

function _orig(x::Vector{Float64})
    s = 0.0
    for v in x; s += v * sin(sqrt(abs(v))); end
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
        f_star = -4189.828872721624,
        f_star_justification = "A",
        n_minimizers_documented = 1,
        description = "DFO Benchmark Suite v2.1.0 core (Julia): C adaptation of schwefel (F = f o T).",
        x0_idx = Float64[-0.2654820680618286, 31376.50022170097, -0.2914730168320001, 64.83117630599907, -20194.546196572097, 0.46200509089322406, 136.24158040302484, -22836.444868294464, 0.9450915534408004, -239.57401651708273],
        x0_alt_idx = Vector{Float64}[Float64[217.13383853435516, 25778.473765748615, 3.962774048768001, 100.78139384169555, 36148.27621720977, -0.6143394878957853, -97.93109083822569, -854.5403178923007, -1.5769906904129178, -99.46807364669556], Float64[-207.78301537036896, 28083.543482905476, 1.4102258094080002, 136.73161137739186, 17616.91602131474, -1.6906840666847949, 451.5328893945125, 25375.075256742173, 4.698016576912643, 40.63786922369172], Float64[325.83349883556366, 29510.491403050193, -2.438854868992, 57.641132798859815, -8401.862435547988, -3.4293945401131944, 217.36021815326194, -29101.818629037112, 2.5424797293580603, 180.74381209407892], Float64[-104.0242487192154, 23912.46494709783, 3.0470980136960013, 93.5913503345563, -26933.222631443015, 4.105017511409873, -16.812453087988615, -2872.203054402642, 0.386942881803479, 320.8497549644661], Float64[113.3750718832016, 28851.90005529109, -0.8019826647039996, 129.5415678702526, 28728.92539497815, 2.3663070379814712, -250.98512432923908, 23357.412520231832, -1.7685939657511027, 460.95569783485337], Float64[-311.5417820215225, 31156.969772447945, 4.748796839936, 34.016704132545044, 10197.565199083117, 1.2899624591924628, 347.4561466206234, 49587.0280948663, 4.506413301574458, -199.54374712554352], Float64[459.2376273870468, 24388.114253812735, 0.8997161615360008, 69.96692166824141, -15821.213257779604, 0.21361788040345298, 64.30618466224848, -9137.57681514529, 2.3508764540198754, -59.43780425515632], Float64[26.909433007240295, 26693.183970969596, -3.273497626624, 105.91713920393778, 48009.02741696997, -1.5250925930249468, -169.86648657900201, 17092.038759489187, 0.1953396064652937, 80.66813861523087], Float64[244.3087536096573, 31632.61907916285, 2.601414987776, 141.8673567396341, 21990.24896010723, -2.6014371718139566, 428.57478437086047, 43321.654334123654, -2.326742637388425, 220.77408148561813]],
        suite_release = "2.1.0",
        base_function = "schwefel",
        m = 0,
        pattern = "C",
        strategy = "moderate",
        scale_factors_C = Float64[1.0, 100.0, 0.01, 1.0, 100.0, 0.01, 1.0, 100.0, 0.01, 1.0],
        lb_orig = Float64[-343.245849609375, -300.1884875273078, -344.09664, -322.0616409829238, -330.47947544566625, -366.20111340639335, -335.504962823721, -309.3069420891491, -344.6511054491658, -301.55848948296364],
        ub_orig = Float64[500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0],
        lb_C_work = Float64[-343.245849609375, -30018.848752730784, -3.4409664, -322.0616409829238, -33047.94754456663, -3.662011134063933, -335.5049628237209, -30930.69420891491, -3.4465110544916584, -301.55848948296364],
        ub_C_work = Float64[500.0, 50000.0, 5.0, 500.0, 50000.0, 5.0, 500.0, 50000.0, 5.0, 500.0],
        ProblemData = _problem_data(),
        func_idx = _objective,
        func_F = xr -> _objective(_cell2idx(xr)),
        idx2cell = _idx2cell,
        cell2idx = _cell2idx,
        decode = xr -> _decode(_cell2idx(xr)))
end

end # module
