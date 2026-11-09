# griewank_c_10D  DFO Benchmark Suite v2.1 core, Julia version (self-contained; Base only).
#GRIEWANK_C_10D  DFO Benchmark Suite v2.1.0 core: griewank, dimension 10, pattern C, strategy extreme.
#
# Copyright (C) 2026 J. F. A. Madeira
#
# Objective F(z) = f(T(z)), where f is the original function (embedded below, translated from the MATLAB original) and T acts on the finite
# coordinates only: uniform grid of 0 levels on the interval of the domain, piecewise-affine map with |T - Id| <=
# half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
# Global optimum of F on this domain: f* = 0 (E: exact by structure); the optimal solution is listed in
# reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
# Variables: x1:C x2:C x3:C x4:C x5:C x6:C x7:C x8:C x9:C x10:C
#
# Julia version of the MATLAB wrapper griewank_c_10D.m: same instance, same values (level indices 1..m as in MATLAB).
# Usage:  include("griewank_c_10D.jl");  PB = griewank_c_10D.problem();  f = PB.func_idx(PB.x0_idx)
# Copyright (C) 2026 J. F. A. Madeira. LGPL-3.0-or-later.
module griewank_c_10D

const NAME = "griewank_c_10D"
const N = 10
const TYPES = Char['C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C']
const CARD = Float64[Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf]
const LB = Float64[-270.0029296875, -256.7433826144388, -400.0, -400.0, -400.0, -400.0, -400.0, -400.0, -278.2615270814498, -287.5916191725778]
const UB = Float64[400.0, 400.0, 264.525312, 285.07788421491045, 275.49101837306193, 285.0082280032212, 291.6565654146861, 262.337520430322, 400.0, 400.0]
const LBW = Float64[-270.0029296875, -256.7433826144388, -400.0, -400.0, -400.0, -40000000000.0, -40000000000.0, -40000000000.0, -27826152708.14498, -28759161917.257782]
const UBW = Float64[400.0, 400.0, 264.525312, 285.07788421491045, 275.49101837306193, 28500822800.322113, 29165656541.468613, 26233752043.032204, 40000000000.0, 40000000000.0]
const TB = Float64[NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN]
const TI = Float64[NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN]
const CW = Float64[0.0, 0.0, 5.684341886080802e-14, -5.684341886080802e-14, 0.0, 0.0, 0.0, 7.62939453125e-06, 0.0, 0.0]
const CX = Float64[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
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
    n = length(x)
    s = 0.0
    for v in x; s += v^2; end
    p = 1.0
    for i in 1:n; p *= cos(x[i] / sqrt(i)); end
    return 1 + 1 / 4000 * s - p
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
        f_star = 0.0,
        f_star_justification = "E",
        n_minimizers_documented = 1,
        description = "DFO Benchmark Suite v2.1.0 core (Julia): C adaptation of griewank (F = f o T).",
        x0_idx = Float64[64.34423542022705, 312.0138037923956, -24.46345568255998, -303.55838198057484, 240.4730767744585, 2341428021.318817, -29666681884.502754, 15061212151.825287, 4361585085.627441, -11668750685.045858],
        x0_alt_idx = Vector{Float64}[Float64[222.3576216697693, 266.0687934860357, -327.48699795456, -383.1655163812246, 89.74367598046945, -6170508539.667961, 19874117091.566452, -32999056680.09938, -15890935074.534527, 349794858.16091156], Float64[-107.40944528579712, 284.98732714159564, 109.50484721664003, -353.20584214442096, -121.88649281109099, -14682445100.654728, 488587057.45241547, -11532715990.252384, 34497870245.50513, 12368340401.367683], Float64[308.7251868247986, 325.52704211779553, -193.51869505535998, -323.2461679076173, -333.5166616026513, -28432496468.402596, -22951432918.437107, 9933624699.59462, 17188738647.01822, 24386885944.574463], Float64[-21.041880130767822, 272.3749713712224, 238.36959571967998, -389.15745122858533, 185.66238557664428, -36944433029.389366, 26589366057.63211, -34650313656.241455, -120392951.46868134, 36405431487.78123], Float64[139.91585493087769, 291.2935050267823, -64.65394655232001, -359.19777699178167, -31.504216669240805, 22639122897.518044, 7203836023.518066, -13183972966.394459, -17429524549.95558, -20253426073.05069], Float64[-189.85121202468872, 310.21203868234227, -393.19526080512, -342.9339538346597, -182.23361746322985, 8889071529.770172, -12181694010.595978, 8282367723.452545, 32959280770.08407, -8234880529.843922], Float64[351.90896940231323, 254.65761445569802, 69.31435634688006, -312.974279597856, -393.86378625479017, 377134968.7834091, -35621713986.4855, -39777901108.472115, 15650149171.597168, 3783665013.36285], Float64[22.141902446746826, 295.19732943189786, -259.22695790592, -378.88556291882406, 125.31526092450542, -13372916398.964458, 13919084989.58371, -18311560418.62512, -1658982426.8897362, 15802210556.569626], Float64[181.1367383003235, 314.1158630874578, 203.28265924608002, -348.9258886820204, -86.31490786705501, -21884852959.951233, -5466445044.530327, 3154780271.2218857, -21911502587.051712, 27820756099.776398]],
        suite_release = "2.1.0",
        base_function = "griewank",
        m = 0,
        pattern = "C",
        strategy = "extreme",
        scale_factors_C = Float64[1.0, 1.0, 1.0, 1.0, 1.0, 100000000.0, 100000000.0, 100000000.0, 100000000.0, 100000000.0],
        lb_orig = Float64[-270.0029296875, -256.7433826144388, -400.0, -400.0, -400.0, -400.0, -400.0, -400.0, -278.2615270814498, -287.5916191725778],
        ub_orig = Float64[400.0, 400.0, 264.525312, 285.07788421491045, 275.49101837306193, 285.0082280032212, 291.6565654146861, 262.337520430322, 400.0, 400.0],
        lb_C_work = Float64[-270.0029296875, -256.7433826144388, -400.0, -400.0, -400.0, -40000000000.0, -40000000000.0, -40000000000.0, -27826152708.14498, -28759161917.257782],
        ub_C_work = Float64[400.0, 400.0, 264.525312, 285.07788421491045, 275.49101837306193, 28500822800.322113, 29165656541.468613, 26233752043.032204, 40000000000.0, 40000000000.0],
        ProblemData = _problem_data(),
        func_idx = _objective,
        func_F = xr -> _objective(_cell2idx(xr)),
        idx2cell = _idx2cell,
        cell2idx = _cell2idx,
        decode = xr -> _decode(_cell2idx(xr)))
end

end # module
