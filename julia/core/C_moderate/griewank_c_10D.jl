# griewank_c_10D  DFO Benchmark Suite v2.1 core, Julia version (self-contained; Base only).
#GRIEWANK_C_10D  DFO Benchmark Suite v2.1.0 core: griewank, dimension 10, pattern C, strategy moderate.
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
const LBW = Float64[-270.0029296875, -25674.338261443885, -4.0, -400.0, -40000.0, -4.0, -400.0, -40000.0, -2.7826152708144978, -287.5916191725778]
const UBW = Float64[400.0, 40000.0, 2.6452531199999996, 285.07788421491045, 27549.101837306196, 2.850082280032211, 291.65656541468616, 26233.752043032204, 4.0, 400.0]
const TB = Float64[NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN]
const TI = Float64[NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN]
const CW = Float64[0.0, 0.0, 0.0, -5.684341886080802e-14, 7.275957614183426e-12, -4.440892098500626e-16, 5.684341886080802e-14, 7.275957614183426e-12, 0.0, 0.0]
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
        x0_idx = Float64[64.34423542022705, 31201.38037923956, -0.24463455682559987, -303.55838198057484, 24047.307677445853, 0.2341428021318812, -296.66681884502754, 15061.21215182529, 0.43615850856274463, -116.68750685045856],
        x0_alt_idx = Vector{Float64}[Float64[222.3576216697693, 26606.879348603576, -3.2748699795456, -383.1655163812246, 8974.367598046949, -0.6170508539667958, 198.7411709156646, -32999.056680099384, -1.5890935074534525, 3.497948581609137], Float64[-107.40944528579712, 28498.732714159563, 1.0950484721664004, -353.20584214442096, -12188.649281109098, -1.4682445100654729, 4.885870574524176, -11532.715990252382, 3.4497870245505124, 123.6834040136768], Float64[308.7251868247986, 32552.70421177956, -1.9351869505535997, -323.2461679076173, -33351.66616026513, -2.8432496468402593, -229.51432918437106, 9933.62469959462, 1.7188738647018225, 243.86885944574465], Float64[-21.041880130767822, 27237.497137122245, 2.3836959571967995, -389.15745122858533, 18566.23855766443, -3.6944433029389363, 265.89366057632105, -34650.31365624146, -0.012039295146867968, 364.0543148778123], Float64[139.91585493087769, 29129.35050267824, -0.6465394655232002, -359.19777699178167, -3150.421666924085, 2.2639122897518043, 72.03836023518073, -13183.97296639446, -1.7429524549955582, -202.5342607305069], Float64[-189.85121202468872, 31021.203868234235, -3.9319526080512, -342.9339538346597, -18223.361746322982, 0.8889071529770174, -121.81694010595976, 8282.367723452546, 3.2959280770084067, -82.3488052984392], Float64[351.90896940231323, 25465.761445569802, 0.6931435634688006, -312.974279597856, -39386.37862547902, 0.037713496878341246, -356.21713986485497, -39777.90110847212, 1.5650149171597167, 37.836650133628495], Float64[22.141902446746826, 29519.732943189792, -2.5922695790591996, -378.88556291882406, 12531.52609245055, -1.337291639896446, 139.19084989583712, -18311.56041862512, -0.16589824268897324, 158.02210556569628], Float64[181.1367383003235, 31411.586308745787, 2.0328265924607996, -348.9258886820204, -8631.490786705497, -2.188485295995123, -54.66445044530326, 3154.7802712218836, -2.1911502587051706, 278.20756099776395]],
        suite_release = "2.1.0",
        base_function = "griewank",
        m = 0,
        pattern = "C",
        strategy = "moderate",
        scale_factors_C = Float64[1.0, 100.0, 0.01, 1.0, 100.0, 0.01, 1.0, 100.0, 0.01, 1.0],
        lb_orig = Float64[-270.0029296875, -256.7433826144388, -400.0, -400.0, -400.0, -400.0, -400.0, -400.0, -278.2615270814498, -287.5916191725778],
        ub_orig = Float64[400.0, 400.0, 264.525312, 285.07788421491045, 275.49101837306193, 285.0082280032212, 291.6565654146861, 262.337520430322, 400.0, 400.0],
        lb_C_work = Float64[-270.0029296875, -25674.338261443885, -4.0, -400.0, -40000.0, -4.0, -400.0, -40000.0, -2.7826152708144978, -287.5916191725778],
        ub_C_work = Float64[400.0, 40000.0, 2.6452531199999996, 285.07788421491045, 27549.101837306196, 2.850082280032211, 291.65656541468616, 26233.752043032204, 4.0, 400.0],
        ProblemData = _problem_data(),
        func_idx = _objective,
        func_F = xr -> _objective(_cell2idx(xr)),
        idx2cell = _idx2cell,
        cell2idx = _cell2idx,
        decode = xr -> _decode(_cell2idx(xr)))
end

end # module
