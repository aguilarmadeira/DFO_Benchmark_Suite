# interBK_nC2_K8x2_s3_b0p36_l0  DFO Benchmark Suite v2.1 core, Julia version (self-contained; Base only).
#INTERBK_NC2_K8X2_S3_B0P36_L0  DFO Benchmark Suite v2.1.0 core, Family BK: f = ||x_C - x_C* - lambda s(kappa)||^2 + b (1 - h(kappa)/nK).
#
# Copyright (C) 2026 J. F. A. Madeira
#
# nC = 2 continuous variables in [0,1]; nK = 8 categorical variables with m = 2 categories each (M_K = 256 combinations)
# lambda = 0, sigma = 0.15, b = 0.36, q = lambda^2 sigma^2 / b = 0, j_min = 0 (minimal number of simultaneous categorical
# changes from an antipodal conditional minimiser that improves f without continuous re-optimisation; Inf: none; 0: lambda = 0)
# h(kappa) = number of categorical coordinates equal to the optimal combination; s(kappa) = s if h = 0, else 0.
# Global optimum f* = 0 (E: by construction); the optimal solution is listed in reference_optima.csv of this
# folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
#
# Julia version of the MATLAB wrapper interBK_nC2_K8x2_s3_b0p36_l0.m: same instance, same values (level indices 1..m as in MATLAB).
# Usage:  include("interBK_nC2_K8x2_s3_b0p36_l0.jl");  PB = interBK_nC2_K8x2_s3_b0p36_l0.problem();  f = PB.func_idx(PB.x0_idx)
# Copyright (C) 2026 J. F. A. Madeira. LGPL-3.0-or-later.
module interBK_nC2_K8x2_s3_b0p36_l0

const NAME = "interBK_nC2_K8x2_s3_b0p36_l0"
const N = 10
const TYPES = Char['C', 'C', 'K', 'K', 'K', 'K', 'K', 'K', 'K', 'K']
const CARD = Float64[Inf, Inf, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0]
const NC = 2
const NK = 8
const MK = 2
const XS = Float64[0.21116886138916016, 0.7993184567995317]
const S = Float64[0.05582156087217378, -0.13922626670924634]
const L2S = Vector{Int}[Int[1, 0], Int[0, 1], Int[1, 0], Int[1, 0], Int[0, 1], Int[0, 1], Int[0, 1], Int[1, 0]]
const KS = Int[0, 1, 0, 0, 1, 0, 0, 1]
const LAMBDA = 0.0
const B = 0.36

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

function _objective(z)
    z = _vector(z)
    u = z[1:NC]
    any(v -> v < 0 || v > 1, u) && _err("continuous coordinates must lie in [0,1].")
    kk = z[NC+1:end]
    any(k -> k != floor(k) || k < 1 || k > MK, kk) && _err("label indices must be integers in 1..$MK.")
    h = 0
    for i in 1:NK
        h += L2S[i][Int(kk[i])] == KS[i] ? 1 : 0
    end
    f = 0.0
    for j in 1:NC
        d = u[j] - XS[j] - (h == 0 ? LAMBDA * S[j] : 0.0)
        f += d * d
    end
    return f + B * (1 - h / NK)
end

_problem_data() = Any[[[0.0, 1.0] for _ in 1:NC]..., [["k$(i)_$j" for j in 1:MK] for i in NC+1:N]...]

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
        cardinality = Float64[Inf, Inf, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0],
        f_star = 0.0,
        f_star_justification = "E",
        n_minimizers_documented = 1,
        description = "DFO Benchmark Suite v2.1.0 core, Family BK (controlled categorical dimension): 8 categorical variables x 2 categories.",
        x0_idx = Float64[0.780517578125, 0.4244779759183051, 2.0, 2.0, 2.0, 1.0, 1.0, 1.0, 1.0, 2.0],
        x0_alt_idx = Vector{Float64}[Float64[0.266357421875, 0.35451912818167963, 1.0, 2.0, 1.0, 2.0, 2.0, 2.0, 1.0, 1.0], Float64[0.524169921875, 0.3723517756439567, 2.0, 2.0, 1.0, 2.0, 1.0, 2.0, 2.0, 1.0], Float64[0.020263671875, 0.4340801707056851, 1.0, 2.0, 2.0, 2.0, 1.0, 1.0, 2.0, 1.0], Float64[0.981201171875, 0.3641213229690596, 1.0, 2.0, 2.0, 1.0, 2.0, 1.0, 1.0, 2.0], Float64[0.471435546875, 0.39292790733119953, 2.0, 2.0, 1.0, 1.0, 2.0, 2.0, 1.0, 2.0], Float64[0.729248046875, 0.4107605547934766, 1.0, 2.0, 1.0, 1.0, 1.0, 1.0, 2.0, 2.0], Float64[0.225341796875, 0.3408017070568511, 2.0, 2.0, 2.0, 1.0, 2.0, 1.0, 2.0, 1.0], Float64[0.858154296875, 0.4025301021185795, 2.0, 2.0, 1.0, 2.0, 2.0, 2.0, 1.0, 1.0], Float64[0.345458984375, 0.43133668648071943, 1.0, 2.0, 1.0, 2.0, 1.0, 1.0, 2.0, 1.0]],
        suite_release = "2.1.0",
        configuration = "interBK_nC2_K8x2",
        lambda = 0.0,
        sigma = 0.15,
        b = 0.36,
        q = 0.0,
        j_min = 0.0,
        obstacle_ratio = 0.0,
        seed = 3.0,
        relabel = 1.0,
        base_function = "familyBK",
        family = "BK",
        pattern = "BK",
        nK = 8,
        m = 2,
        lb_orig = Float64[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        ub_orig = Float64[1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
        lb_C_work = Float64[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        ub_C_work = Float64[1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
        ProblemData = _problem_data(),
        func_idx = _objective,
        func_F = xr -> _objective(_cell2idx(xr)),
        idx2cell = _idx2cell,
        cell2idx = _cell2idx
        )
end

end # module
