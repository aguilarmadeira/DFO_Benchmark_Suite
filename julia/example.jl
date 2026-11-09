# Minimal use of the Julia version of the DFO Benchmark Suite: evaluate one instance and run the quasi-random (Halton)
# search baseline of benchmark/solver_random_search.m with the protocol of benchmark/PROTOCOL.md (x0 first, budget 1000n).
include(joinpath(@__DIR__, "core", "CDK_baseline", "ackley_cdk_10D.jl"))

function primes_first(k)
    p = Int[]; c = 2
    while length(p) < k
        all(q -> c % q != 0, filter(q -> q * q <= c, p)) && push!(p, c); c += 1
    end
    return p
end
function halton_point(idx, d)   # point idx (1-based) of the Halton sequence, bases = first d primes (as benchmark/halton_point.m)
    x = zeros(d)
    for (j, b) in enumerate(primes_first(d))
        i = idx; f = 1 / b; v = 0.0
        while i > 0
            v = v + f * (i % b); i = i ÷ b; f = f / b
        end
        x[j] = v
    end
    return x
end

PB = ackley_cdk_10D.problem(); n = PB.nvar; seed = 1; budget = 1000n
println(PB.name, " ", n, " ", PB.var_types, " f* = ", PB.f_star)
best = PB.func_idx(PB.x0_idx)                     # protocol: x0 is the first evaluation
for t in 1:budget-1
    u = halton_point(1048576 * seed + t, n)
    z = [PB.var_types[i] == "C" ? PB.lb_C_work[i] + u[i] * (PB.ub_C_work[i] - PB.lb_C_work[i]) :
         min(1 + floor(u[i] * PB.cardinality[i]), PB.cardinality[i]) for i in 1:n]
    global best = min(best, PB.func_idx(z))
end
println("best value after $budget evaluations: $best")
