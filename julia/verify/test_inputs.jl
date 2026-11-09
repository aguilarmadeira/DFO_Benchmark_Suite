# Invalid-input test of all Julia files (as test_inputs.py for Python). Usage: julia test_inputs.jl <core_julia>
core = ARGS[1]; ninst = 0; nD = 0; checks = 0; fails = Any[]
function expect_error(f, a)
    global checks += 1
    try Base.invokelatest(f, a); catch e; return e isa ArgumentError; end
    return false
end
for fo in sort(readdir(core)), f in sort(filter(x -> endswith(x, ".jl"), readdir(joinpath(core, fo))))
    name = f[1:end-3]; m = Module(); Base.include(m, joinpath(core, fo, f))
    PB = Base.invokelatest(getfield(Base.invokelatest(getfield, m, Symbol(name)), :problem)); global ninst += 1
    il = Base.invokelatest; z0 = collect(PB.x0_idx); xr = il(PB.idx2cell, z0)
    (il(PB.func_F, xr) == il(PB.func_idx, z0) && il(PB.cell2idx, xr) == z0) || push!(fails, (name, "roundtrip"))
    hasD = false
    for (i, t) in enumerate(PB.var_types)
        bads = t == "K" ? Any[NaN, "zz", nothing, 3] : Any[NaN, Inf, "a", nothing, true]
        for b in bads
            x = Any[xr...]; x[i] = b
            expect_error(PB.func_F, x) || push!(fails, (name, t, i, b, "func_F accepted"))
        end
        z = copy(z0); z[i] = NaN
        expect_error(PB.func_idx, z) || push!(fails, (name, t, i, "func_idx accepted NaN"))
        if t != "C"
            for b in (0.0, PB.cardinality[i] + 1, 1.5)
                z = copy(z0); z[i] = b
                expect_error(PB.func_idx, z) || push!(fails, (name, t, i, b, "bad index accepted"))
            end
        end
        t == "D" && (hasD = true)
    end
    global nD += hasD
    expect_error(PB.func_idx, z0[1:end-1]) || push!(fails, (name, "short vector accepted"))
end
println("$ninst instances ($nD with D variables), $checks invalid-input checks, failures: $(length(fails))")
for x in first(fails, 10); println("   ", x); end
exit(isempty(fails) ? 0 : 1)
