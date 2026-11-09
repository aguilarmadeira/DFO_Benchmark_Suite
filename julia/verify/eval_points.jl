# Evaluates the Julia files at the points of points.csv and writes values_julia.csv (same format as values_python.csv).
# Usage: julia eval_points.jl <core_julia> <points.csv> <values_julia.csv>
core, pts, out = ARGS
loaded = Dict{String,Any}()
open(out, "w") do io
    for l in eachline(pts)
        fo, name, j, z = split(l, ',')
        key = fo * "/" * name
        if !haskey(loaded, key)
            m = Module()
            Base.include(m, joinpath(core, fo, name * ".jl"))
            loaded[key] = Base.invokelatest(getfield(Base.invokelatest(getfield, m, Symbol(name)), :problem))
        end
        PB = loaded[key]
        x = [parse(Float64, s) for s in split(z, ';')]
        println(io, fo, ",", name, ",", j, ",", repr(Base.invokelatest(PB.func_idx, x)))
    end
end
