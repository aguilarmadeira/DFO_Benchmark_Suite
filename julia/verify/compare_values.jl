# MATLAB-Julia check of the DFO Benchmark Suite (julia/verify).
#   julia compare_values.jl            compare values_julia.csv with the MATLAB values of ../../python/verify
#                                      (values_matlab.csv from eval_points.m if present, else values_matlab_reference.csv)
#   julia compare_values.jl --julia    first recompute values_julia.csv with the Julia files of ../core
# Exit code 1 when the point sets differ, when finiteness differs, or when a relative difference
# |j - m| / max(1, |m|) exceeds 1e-12.
const TOL = 1e-12
here = @__DIR__; pv = joinpath(here, "..", "..", "python", "verify")
if "--julia" in ARGS
    run(`$(Base.julia_cmd()) --startup-file=no $(joinpath(here, "eval_points.jl")) $(joinpath(here, "..", "core")) $(joinpath(pv, "points.csv")) $(joinpath(here, "values_julia.csv"))`)
end
function load(f)
    d = Dict{NTuple{3,String},Float64}()
    for l in eachline(f)
        c = split(l, ','); d[(c[1], c[2], c[3])] = parse(Float64, c[end])
    end
    return d
end
mfile = isfile(joinpath(pv, "values_matlab.csv")) ? joinpath(pv, "values_matlab.csv") : joinpath(pv, "values_matlab_reference.csv")
m = load(mfile); j = load(joinpath(here, "values_julia.csv"))
keys_pts = Set((c = split(l, ','); (c[1], c[2], c[3])) for l in eachline(joinpath(pv, "points.csv")))
fail = String[]
(Set(keys(m)) == keys_pts && Set(keys(j)) == keys_pts) || push!(fail, "point sets differ")
rows = Tuple{Float64,NTuple{3,String}}[]
for k in intersect(keys_pts, keys(m), keys(j))
    a, b = j[k], m[k]
    d = isfinite(a) != isfinite(b) ? Inf : (a == b ? 0.0 : (isfinite(a) ? abs(a - b) / max(1.0, abs(b)) : Inf))
    push!(rows, (d, k))
end
sort!(rows, rev = true)
nsame = count(r -> r[1] == 0, rows); nbad = count(r -> r[1] > TOL, rows)
println("$(basename(mfile)) vs values_julia.csv: $(length(rows)) points, identical $nsame ($(round(100nsame / length(rows), digits = 1))%), max relative difference $(rows[1][1]), above $TOL: $nbad")
for (d, k) in first(rows, 5); println("   $d  $(k[1])/$(k[2]) point $(k[3])"); end
nbad > 0 && push!(fail, "$nbad points above the tolerance")
if isempty(fail); println("PASS"); else; println("FAIL: ", join(fail, "; ")); exit(1); end
