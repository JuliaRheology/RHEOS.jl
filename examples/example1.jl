#!/usr/bin/env julia
include("../src/RHEOS.jl")
using RHEOS

filedir = "../data/rheologyData1.csv"

data_raw = fileload(filedir, ["time","stress","strain"], "strlx")

#=
resampling examples
=#

data_beat = fixed_resample(data_raw, [1,450],[8],["up"])
data_resampled = fixed_resample(data_beat, [1,length(data_beat.t)],[2],["down"])
data_resampled.insight = true

# data_smoothed = smooth(data_resampled, 5.0)
# data_varresampled = var_resample(data_smoothed, :σ, 0.1, _mapback = false)
#
# data_varresampled.insight = true
#
# data_mappedback = mapbackdata(data_varresampled, data_resampled)
# data_mappedback2 = mapbackdata(data_varresampled, data_raw)

#=
fitting examples
=#

# SLS fit
p0 = [1000.0, 1000.0, 100.0]
lb = [0.0, 0.0, 0.0]
ub = [1e5, 1e5, 1e5]
modelfit!(data_resampled, "SLS", p0, lb, ub)

# # Spring-pot fit
# p0 = [1000.0, 0.5]
# lb = [0.0, 0.0]
# ub = [1e5, 1.0]
# modelfit!(data_resampled, "springpot", p0, lb, ub)

# # Fract Special fit
# p0 = [247.0, 6.48e2, 0.25, 4.26e3]
# lb = [0.0, 0.0, 0.02, 0.0]
# ub = [1e3, 1e4, 0.98, 1e5]
# modelfit!(data_resampled, "fractspecial", p0, lb, ub)

# fiteval(data_resampled, "SLS")
# fiteval(data_resampled, "springpot")
# fiteval(data_resampled, "fractspecial")

saveresult(data_resampled; include_data = true)
saveresult(data_resampled; include_data = false)

alldatasaved = loadresult("../data/rheologyData1_RheologyData.jld")
allmetasaved = loadresult("../data/rheologyData1_RheologyMetadata.jld")

println(" ")
for n in fieldnames(data_resampled)

    println(String(n), " field is same as datasaved: ", getfield(data_resampled, n)==getfield(alldatasaved, n))
    println(String(n), " field is same as METAsaved: ", getfield(data_resampled, n)==getfield(allmetasaved, n))

    if getfield(data_resampled, n)!=getfield(allmetasaved, n)
        println(String(n), " field = ", getfield(allmetasaved, n))
    end

    println(" ")

end