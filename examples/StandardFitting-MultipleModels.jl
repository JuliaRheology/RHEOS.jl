#!/usr/bin/env julia
include("../src/RHEOS.jl")
using PyPlot
using RHEOS

filedir = "../data/rheologyData1.csv"

data_raw = fileload(filedir, ["time","stress","strain"])

println(data_raw.log)
println(data_raw.sampling)
plot(data_raw.t, data_raw.σ)
plot(data_raw.t, data_raw.ϵ)
show()

data_partial = constructRheologyData(["stress", "time"], data_raw.σ, data_raw.t)
println(data_raw.log)
println(data_raw.sampling)
plot(data_raw.t, data_raw.σ)
show()

data_partial = constructRheologyData(["strain", "time"], data_raw.ϵ, data_raw.t)
println(data_raw.log)
println(data_raw.sampling)
plot(data_raw.t, data_raw.ϵ)
show()

# data_resampled = fixed_resample(data_raw, [1, 450],[8],["up"])
# data_resampled = smooth(data_resampled, 5.0)
# data_resampled = var_resample(data_resampled, :σ, 0.1, _mapback = false)
# data_resampled = mapbackdata(data_resampled, data_raw)

# # SLS fit
# p0 = [1000.0, 1000.0, 100.0]
# lb = [0.0, 0.0, 0.0]
# ub = [1e5, 1e5, 1e5]
# modelfit!(data_resampled, "SLS", p0, lb, ub)

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

# fiteval(data_resampled)

# saveresult(data_resampled; include_data = true)
# saveresult(data_resampled; include_data = false)

# alldatasaved = loadresult("../data/rheologyData1_RheologyData.jld")
# allmetasaved = loadresult("../data/rheologyData1_RheologyMetadata.jld")

# println(" ")
# for n in fieldnames(data_resampled)

#     println(String(n), " field is same as datasaved: ", getfield(data_resampled, n)==getfield(alldatasaved, n))
#     println(String(n), " field is same as METAsaved: ", getfield(data_resampled, n)==getfield(allmetasaved, n))

#     if getfield(data_resampled, n)!=getfield(allmetasaved, n)
#         println(String(n), " field = ", getfield(allmetasaved, n))
#     end

#     println(" ")

# end
