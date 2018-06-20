#!/usr/bin/env julia
include("../src/RHEOS.jl")
using RHEOS
using PyPlot

filedir = "../data/rheologyData1_incomplete.csv"

data_partial = fileload(["time", "stress"], filedir)

# data_resampled = var_resample(data_partial, :σ, 0.1; _mapback = false)
data_resampled = downsample(data_partial, [1, 450], [3])

println(data_resampled.sampling)
println(data_resampled.log)

plot(data_partial.t, data_partial.σ, "x")
plot(data_resampled.t, data_resampled.σ, "o")
show()

# # plot init
# fig, ax = subplots()

# plot stress data to check loaded in properly
# ax[:plot](data_partial.t, data_partial.σ)

# # get SLS data and plot strain data
# modelcomplete!(data_partial, "SLS", [844.152, 2043.98, 5.15527])
# ax[:plot](data_partial.t, data_partial.measured, label="SLS")

# # get springpot data and plot strain data
# modelcomplete!(data_partial, "springpot", [2997.7, 0.281836])
# ax[:plot](data_partial.t[2:end], data_partial.measured, label="springpot")

# # legend and plot
# ax[:legend](loc="best")
# show()
