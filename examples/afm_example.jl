#!/usr/bin/env julia
include("../src/RHEOS.jl")
using RHEOS
using PyPlot

fdir = "../data/afmData1.txt"

data = AFMfileload(fdir, "strlx"; cpfind = "hertz", R = 25e-6)

data.insight = true

# data_resampled = var_resample(data, :δ, 0.05; _mapback = true)
# data_resampled = downsample(data, [1, length(data.t)], [10])
# data_resampled = smooth(data, 1.5; to_smooth = [:f])
# println(length(data.t))
# data_resampled = fixed_resample(data, [1, length(data.t)], [3], ["up"])
# println(length(data_resampled.t))
# data_resampled_mapped = mapbackdata(data_resampled, data)
# println(length(data_resampled_mapped.t))

p0 = [1000.0, 1000.0, 100.0]
lb = [0.0, 0.0, 0.0]
ub = [1e5, 1e5, 1e5]
# modelfit!(data_resampled, "SLS", p0, lb, ub)