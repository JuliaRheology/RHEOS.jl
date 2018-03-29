#!/usr/bin/env julia
include("../src/RHEOS.jl")
using RHEOS
using PyPlot

fdir = "../data/afmData1.txt"

data = AFMfileload(fdir, "strlx"; cpfind = "hertz", R = 25e-6)

# data_resampled = var_resample(data, :δ, 0.05; _mapback = true)
# data_resampled = downsample(data, [1, length(data.t)], [10])
# data_resampled = smooth(data, 1.5; to_smooth = [:f])
# println(length(data.t))
# data_resampled = fixed_resample(data, [1, length(data.t)], [3], ["up"])
# println(length(data_resampled.t))
# data_resampled_mapped = mapbackdata(data_resampled, data)
# println(length(data_resampled_mapped.t))

data.insight = true

# p0 = [12363.7, 5201.04, 3.38386]
# lb = [12362.0, 5200.0, 3.0]
# ub = [12364.0, 5202.0, 3.5]
# modelfit!(data, "SLS", p0, lb, ub)
# fiteval(data, "SLS")

p0 = [12363.7, 0.5]
lb = [0.0, 0.01]
ub = [20000.0, 0.99]
modelfit!(data, "springpot", p0, lb, ub)
fiteval(data, "springpot")