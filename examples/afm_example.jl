#!/usr/bin/env julia
include("../src/RHEOS.jl")
using RHEOS
using PyPlot

fdir = "../data/afmData1.txt"

data = AFMfileload(fdir, "strlx"; cpfind = "hertz", R = 25e-6)

# data.insight = true

data_resampled = fixed_resample(data, [1,length(data.t)], [10], ["down"])

p0 = [1000.0, 1000.0, 100.0]
lb = [0.0, 0.0, 0.0]
ub = [1e5, 1e5, 1e5]
modelfit!(data_resampled, "SLS", p0, lb, ub)