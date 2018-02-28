#!/usr/bin/env julia
include("../src/RHEOS.jl")
using RHEOS

filedir = "../data/rheologyData1.csv"

data_raw = fileload(filedir, ["time","stress","strain"], "strlx")

#=
resampling examples
=#

data_resampled = fixed_resample(data_raw, [1,450],[8],["up"])
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
lb = [0.0, 0.0, 0.0]
ub = [1e5, 1e5, 1e5]
p0 = [1000.0, 1000.0, 100.0]
modelfit!(data_resampled, "SLS", p0, lb, ub)

# Spring-pot fit
lb = [0.0, 0.0]
ub = [1e5, 1.0]
p0 = [1000.0, 0.5]
modelfit!(data_resampled, "springpot", p0, lb, ub)

# Fract Special fit
lb = [0.0, 0.0, 0.02, 0.0]
ub = [1e3, 1e3, 0.98, 1e5]
p0 = [247.0, 6.48e2, 0.25, 4.26e3]
modelfit!(data_resampled, "fractspecial", p0, lb, ub)

fiteval(data_resampled, "SLS")
fiteval(data_resampled, "springpot")
fiteval(data_resampled, "fractspecial")
