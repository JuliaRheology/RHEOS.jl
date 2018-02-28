#!/usr/bin/env julia
include("../src/RHEOS.jl")
using RHEOS

# load in data
filedir = "../data/rheologyData1.csv"
data_raw = fileload(filedir, ["time","stress","strain"], "strlx")

# resample
# data_resampled = fixed_resample(data_raw, [1,350],[8],["up"])
data_resampled = fixed_resample(data_raw, [1,350,550],[8,2],["up","down"])
println(data_resampled.sampling)

# embed previously found fit parameters
# data_resampled.fittedmodels["SLS"] = [844.152, 2043.98, 5.15527]
data_resampled.fittedmodels["springpot"] = [2997.7, 0.281836]


# check plots
# fiteval(data_resampled, "SLS")
fiteval(data_resampled, "springpot")
