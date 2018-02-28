#!/usr/bin/env julia
include("../src/RHEOS.jl")
using RHEOS
using PyPlot

filedir = "../data/rheologyData1.csv"

data_raw = fileload(filedir, ["time","stress","strain"], "strlx")

data_partial = RheologyData(data_raw.ϵ, data_raw.t, "strlx")

modelcomplete!(data_partial, "SLS", [844.152, 2043.98, 5.15527])
plot(data_partial.t, data_partial.σ)

modelcomplete!(data_partial, "springpot", [2997.7, 0.281836])
plot(data_partial.t[2:end], data_partial.σ)

show()
#= think about extending what is stored in fittedmodels dict, could be
params, convolved with those params, cost and more. Then modelcomplete!
should place this convolved there, not in main self =#
