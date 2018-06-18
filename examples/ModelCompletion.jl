#!/usr/bin/env julia
include("../src/RHEOS.jl")
using RHEOS
using PyPlot

filedir = "../data/rheologyData1.csv"

data_raw = fileload(filedir, ["time","stress","strain"], "strlx")

data_partial = RheologyData(data_raw.ϵ[1:450], data_raw.t[1:450], "strlx", filedir)

# plot init
fig, ax = subplots()

# get SLS data and plot
modelcomplete!(data_partial, "SLS", [844.152, 2043.98, 5.15527])
ax[:plot](data_partial.t, data_partial.measured, label="SLS")

# get springpot data and plot
modelcomplete!(data_partial, "springpot", [2997.7, 0.281836])
ax[:plot](data_partial.t[2:end], data_partial.measured, label="springpot")

# legend and plot
ax[:legend](loc="best")
show()
