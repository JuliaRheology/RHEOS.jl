#!/usr/bin/env julia
include("../src/RHEOS.jl")
using RHEOS
using PyPlot

# get step data
data_step = stepdata_generate(1000.0, 250.0, 750.0, 10.0, 1000.0, "strlx")

# plot init
fig, ax = subplots()

# get SLS data and plot
modelcomplete!(data_step, "SLS", [844.152, 2043.98, 5.15527])
ax[:plot](data_step.t, data_step.measured, label="SLS")

# get springpot data and plot
modelcomplete!(data_step, "springpot", [2997.7, 0.281836])
ax[:plot](data_step.t[2:end], data_step.measured, label="springpot")

# legend and plot
ax[:legend](loc="best")
show()
