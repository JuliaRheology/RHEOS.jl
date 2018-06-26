#!/usr/bin/env julia
include("../src/RHEOS.jl")
using RHEOS
using PyPlot

filedir = "../data/rheologyData1_incomplete.csv"

data_partial = fileload(["time", "strain"], filedir)

println(data_partial.t[5])
println(data_partial.σ[5])
println(data_partial.ϵ[5])

data_resampled = fixed_resample(data_partial, [1, 450], [8], ["up"])

sls_fit = RheologyModel(G_SLS, [843.149, 2024.2, 5.22901])

sls_predicted = modelpredict(data_resampled, sls_fit)

fig, ax = subplots()
ax[:plot](data_resampled.t, data_resampled.σ)
ax[:plot](sls_predicted.t, sls_predicted.σ, label="SLS")
ax[:legend](loc="best")
show()
