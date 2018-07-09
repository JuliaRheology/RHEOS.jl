#!/usr/bin/env julia
using RHEOS
using PyPlot

filedir = "../data/rheologyData1_incomplete.csv"

data_partial = fileload(["time", "strain"], filedir)

data_resampled = fixed_resample(data_partial, [0.0, 40.0], [8], ["up"])

sls_fit = RheologyModel(G_SLS, [843.149, 2024.2, 5.22901])

sls_predicted = modelpredict(data_resampled, sls_fit)

fig, ax = subplots()
ax[:plot](sls_predicted.t, sls_predicted.σ, "--", label="SLS")
ax[:legend](loc="best")
show()
