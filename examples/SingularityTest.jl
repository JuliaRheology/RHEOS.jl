#!/usr/bin/env julia
using PyPlot
using RHEOS

filedir = "../data/rheologyData1.csv"

data_raw = fileload(["time","stress","strain"], filedir)

data_resampled = fixed_resample(data_raw, [0.0, 40.0], [8], ["up"])

lb = [0.0, 0.0]
ub = [Inf, 1.0]
springpot_fit = modelfit(data_resampled, SpringPot(), :G; lo=lb, hi=ub)

println(springpot_fit.parameters)