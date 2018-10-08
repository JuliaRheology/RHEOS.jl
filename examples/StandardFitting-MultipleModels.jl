#!/usr/bin/env julia
# using PyPlot
# using RHEOS

filedir = "../data/rheologyData1.csv"

data_raw = fileload(["time","stress","strain"], filedir)

# data_resampled = var_resample(data_raw, :σ, 0.1; _mapback = false)
# data_resampled = downsample(data_raw, [0.0, 40.0], [3])
# data_resampled = fixed_resample(data_raw, [0.0, 20.0, 41.0], [8, 25], ["up", "down"])
data_resampled = fixed_resample(data_raw, [0.0, 40.0], [8], ["up"])
# data_resampled = smooth(data_raw, 5.0)
# data_resampled = mapbackdata(data_resampled, data_raw)

# SLS fit
sls_fit = modelfit(data_resampled, SLS(), :G)

# Spring-pot fit
lb = [0.0, 0.0]
ub = [Inf, 1.0]
springpot_fit = modelfit(data_resampled, SpringPot(), :G; lo=lb, hi=ub)

# # get curves based on models fitted
sls_predicted = modelpredict(data_resampled, sls_fit, :G)
springpot_predicted = modelpredict(data_resampled, springpot_fit, :G)

# plot all data
plot(data_resampled.t, data_resampled.σ)
plot(sls_predicted.t, sls_predicted.σ)
plot(springpot_predicted.t, springpot_predicted.σ, "--")
show()