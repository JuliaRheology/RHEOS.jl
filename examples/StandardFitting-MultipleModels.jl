#!/usr/bin/env julia
include("../src/RHEOS.jl")
using PyPlot
using RHEOS

filedir = "../data/rheologyData1.csv"

data_raw = fileload(["time","stress","strain"], filedir)

# data_resampled = var_resample(data_raw, :σ, 0.1; _mapback = false)
# data_resampled = downsample(data_raw, [1, 450], [3])
# data_resampled = fixed_resample(data_raw, [1, 200, 450], [8, 25], ["up", "down"])
data_resampled = fixed_resample(data_raw, [1, 450], [8], ["up"])
# data_resampled = smooth(data_raw, 5.0)
# data_resampled = mapbackdata(data_resampled, data_raw)

# SLS fit
p0 = [1000.0, 1000.0, 100.0]
lb = [0.0, 0.0, 0.0]
ub = [1e5, 1e5, 1e5]
sls_fit = modelfit(data_resampled, G_SLS; p0=p0, lo=lb, hi=ub, verbose=false)
println(sls_fit.parameters)

# Spring-pot fit
p0 = [1000.0, 0.5]
lb = [0.0, 0.0]
ub = [1e5, 1.0]
springpot_fit = modelfit(data_resampled, G_springpot; p0=p0, lo=lb, hi=ub, verbose=false)
println(springpot_fit.parameters)

# save and load springpot model (RheologyModel type) for demonstration purposes (contains model name, params and log)
savemodel(springpot_fit)
loaded_springpot = loadmodel("../data/rheologyData1.csv_RHEOS.G_springpot.jld")

# get curves based on models fitted
sls_predicted = modelpredict(data_resampled, sls_fit)
springpot_predicted = modelpredict(data_resampled, loaded_springpot)

# # save and load RheologyData for demonstration purposes (contains σ, ϵ, t, sampling and log)
# savedata(sls_predicted)
# loaded_sls = loaddata(string(filedir, "_RheologyData.jld"))

# # plot all data
# plot(loaded_sls.t, loaded_sls.σ)
# plot(springpot_predicted.t, springpot_predicted.σ)
# plot(data_resampled.t, data_resampled.σ)
# show()

# finally, export SLS data as a CSV for example plotting in different software
exportdata(sls_predicted)