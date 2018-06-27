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

# Spring-pot fit
p0 = [1000.0, 0.5]
lb = [0.0, 0.0]
ub = [1e5, 1.0]
springpot_fit = modelfit(data_resampled, G_springpot; p0=p0, lo=lb, hi=ub, verbose=false)

savemodel(springpot_fit)

loaded_springpot = loadmodel("../data/rheologyData1.csv_RHEOS.G_springpot.jld")

sls_predicted = modelpredict(data_resampled, sls_fit)
springpot_predicted = modelpredict(data_resampled, loaded_springpot)

savedata(sls_predicted)

loaded_sls = loaddata(string(filedir, "_RheologyData.jld"))

plot(loaded_sls.t, loaded_sls.σ)
plot(springpot_predicted.t, springpot_predicted.σ)
plot(data_resampled.t, data_resampled.σ)
show()

