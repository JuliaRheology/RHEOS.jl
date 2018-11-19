#!/usr/bin/env julia
# using PyPlot
# using RHEOS

# data_resampled = var_resample(data_raw, :σ, 0.1; _mapback = false)
# data_resampled = downsample(data_raw, [0.0, 40.0], [3])
# data_resampled = fixed_resample(data_raw, [0.0, 20.0, 41.0], [8, 25], ["up", "down"])
# data_resampled = fixed_resample(data_raw, [0.0, 40.0], [8], ["up"])
# data_resampled = smooth(data_raw, 5.0)
# data_resampled = mapbackdata(data_resampled, data_raw)