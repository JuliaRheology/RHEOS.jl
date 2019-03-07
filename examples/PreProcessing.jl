#!/usr/bin/env julia
using Revise
using RHEOS

filedir = "DataIncomplete.csv"

# repeated step loading generated with FractionalSLS([2.0, 0.5, 0.5, 0.7])
data = importdata(filedir; t_col = 2, σ_col = 1)

# data_smoothed = smooth(data, 5.0)

data_resampled = var_resample(data, :σ, 0.1; _mapback = false)

data_resampled = fixedresample(data,-4)
data_resampled = fixedresample(data,[1,-4],time_boundaries=[1,15,150])



# data_resampled = mapbackdata(data_resampled, data_raw)
# data_resampled = downsample(data_raw, [0.0, 40.0], [3])
