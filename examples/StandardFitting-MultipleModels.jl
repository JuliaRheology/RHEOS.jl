#!/usr/bin/env julia
using PyPlot
using Revise
using RHEOS

filedir = "Epi_relax.csv"

# repeated step loading generated with FractionalSLS([2.0, 0.5, 0.5, 0.7])
data = importdata(filedir; t_col =1, ϵ_col = 3, σ_col = 2)
data_resampled = fixedresample(data,[1],time_boundaries=[0.0,80])

# SLS fit
sls_fit = modelfit(data, SLS(); modtouse = :G, verbose=true)
sls_fit2 = modelfit(data_resampled, SLS());

# Spring-pot fit: cₐ, a, kᵦ, kᵧ
p0 = [1e4, 1e3, 0.3, 4e2]
lb = [0.1, 0.1, 0.1, 0.1]
ub = [Inf, Inf, 0.99, Inf]
fractspecial_fit = modelfit(data_resampled, FractionalSpecial(); p0 = p0, lo=lb, hi=ub, verbose=true)
fracspecial_predicted = modelpredict(data, fractspecial_fit)

fractspecial_fit = modelstepfit(data_resampled, FractionalSpecial(); p0 = p0, lo=lb, hi=ub, verbose=true)
fracspecial_predicted = modelpredict(data, fractspecial_fit)

# # get curves based on models fitted
sls_predicted = modelpredict(data, sls_fit; modtouse = :G)
sls_predicted = modelpredict(data, sls_fit2)
fractsls_predicted = modelpredict(data, fractsls_fit, :G)

# plot all data
fig, ax = subplots()
ax[:loglog](data.t, data.σ, label="Data", color="black")
#ax[:loglog](sls_predicted.t, sls_predicted.σ, label="SLS")
ax[:plot](fracspecial_predicted.t, fracspecial_predicted.σ, "--", label="Fractional SLS")
ax[:legend](loc="best")
show()
