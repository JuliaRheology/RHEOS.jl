#!/usr/bin/env julia
using PyPlot
using RHEOS

filedir = "DataComplete.csv"

# repeated step loading generated with FractionalSLS([2.0, 0.5, 0.5, 0.7])
data = fileload(["stress","strain", "time"], filedir)

# SLS fit
sls_fit = modelfit(data, SLS(), :G)

# Spring-pot fit: cₐ, a, kᵦ, kᵧ
lb = [0.1, 0.01, 0.1, 0.1]
ub = [Inf, 0.99, Inf, Inf]
fractsls_fit = modelfit(data, FractionalSLS([2.0, 0.5, 0.5, 0.7]), :G; lo=lb, hi=ub, rel_tol=1e-5, verbose=true)

# # get curves based on models fitted
sls_predicted = modelpredict(data, sls_fit, :G)
fractsls_predicted = modelpredict(data, fractsls_fit, :G)

# plot all data
fig, ax = subplots()
ax[:plot](data.t, data.σ, label="Data", color="black")
ax[:plot](sls_predicted.t, sls_predicted.σ, label="SLS")
ax[:plot](fractsls_predicted.t, fractsls_predicted.σ, "--", label="Fractional SLS")
ax[:legend](loc="best")
show()