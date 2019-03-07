#!/usr/bin/env julia
using PyPlot
using Revise
using RHEOS

filedir = "DataRelaxation.csv"

# repeated step loading generated with FractionalSLS([2.0, 0.5, 0.5, 0.7])
data = importdata(filedir; t_col =3, ϵ_col = 2, σ_col = 1)

# SLS fit
sls_fit = modelfit(data, SLS(); modtouse = :G, verbose=true)
sls_fit2 = modelfit(data, SLS());

# Spring-pot fit: cₐ, a, kᵦ, kᵧ
lb = [0.1, 0.01, 0.1, 0.1]
ub = [Inf, 0.99, Inf, Inf]
fractsls_fit = modelfit(data, FractionalSLS(); modtouse=:G, lo=lb, hi=ub, verbose=true)

# # get curves based on models fitted
sls_predicted = modelpredict(data, sls_fit; modtouse = :G)
sls_predicted = modelpredict(data, sls_fit2)
fractsls_predicted = modelpredict(data, fractsls_fit, :G)

# plot all data
fig, ax = subplots()
ax[:plot](data.t, data.σ, label="Data", color="black")
ax[:plot](sls_predicted.t, sls_predicted.σ, label="SLS")
ax[:plot](fractsls_predicted.t, fractsls_predicted.σ, "--", label="Fractional SLS")
ax[:legend](loc="best")
show()
