#!/usr/bin/env julia
using Revise
using RHEOS
using PyPlot

filedir = "DataIncomplete.csv"

data = importdata(filedir; t_col =2, ϵ_col = 1)

model = SLS([500.0, 100.0, 50.0])

predicted = modelpredict(data, model)

fig, ax = subplots()
ax[:plot](predicted.t, predicted.σ, label="SLS")
ax[:legend](loc="best")
show()
