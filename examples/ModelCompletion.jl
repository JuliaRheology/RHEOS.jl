#!/usr/bin/env julia
using RHEOS
using PyPlot

filedir = "DataIncomplete.csv"

data = importdata(["strain", "time"], filedir)

model = SLS([500.0, 100.0, 50.0])

predicted = modelpredict(data, model, :G)

fig, ax = subplots()
ax[:plot](predicted.t, predicted.σ, label="SLS")
ax[:legend](loc="best")
show()
