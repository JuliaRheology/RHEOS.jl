#!/usr/bin/env julia
include("../src/RHEOS.jl")
using RHEOS
using PyPlot

# get data
sine = oscillatordata(800.0, 1/50; phase = -90.0)
ramp = rampdata(800.0, 10.0, 400.0) - rampdata(825.0, 400.0, 800.0)
rampedsine = ramp*sine

data = RheologyData(rampedsine)

sls_fit = RheologyModel(G_SLS, [844.152, 2043.98, 5.15527])
springpot_fit = RheologyModel(G_springpot, [2997.7, 0.281836])

sls_predicted = modelpredict(data, sls_fit)
springpot_predicted = modelpredict(data, springpot_fit)

fig, ax = subplots()
ax[:plot](data.t, data.σ, "-", label="original")
ax[:plot](sls_predicted.t, sls_predicted.σ*1e-3, "-", label="SLS")
ax[:plot](springpot_predicted.t, springpot_predicted.σ*1e-3, "--", label="springpot")
ax[:legend](loc="best")
show()

# get data
step = stepdata(50.0, 25.0; amplitude = 1.0)
repeatedstep = repeatdata(step, 5)

data = RheologyData(repeatedstep)

sls_fit = RheologyModel(G_SLS, [844.152, 2043.98, 5.15527])
springpot_fit = RheologyModel(G_springpot, [2997.7, 0.281836])

sls_predicted = modelpredict(data, sls_fit)
springpot_predicted = modelpredict(data, springpot_fit)

fig, ax = subplots()
ax[:plot](data.t, data.σ, "-", label="original")
ax[:plot](sls_predicted.t, sls_predicted.σ*1e-3, "-", label="SLS")
ax[:plot](springpot_predicted.t, springpot_predicted.σ*1e-3, "--", label="springpot")
ax[:legend](loc="best")
show()
