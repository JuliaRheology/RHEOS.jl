#!/usr/bin/env julia
include("../src/RHEOS.jl")
using RHEOS
using PyPlot

# get data
sine = singen(800.0, 1/50; phase = -90.0)
ramp = rampgen(800.0, 10.0, 400.0) - rampgen(825.0, 400.0, 800.0)
rampedsine = ramp*sine

data = RheologyData(rampedsine)

sls_model = RheologyModel(G_SLS, [844.152, 2043.98, 5.15527])
springpot_model = RheologyModel(G_springpot, [2997.7, 0.281836])

sls_predicted = modelpredict(data, sls_model)
springpot_predicted = modelpredict(data, springpot_model)

fig, ax = subplots()
ax[:plot](data.t, data.σ, "-", label="original")
ax[:plot](sls_predicted.t, sls_predicted.σ*1e-3, "-", label="SLS")
ax[:plot](springpot_predicted.t, springpot_predicted.σ*1e-3, "--", label="springpot")
ax[:legend](loc="best")
show()

# get data
step = stepgen(50.0, 25.0; stepsize = 0.05, t_trans = 2.5)
repeatedstep = repeatdata(step, 5; t_trans = 2.5)

data = RheologyData(repeatedstep)

sls_model = RheologyModel(G_SLS)
springpot_model = RheologyModel(G_springpot)
fractzener_model = RheologyModel(G_fractzener)

sls_predicted = modelpredict(data, sls_model)
springpot_predicted = modelpredict(data, springpot_model)
fractzener_predicted = modelpredict(data, fractzener_model)

fig, ax = subplots()
ax[:plot](data.t, data.σ, "-", label="original")
ax[:plot](sls_predicted.t, sls_predicted.σ, "-", label="SLS")
ax[:plot](springpot_predicted.t, springpot_predicted.σ, "--", label="springpot")
ax[:plot](fractzener_predicted.t, fractzener_predicted.σ, "--", label="fract zener")
ax[:legend](loc="best")
show()
