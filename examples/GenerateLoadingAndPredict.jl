#!/usr/bin/env julia
using RHEOS
using PyPlot

# get data
sine = singen(800.0, 1/50; phase = -90.0)
ramp = rampgen(800.0, 10.0, 400.0) - rampgen(825.0, 400.0, 800.0)
data = ramp*sine

sls_predicted = modelpredict(data, SLS([844.152, 2043.98, 5.15527]), :G)
springpot_predicted = modelpredict(data, SpringPot([2997.7, 0.281836]), :G)

fig, ax = subplots()
ax[:plot](data.t, data.σ, "-", label="original")
ax[:plot](sls_predicted.t, sls_predicted.σ*1e-3, "-", label="SLS")
ax[:plot](springpot_predicted.t, springpot_predicted.σ*1e-3, "--", label="springpot")
ax[:legend](loc="best")
show()

# get data
step = stepgen(50.0, 25.0; stepsize = 0.01)
data = repeatdata(step, 3)

sls_predicted = modelpredict(data, SLS(), :G)
fractSLS_predicted = modelpredict(data, FractionalSLS([2.0, 0.5, 0.5, 0.7]), :G)

fig, ax = subplots()
ax[:plot](data.t, data.σ, "-", label="original")
ax[:plot](sls_predicted.t, sls_predicted.σ, "-", label="SLS")
ax[:plot](fractSLS_predicted.t, fractSLS_predicted.σ, "--", label="fractional SLS")
ax[:legend](loc="best")
show()

# get data
step = stepgen(50.0, 25.0; stepsize = 0.01)
data = repeatdata(step, 5)

springpot_predicted = modelpredict(data, SpringPot(), :J)
fractionalSpecial_predicted = modelpredict(data, FractionalSpecial(), :J)

fig, ax = subplots()
ax[:plot](springpot_predicted.t, springpot_predicted.ϵ, "--", label="springpot")
ax[:plot](fractionalSpecial_predicted.t, fractionalSpecial_predicted.ϵ, "--", label="fractional SLS")
ax[:legend](loc="best")
show()
