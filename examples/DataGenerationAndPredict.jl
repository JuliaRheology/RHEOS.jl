#!/usr/bin/env julia
include("../src/RHEOS.jl")
using RHEOS
using PyPlot

# get data
sine = oscillatordata(800.0, 1/50; phase = -90.0)
ramp = rampdata(800.0, 10.0, 400.0) - rampdata(825.0, 400.0, 800.0)
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
step = stepdata(50.0, 25.0; stepsize = 0.05, t_trans = 2.5)
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

######################################
#### TEMPORARY EXAMPLES FOR DEBUG ####
######################################

# get data
# step = stepdata(50.0, 25.0; stepsize = 0.05, t_trans = 2.5)
# repeatedstep = repeatdata(step, 5; t_trans = 2.5)

# data = RheologyData(repeatedstep)

# sls_model = RheologyModel(G_SLS)
# springpot_model = RheologyModel(G_springpot)
# fractKV_model = RheologyModel(G_fractKV)

# fractzener_model = RheologyModel(G_fractzener)

# sls_predicted = modelpredict(data, sls_model)
# springpot_predicted = modelpredict(data, springpot_model)
# fractzener_predicted = modelpredict(data, fractzener_model)
# fractKV_predicted = modelpredict(data, fractKV_model)

# fig, ax = subplots()
# ax[:plot](data.t, data.σ, "-", label="original")
# ax[:plot](sls_predicted.t, sls_predicted.σ, "-", label="SLS")
# ax[:plot](springpot_predicted.t, springpot_predicted.σ, "--", label="springpot")
# ax[:plot](fractKV_predicted.t, fractKV_predicted.σ, "--", label="fractKV")
# ax[:plot](fractzener_predicted.t, fractzener_predicted.σ, "--", label="fract zener")
# ax[:legend](loc="best")
# show()

# get data 
# stepUP = stepdata(150.0, 75.0; stepsize = 0.05, t_trans = 10.0)
# stepDOWN = stepdata(150.0, 75.0; amplitude = -1.0, stepsize = 0.05, t_trans = 10.0)

# dataUP = RheologyData(stepUP)
# dataDOWN = RheologyData(stepDOWN)

# springpot_model = RheologyModel(G_springpot)

# sp_predUP = modelpredict(dataUP, springpot_model)
# sp_predDOWN = modelpredict(dataDOWN, springpot_model)

# fig, ax = subplots()
# ax[:plot](dataUP.t, dataUP.ϵ, "-", label="original")
# ax[:plot](dataDOWN.t, dataDOWN.ϵ, "-", label="original")
# ax[:plot](sp_predUP.t, sp_predUP.ϵ, "-", label="springpot UP")
# ax[:plot](sp_predDOWN.t, sp_predDOWN.ϵ, "-", label="springpot DOWN")
# ax[:plot](sp_predDOWN.t, -sp_predDOWN.ϵ, "--", label="springpot DOWN flipped!")
# ax[:legend](loc="best")
# show()

# get data
# step = stepdata(50.0, 25.0; amplitude = +1.0, stepsize = 1.0, t_trans = 2.0)
# repeatedstep = repeatdata(step, 5; t_trans = 2.0)

# data = RheologyData(repeatedstep)

# springpot_model = RheologyModel(J_springpot)

# springpot_predicted = modelpredict(data, springpot_model)

# fig, ax = subplots()
# ax[:plot](data.t, data.ϵ, "-", label="original")
# # ax[:plot](springpot_predicted.t, springpot_predicted.ϵ, "--", label="springpot")
# ax[:legend](loc="best")
# show()

# get data
# stepUP = stepdata(50.0, 2.5; amplitude = +1.0, stepsize = 0.05, t_trans = 2.5)
# stepDOWN = stepdata(50.0, 27.5; amplitude = -2.0, stepsize = 0.05, t_trans = 2.5)
# step = stepUP + stepDOWN
# repeatedstep = repeatdata(step, 5; t_trans = 2.5)

# data = RheologyData(repeatedstep)

# springpot_model = RheologyModel(G_springpot)

# springpot_predicted = modelpredict(data, springpot_model)

# fig, ax = subplots()
# ax[:plot](data.t, data.ϵ, "-", label="original")
# ax[:plot](springpot_predicted.t, springpot_predicted.ϵ, "--", label="springpot")
# ax[:legend](loc="best")
# show()