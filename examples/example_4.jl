#!/usr/bin/env julia
include("../src/RHEOS.jl")
using RHEOS
using PyPlot

# get step data
data_step = stepdata_generate(250.0, 0.0, 250.0, 10.0, 170.0, "creep",step_size = 0.1)
plot(data_step.t,data_step.σ)

data_step.σ[1:end] = 0.0;
data_step.σ[1:40] = 4;
  data_step.σ[41] = 92;
  data_step.σ[42] = 52;
  data_step.σ[43] = -15;

  data_step.dcontrolled[1:end] = 0.0;
  data_step.dcontrolled[1] = 4/0.1;
    data_step.dcontrolled[41] = 92/0.1;
    data_step.dcontrolled[42] = 52/0.1;
    data_step.dcontrolled[43] = -15/0.1;



# plot init
fig, ax = subplots()

# get SLS data and plot
modelcomplete!(data_step, "SLS", [844.152, 2043.98, 5.15527])
ax[:plot](data_step.t, data_step.measured, label="SLS")

# get springpot data and plot
modelcomplete!(data_step, "springpot", [2997.7, 0.281836])
ax[:loglog](data_step.t[1:end], data_step.measured, label="springpot")

# get fractspecial data and plot
tic()
modelcomplete!(data_step, "fractspecial_slow", [840.0,1.7e3 , 0.28, 1.46e4])
a = toc()
ax[:loglog](data_step.t[2:end]-data_step.t[44], data_step.measured, label="fractspecial_slow")

tic()
modelcomplete!(data_step, "fractspecial_fast", [840.0,1.7e3 , 0.28, 1.46e4])
toc()
ax[:loglog](data_step.t[2:end]-data_step.t[44], data_step.measured, label="fractspecial_fast")

# legend and plot
ax[:legend](loc="best")
show()
