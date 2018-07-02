#!/usr/bin/env julia
include("../src/RHEOS.jl")
using RHEOS
using PyPlot

# step test
foo = stepdata(800.0, 100.0; baseval = 0.0, stepsize = 0.5)
bar = stepdata(850.0, 500.0; amplitude = -1.0)
baz = foo + bar
plot(baz.t, baz.data)

# ramp test
foo = rampdata(800.0, 100.0, 200.0; baseval = 0.0)
bar = rampdata(825.0, 500.0, 700.0; amplitude = -1.0)
baz = foo - bar
plot(baz.t, baz.data)

# oscillator test
foo = oscillatordata(800.0, 1/50; phase = -90.0)
bar = rampdata(800.0, 10.0, 400.0) - rampdata(825.0, 400.0, 800.0)
baz = foo*bar
plot(baz.t, baz.data)

# repeat test
foo = stepdata(170.0, 125.0; amplitude = -1.0, t_trans = 1.0)
bar = repeatdata(foo, 5; t_trans = 1.0)
plot(bar.t, bar.data, "--")
show()