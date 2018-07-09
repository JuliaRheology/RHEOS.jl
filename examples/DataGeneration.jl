#!/usr/bin/env julia
include("../src/RHEOS.jl")
using RHEOS
using PyPlot

# step test
foo = stepgen(1000.0, 100.0)
bar = stepgen(1000.0, 500.0)
baz = foo - bar
plot(baz.t, baz.data)

# ramp test
foo = rampgen(1000.0, 100.0, 200.0)
bar = rampgen(1000.0, 500.0, 700.0)
baz = foo*2 - 2*bar
plot(baz.t, baz.data)

# oscillator test
foo = singen(1000.0, 1/50; phase = -pi/2)
bar = rampgen(1000.0, 10.0, 400.0) - rampgen(1000.0, 400.0, 800.0)
baz = foo*bar
plot(baz.t, baz.data)

# repeat test
foo = stepgen(170.0, 125.0; amplitude = -1.0, stepsize = 0.001)
bar = repeatdata(foo, 5)
plot(bar.t, bar.data, "--")
show()

# complicated test with noise
stepup = stepgen(50.0, 25.0; stepsize = 0.05, t_trans = 2.5)
osci = singen(50.0, 0.4; stepsize = 0.05, amplitude = 0.1)
rampup = rampgen(50.0, 25.0, 37.5; stepsize = 0.05)
rampdown = rampgen(50.0, 37.5, 48.0; stepsize = 0.05, amplitude = -1.0)
combined = osci*(rampup + rampdown) + stepup
repeated = repeatdata(combined, 3)
noisyrepeated = addnoise(repeated; amplitude = 0.05, seed = 1)
plot(repeated.t, repeated.data)
plot(noisyrepeated.t, noisyrepeated.data)
show()
