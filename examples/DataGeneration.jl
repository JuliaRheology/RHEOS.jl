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
foo = stepgen(170.0, 125.0; amplitude = -1.0, t_trans = 1.0)
bar = repeatdata(foo, 5; t_trans = 1.0)
plot(bar.t, bar.data, "--")
show()

# complicated test with noise
stepup = stepgen(50.0, 25.0; stepsize = 0.05, t_trans = 2.5)
osci = singen(50.0, 0.4; stepsize = 0.05, amplitude = 0.1)
rampup = rampgen(50.0, 25.0, 37.5; stepsize = 0.05)
rampdown = rampgen(50.0, 37.5, 48.0; stepsize = 0.05, amplitude = -1.0)
combined = osci*(rampup + rampdown) + stepup
repeated = repeatdata(combined, 3)
plot(repeated.t, repeated.data)
show()

## To add

# direct assignment
# foo = zerosgen(1000) (or flat gen, e.g. constantgen(1.0) for constant line at 1.0)

# add constant to all data elements
# baz = foo+1

# shift (and mirror?)
# baz = shiftgen, 100)
# baz = mirrorgen)