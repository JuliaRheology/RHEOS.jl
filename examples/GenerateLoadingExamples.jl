#!/usr/bin/env julia
# using RHEOS
# using PyPlot

# NEW EXAMPLES FROM ALE/ALEXANDRE? COPY-PASTED FROM SOURCE FILE
# timeline(end; start=0; stepsize=(end-start)/250)
# → RheolTimeData with time only defined

# strainfunction(RheolData, Function)
# stressfunction(RheolData, Function)

# d=timeline(100; stepsize=0.1)
# strainfunction!(d,t->step(t,ts=10)+rand())

# foo = stepgen(1000.0, 100.0)
# bar = stepgen(1000.0, 500.0)
# baz = foo - bar

# fig, ax = subplots()
# ax[:plot](foo.t, foo.ϵ, label="foo", linewidth=5, color="black")
# ax[:plot](bar.t, bar.ϵ, "--", label="bar", linewidth=3, color="orange")
# ax[:plot](baz.t, baz.ϵ, "-.", label="baz", linewidth=3.5, color="green")
# ax[:legend](loc="best")
# show()

# ramp test
foo = rampgen(1000.0, 100.0, 200.0)
bar = rampgen(1000.0, 500.0, 700.0)
baz = 2*foo - 2*bar
plot(baz.t, baz.ϵ)

# oscillator test
foo = singen(1000.0, 1/50; phase = -pi/2)
bar = rampgen(1000.0, 10.0, 400.0) - rampgen(1000.0, 400.0, 800.0)
baz = foo*bar
plot(baz.t, baz.ϵ)

# repeat test
foo = -stepgen(170.0, 125.0; stepsize = 0.001)
bar = repeatdata(foo, 5)
plot(bar.t, bar.ϵ, "--")
show()

# complicated test with noise

# generate a single step at 25 seconds
stepup = stepgen(50.0, 25.0; stepsize = 0.05, t_trans = 2.5)

# generate an oscillation which starts fading in at 25.5 seconds and has faded out by 49.5 seconds
osci = 0.1*singen(50.0, 0.2; stepsize = 0.05)
rampup = rampgen(50.0, 25.5, 37.5; stepsize = 0.05)
rampdown = -rampgen(50.0, 37.5, 49.5; stepsize = 0.05)

# combine the step and faded oscillation
combined = osci*(rampup + rampdown) + stepup

# repeat this three times
repeated = repeatdata(combined, 3)

# add some white noise with amplitude of 0.01
noisyrepeated = repeated + 0.01*noisegen(150.0; seed = 1, stepsize = 0.05)

plot(repeated.t, repeated.ϵ)
plot(noisyrepeated.t, noisyrepeated.ϵ, alpha = 0.7)
show()
