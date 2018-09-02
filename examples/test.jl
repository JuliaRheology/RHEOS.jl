using RHEOS
using PyPlot

t = collect(0.1:0.1:0.5)
println(t)

ϵ = [0.1 for i in t]
println(ϵ)

ϵdot = deriv(ϵ)
println(ϵdot)

params = [1.0]

out = boltzintegral_nonsing(G_spring, t, params, ϵdot)

plot(t, ϵdot)
plot(t, out)
show()