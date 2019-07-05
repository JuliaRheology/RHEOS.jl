#!/usr/bin/env julia
using PyPlot
using RHEOS

frac_spec = RheoModel(FractionalSpecial, (η = 15e3, cᵦ = 1e3, β = 0.2, k = 300))
frequency = frequency_spec(ω_start = 1.e-3)
fracspecial_spect = dynamicmodelpredict(frequency,frac_spec)

p0 = (η = 1e4, cᵦ=1e3, β=0.30, k=3e2)
lb = (η = 1e4, cᵦ=8e2, β=0.20, k=2e2)
ub = (η = 2e4, cᵦ=1.2e3, β=0.40, k=4e2)
fit_frac_spec = dynamicmodelfit(fracspecial_spect,FractionalSpecial; weights ="log", p0 =p0, lo = lb, hi = ub, verbose = true)
predict_frac = dynamicmodelpredict(frequency,fit_frac_spec)

fig, ax = subplots()
ax[:loglog](fracspecial_spect.ω, fracspecial_spect.Gp, label="Data", color="blue", "o")
ax[:loglog](fracspecial_spect.ω, fracspecial_spect.Gpp, label="Data", color="red", "o")
ax[:loglog](predict_frac.ω, predict_frac.Gp, label="Fit", color="blue")
ax[:loglog](predict_frac.ω, predict_frac.Gpp, label="Fit", color="red")
