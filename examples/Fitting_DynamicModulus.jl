#!/usr/bin/env julia
using PyPlot
using Revise
using RHEOS

filedir = "FrequencyData.csv"

# ------------------ IMPORT DATA example -------------------------------

data = importdata(filedir; ω_col =3, Gp_col = 1, Gpp_col = 2)
data_ext = extract(data,0)
data_ext = extract(data,frec_only)

#------------------------------GENERATE DATA AND FIT----------------------------

frac_spec = RheoModel(FractionalSpecial, (η = 15e3, cᵦ = 1e3, β = 0.2, k = 300))
frequency = frequency_spec(ω_start = 1.e-3)
fracspecial_spect = dynamicmodelpredict(frequency,frac_spec)

p0 = (η = 1e4, cᵦ=1e3, β=0.30, k=3e2)
lb = (η = 1e4, cᵦ=8e2, β=0.20, k=2e2)
ub = (η = 2e4, cᵦ=1.2e3, β=0.40, k=4e2)
fit_log = dynamicmodelfit(fracspecial_spect,FractionalSpecial; weights ="log", p0 =p0, lo = lb, hi = ub, verbose = true)
predict_log = dynamicmodelpredict(frequency,fit_log)

fit_global = dynamicmodelfit(fracspecial_spect,FractionalSpecial; weights ="global", p0 =p0, lo = lb, hi = ub, verbose = true)
predict_global = dynamicmodelpredict(frequency,fit_global)

fit_linear = dynamicmodelfit(fracspecial_spect,FractionalSpecial; weights ="linear", p0 =p0, lo = lb, hi = ub, verbose = true)
predict_linear = dynamicmodelpredict(frequency,fit_linear)


fig, ax = subplots()
ax[:loglog](fracspecial_spect.ω, fracspecial_spect.Gp, label="Data", color="blue", "o")
ax[:loglog](fracspecial_spect.ω, fracspecial_spect.Gpp, label="Data", color="red", "o")
ax[:loglog](predict_frac.ω, predict_frac.Gp, label="Fit", color="blue")
ax[:loglog](predict_frac.ω, predict_frac.Gpp, label="Fit", color="red")
