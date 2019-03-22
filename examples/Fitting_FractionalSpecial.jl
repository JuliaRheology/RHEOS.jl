#!/usr/bin/env julia
using PyPlot
using Revise
using RHEOS

filedir = "Epi_relax.csv"

# ------------------ IMPORT DATA AND RESAMPLING -------------------------------
# repeated step loading generated with FractionalSLS([2.0, 0.5, 0.5, 0.7])
data = importdata(filedir; t_col =1, σ_col = 2, ϵ_col = 3)
data_resampled_fix = fixedresample(data,-2,time_boundaries=[-0.6, 80.0])
data_resampled_var = fixedresample(data,[1,-12],time_boundaries=[0.0, 8.0,80.0])
data_cut = cutting(data,0,50.0)

# ----------------- EXTRACTING data for prediction -----------------------------
data_predict = extract(data,strain_only)
data_predict_fix = extract(data_resampled_fix,strain_only)
data_predict_var = extract(data_resampled_var,strain_only)

# -------------- FITTING ------------------------------------------------------
# Fractional model: η, cβ, β, k - Singularity - constant sampling
p0 = (η = 1e4, cᵦ=1e3, β=0.30, k=4e2)
lb = (η = 0.0, cᵦ=0.0, β=0.01, k=0.0)
ub = (η = Inf, cᵦ=Inf, β=0.99, k=Inf)
fractspecial_fit_fix = modelfit(data_resampled_fix, FractionalSpecial, strain_imposed; p0 = p0, lo=lb, hi=ub, verbose=true)
fracspecial_predicted_fix = modelpredict(data_predict, fractspecial_fit_fix)

# Fractional model - Singularity - variable sampling
fractspecial_fit_var = modelfit(data_resampled_var, FractionalSpecial, strain_imposed; p0 = p0, lo=lb, hi=ub, verbose=true)
fracspecial_predicted_var = modelpredict(data_predict_var, fractspecial_fit_var)

# Step Fit
fractspecial_fit_step = modelstepfit(data_resampled_fix, FractionalSpecial,strain_imposed; step = 0.3, p0 = p0)
fracspecial_predicted_step = modelsteppredict(data_predict, fractspecial_fit_step)

# plot all data
fig, ax = subplots()
ax[:loglog](data.t, data.σ, label="Data", color="black")
ax[:plot](fracspecial_predicted_fix.t, fracspecial_predicted_fix.σ, "--", label="Fractional special fix")
ax[:plot](fracspecial_predicted_var.t, fracspecial_predicted_var.σ, "--", label="Fractional special var")
ax[:plot](fracspecial_predicted_step.t, fracspecial_predicted_step.σ, "--", label="Fractional special step")
ax[:legend](loc="best")
show()

# ------------------ CREEP PREDICTION -----------------------------------------
# Generate step in stress
time = time_line(t_end = 200.0)
step_creep = stressfunction(time,hstep(amp=0.3))
# Predict creep
fracspecial_creep = modelpredict(step_creep, fractspecial_fit_step)
fracspecial_creep = modelsteppredict(step_creep, fractspecial_fit_step)

fig, ax = subplots()
ax[:loglog](fracspecial_creep.t, fracspecial_creep.ϵ, "--", label="Fractional special fix")

# ------------------ FREQUENCY RESPONSE PREDICTION------------------------------
frequency = frequency_spec(ω_start = 1.e-3)
fracspecial_spect = dynamicmodelpredict(frequency,fractspecial_fit_step)

fig, ax = subplots()
ax[:loglog](fracspecial_spect.ω, fracspecial_spect.Gp, "--", color = "blue", label="Fractional special fix")
ax[:loglog](fracspecial_spect.ω, fracspecial_spect.Gpp, "--", color = "red", label="Fractional special fix")
