#!/usr/bin/env julia
using PyPlot
using Revise
using RHEOS

filedir = "Epi_relax.csv"

# repeated step loading generated with FractionalSLS([2.0, 0.5, 0.5, 0.7])
data = importdata(filedir; t_col =1, σ_col = 2, ϵ_col = 3)
data_resampled_fix = fixedresample(data,-2,time_boundaries=[-0.6, 80.0])
data_resampled_var = fixedresample(data,[1,-12],time_boundaries=[0.0, 8.0,80.0])

# variable for prediction
data_predict = importdata(filedir; t_col =1, ϵ_col = 3)
data_predict_fix = fixedresample(data_predict,-2,time_boundaries=[-0.6, 80.0])
data_predict_var = fixedresample(data_predict,[1,-12],time_boundaries=[0.0, 8.0,80.0])

# SLS fit - No singularity - constant sampling
sls_fit_fix = modelfit(data_resampled_fix, SLS(), strain_imposed, verbose=true)
# OR sls_fit_fix2 = modelfit(data_resampled_fix, SLS(), 1);
sls_predicted_fix = modelpredict(data_predict, sls_fit_fix)

# SLS fit - No singularity - variable sampling
sls_fit_var = modelfit(data_resampled_var, SLS(), strain_imposed, verbose=true)
sls_predicted_var = modelpredict(data_predict_var, sls_fit_var)

# Fractional model: η, cβ, β, k - Singularity - constant sampling
p0 = [1e4, 1e3, 0.3, 4e2]
lb = [0.1, 0.1, 0.1, 0.1]
ub = [Inf, Inf, 0.99, Inf]
fractspecial_fit_fix = modelfit(data_resampled_fix, FractionalSpecial(), strain_imposed; p0 = p0, lo=lb, hi=ub, verbose=true)
fracspecial_predicted_fix = modelpredict(data_predict, fractspecial_fit_fix)

# Fractional model - Singularity - variable sampling
fractspecial_fit_var = modelfit(data_resampled_var, FractionalSpecial(), strain_imposed; p0 = p0, lo=lb, hi=ub, verbose=true)
fracspecial_predicted_var = modelpredict(data_predict_var, fractspecial_fit_var)


# Step Fit
fractspecial_fit_step = modelstepfit(data_resampled_fix, FractionalSpecial(),strain_imposed; step = 0.3)
fracspecial_predicted_step = modelsteppredict(data_predict, fractspecial_fit_step)

# plot all data
fig, ax = subplots()
ax[:loglog](data.t, data.σ, label="Data", color="black")
ax[:loglog](sls_predicted_fix.t, sls_predicted_fix.σ, label="SLS_fix")
ax[:loglog](sls_predicted_var.t, sls_predicted_var.σ, label="SLS_var")
ax[:plot](fracspecial_predicted_fix.t, fracspecial_predicted_fix.σ, "--", label="Fractional special fix")
ax[:plot](fracspecial_predicted_var.t, fracspecial_predicted_var.σ, "--", label="Fractional special var")
ax[:plot](fracspecial_predicted_step.t, fracspecial_predicted_step.σ, "--", label="Fractional special step")
ax[:legend](loc="best")
show()


prova(strain_imposed)
