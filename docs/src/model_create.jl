# # Create Your Model
#md # [![](https://mybinder.org/badge_logo.svg)](@__BINDER_ROOT_URL__/notebooks/model_create.ipynb)
#md # [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/notebooks/model_create.ipynb)

# RHEOS provide an extensive model library (see Models section). However, if you would like to use a model that has not already been implemented in RHEOS, you can quickly define your own model as described below. Note that currently RHEOS can only handle linear viscoelastic models. Therefore, any added model must satisfy this condition. 

# The [`RheoModelClass`](@ref) constructor does not need to contain all moduli functions, depending on what variables are fitted and predicted down the line (relaxation modulus G for stress, creep function J for strain, and complex moduli Gp and Gpp for frequency domain data). As an example, we will assume that only the relaxation modulus is known.

using RHEOS
using Plots
#-

PowerLawEmpirical =  RheoModelClass(
                        ## Model name
                        name = "power_empirical",
                        ## Model parameters,
                        p = (:A, :α),
                        ## Relaxation modulus
                        G = quote
                                A*t^(-α)
                            end,
                        ## Network
                        info = "Empirical model"
                        )
#-

## generate timeline
dϵ = timeline(t_start = 1e-2, t_end = 1e1)
## calculates strain data by applying a function of time
dϵ = strainfunction(dϵ, hstep())

## Fix model parameters
powermodel = RheoModel(PowerLawEmpirical, A = 1, α = 0.8)

## Evaluate relaxation response
dpower = modelpredict(dϵ, powermodel)

# Create a log-log plot
p = plot(
    dpower.t, dpower.σ,
    xscale = :log10,
    yscale = :log10,
    xlabel = "Time",
    ylabel = "Stress",
    ylims = (0.1, 20),
    grid = true,
    legend = false,
    size = (400, 400), 
    framestyle = :box, 
    linewidth = 4
)
#!nb p #hide

# In a similar way, it is possible to define a model with only creep modulus. In this case, only strain can be fitted/predicted.
