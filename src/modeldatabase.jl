#!/usr/bin/env julia

models_directory = joinpath(@__DIR__, "models")

# spring-pot, spring and dash-pot
include(joinpath(models_directory, "elements.jl"))

# fractional Maxwell model and specialized forms
include(joinpath(models_directory, "maxwell.jl"))

# fractional Kelvin-Voigt model and specialized forms
include(joinpath(models_directory, "kelvinvoigt.jl"))

# fractional Zener model and specialized forms (equivalent to Standard Linear Solid in Maxwell form)
include(joinpath(models_directory, "zener.jl"))

# poynting-thomson model and specialized forms (equivalent to Standard Linear Solid in Kelvin form)
include(joinpath(models_directory, "poynting-thomson.jl"))

# Four-parameters model
include(joinpath(models_directory, "burgers.jl"))

# Empirical
include(joinpath(models_directory, "empirical.jl"))

# Generalized
include(joinpath(models_directory, "generalized.jl"))
