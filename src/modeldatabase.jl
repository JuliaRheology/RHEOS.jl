#!/usr/bin/env julia

#
# Some wrapper function for external maths functions
#
@inline function mittleff(α,z)
    α>0 ? mittlefforiginal(α,z) : 10000000000.
end

@inline function mittleff(α,β,z)
    α>0 ? mittlefforiginal(α,β,z) : 10000000000.
end


#
#  These are not useful at the moment
#
#invLaplace(f::Function, t::Vector{RheoFloat}) = InverseLaplace.talbotarr(f, t)
#invLaplace(f::Function, t::RheoFloat) = InverseLaplace.talbot(f, t)



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
