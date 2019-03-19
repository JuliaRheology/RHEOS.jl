#!/usr/bin/env julia
export rheoconv,invLaplace,fun_gen

rheoconv(t::Real) = RheoFloat(t)
rheoconv(t::Array{T,1}) where T<:Real = convert(Vector{RheoFloat},t)

invLaplace(f::Function, t::Array{RheoFloat}) = InverseLaplace.talbotarr(f, t)
invLaplace(f::Function, t::RheoFloat) = InverseLaplace.talbot(f, t)


function fun_gen(;name,p,G,J,Gp,Gpp,info)
    parameters = begin
            string(join(string.(p), ","), "=params")
            end
    gs=Symbol("G_"*name)
    js=Symbol("J_"*name)
    gps=Symbol("Gp_"*name)
    gpps=Symbol("Gpp_"*name)
    mod_G = @eval $gs(t::Union{Array{RheoFloat,1},RheoFloat}, params::Array{RheoFloat,1}) = ($parameters,$G);
    mod_J = @eval $js(t::Union{Array{RheoFloat,1},RheoFloat}, params::Array{RheoFloat,1}) = ($parameters,$J);
    mod_Gp = @eval $gps(ω::Union{Array{RheoFloat,1},RheoFloat}, params::Array{RheoFloat,1}) = ($parameters,$Gp);
    mod_Gpp = @eval $gpps(ω::Union{Array{RheoFloat,1},RheoFloat}, params::Array{RheoFloat,1}) = ($parameters,$Gpp);
    param = p
    info_model = string("Model $name \n",info)
    @eval $gs(t::Union{Array{T1,1},T1}, params::Array{T2,1}) where {T1<:Real, T2<:Real} = $gs(rheoconv(t),rheoconv(params))
    @eval $js(t::Union{Array{T1,1},T1}, params::Array{T2,1}) where {T1<:Real, T2<:Real} = $js(rheoconv(t),rheoconv(params))
    return RheoModelClass(mod_G,mod_J,mod_Gp,mod_Gpp,param,info_model)
end



models_directory = joinpath(@__DIR__, "models")

# spring-pot, spring and dash-pot
include(joinpath(models_directory, "elements.jl"))

# fractional Maxwell model and specialized forms
include(joinpath(models_directory, "maxwell.jl"))

# fractional Kelvin-Voigt model and specialized forms
include(joinpath(models_directory, "kelvinvoigt.jl"))

# fractional Zener model and specialized forms (equivalent to Standard Linear Solid in Maxwell form)
include(joinpath(models_directory, "zener.jl"))

# fractional Special model
include(joinpath(models_directory, "special.jl"))

# poynting-thomson model and specialized forms (equivalent to Standard Linear Solid in Kelvin form)
include(joinpath(models_directory, "poynting-thomson.jl"))

# Standard Linear Solid model 2 time-scales
include(joinpath(models_directory, "sls2.jl"))

# Plateau-d Power Law
function G_platpow(t::Array{T,1}, params::Array{T,1}) where T<:Real
    Gᵩ, G₀, τ, α = params

    G = Gᵩ + (G₀ - Gᵩ)./(1 + t/τ).^(α)
end

PowerLawPlateau() = RheologyModel(G_platpow, null_modulus, null_modulus, null_modulus, [1.0, 2.0, 1.0, 0.2], ["model created with default parameters"])
PowerLawPlateau(params::Array{T, 1}) where T<:Real = RheologyModel(G_platpow, null_modulus, null_modulus, null_modulus, params, ["model created by user with parameters $params"])
