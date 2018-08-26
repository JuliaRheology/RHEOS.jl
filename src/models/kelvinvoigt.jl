#!/usr/bin/env julia

# Fraction Kelvin-Voigt model
function G_fractKV(t::Array{T,1}, params::Array{T,1}) where T<:Real
    cₐ, a, cᵦ, β = params

    G = cₐ*t.^(-a)/gamma(1 - a) + cᵦ*t.^(-β)/gamma(1 - β) 
end

function J_fractKV(t::Array{T,1}, params::Array{T,1}) where T<:Real
    cₐ, a, cᵦ, β = params

    J = cₐ*t.^(a)*mittleff.(a - β, 1 + a, -cᵦ*t.^(a - β)/cₐ)
end

FractionalKelvinVoigt() = RheologyModel(G_fractKV, J_fractKV, null_modulus, null_modulus, [2.0, 0.2, 1.0, 0.5], ["model created with default parameters"])
FractionalKelvinVoigt(params::Array{T, 1}) where T<:Real = RheologyModel(G_fractKV, J_fractKV, null_modulus, null_modulus, params, ["model created by user with parameters $params"])
