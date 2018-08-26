#!/usr/bin/env julia

# Fractional Maxwell Model
function G_fractmaxwell(t::Array{T,1}, params::Array{T,1}) where T<:Real
    cₐ, a, cᵦ, β = params

    G = cᵦ*t.^(-β)*mittleff.(a - β, 1 - β, -cᵦ*t.^(a - β)/cₐ)
end

function J_fractmaxwell(t::Array{T,1}, params::Array{T,1}) where T<:Real
    cₐ, a, cᵦ, β = params

    J = t.^(a)/(cₐ*gamma(1 + a)) + t.^(β)/(cᵦ*gamma(1 + β))
end

FractionalMaxwell() = RheologyModel(G_fractmaxwell, J_fractmaxwell, null_modulus, null_modulus, [2.0, 0.2, 1.0, 0.5], ["model created with default parameters"])
FractionalMaxwell(params::Array{T, 1}) where T<:Real = RheologyModel(G_fractmaxwell, J_fractmaxwell, null_modulus, null_modulus, params, ["model created by user with parameters $params"])
