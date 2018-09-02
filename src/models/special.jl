#!/usr/bin/env julia

# Fractional Special Model
function G_fractspecial(t::Array{T,1}, params::Array{T,1}) where T<:Real
    k, cᵦ, β, η = params

    G = k + cᵦ*t.^(-β).*mittleff.(1 - β, 1 - β, -cᵦ*(t.^(1 - β))/η)
end

function J_fractspecial(t::Array{T,1}, params::Array{T,1}) where T<:Real
    k, cᵦ, β, η = params

    a = η/cᵦ
    b = k*η/cᵦ

    J = InverseLaplace.ILt( s -> 1/s^2 * ((1+a*s^(1-β)) / (η+k/s+b*s^(-β)) ))

    return [J(t_val) for t_val in t]
end
                
FractionalSpecial() = RheologyModel(G_fractspecial, J_fractspecial, null_modulus, null_modulus, [1.0, 1.0, 0.2, 1.0], ["model created with default parameters"])
FractionalSpecial(params::Array{T, 1}) where T<:Real = RheologyModel(G_fractspecial, J_fractspecial, null_modulus, null_modulus, params, ["model created by user with parameters $params"])
