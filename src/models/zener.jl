#!/usr/bin/env julia

# Fractional Standard Linear Solid in (Fractional) Maxwell Form
function G_fracsls(t::Array{T,1}, params::Array{T,1}) where T<:Real
    k₀, k₁, cₐ, a = params

    G = k₀ + k₁*mittleff.(a, -(k₁/cₐ)*t.^a)
end

function J_fracsls(t::Array{T,1}, params::Array{T,1}) where T<:Real

end

function Gp_fracsls(ω::Array{T,1}, params::Array{T,1}) where T<:Real

end

function Gpp_fracsls(ω::Array{T,1}, params::Array{T,1}) where T<:Real

end

FractionalSLS() = RheologyModel(G_fracsls, null_modulus, null_modulus, null_modulus, [0.5, 0.5, 1.0, 0.2], ["model created with default parameters"])
FractionalSLS(params::Array{T, 1}) where T<:Real = RheologyModel(G_fracsls, null_modulus, null_modulus, null_modulus, params, ["model created by user with parameters $params"])

# Standard Linear Solid in Maxwell Form
function G_sls(t::Array{T,1}, params::Array{T,1}) where T<:Real
    k₀, k₁, η₁ = params

    G = k₀ + k₁*exp.(-t*k₁/η₁)
end

function J_sls(t::Array{T,1}, params::Array{T,1}) where T<:Real
    k₀, k₁, η₁ = params

    c₀ = 1/k₀
    c₁ = k₁/(k₀*(k₀ + k₁))
    τᵣ = η₁*(k₀ + k₁)/(k₀*k₁)

    J = c₀ - c₁*exp.(-t/τᵣ)
end

function Gp_sls(ω::Array{T,1}, params::Array{T,1}) where T<:Real

end

function Gpp_sls(ω::Array{T,1}, params::Array{T,1}) where T<:Real

end

SLS() = RheologyModel(G_sls, J_sls, null_modulus, null_modulus, [1.0, 0.5, 1.0], ["model created with default parameters"])
SLS(params::Array{T, 1}) where T<:Real = RheologyModel(G_sls, J_sls, null_modulus, null_modulus, params, ["model created by user with parameters $params"])
