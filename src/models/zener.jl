#!/usr/bin/env julia

# Zener Model (Standard Linear Solid in Maxwell Form)
function G_zener(t::Array{T,1}, params::Array{T,1}) where T<:Real
    k₀, k₁, η₁ = params

    G = k₀ + k₁*exp.(-t*k₁/η₁)
end

function J_zener(t::Array{T,1}, params::Array{T,1}) where T<:Real
    k₀, k₁, η₁ = params

    c₀ = 1/k₀
    c₁ = k₁/(k₀*(k₀ + k₁))
    τᵣ = η₁*(k₀ + k₁)/(k₀*k₁)

    J = c₀ - c₁*exp.(-t/τᵣ)
end

function Gp_zener(ω::Array{T,1}, params::Array{T,1}) where T<:Real

end

Zener() = RheologyModel(G_zener, J_zener, null_modulus, null_modulus, [1.0, 0.5, 1.0], ["model created with default parameters"])
Zener(params::Array{T, 1}) where T<:Real = RheologyModel(G_zener, J_zener, null_modulus, null_modulus, params, ["model created by user with parameters $params"])
