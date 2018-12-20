#!/usr/bin/env julia

# Fractional Zener Model
function G_jeffreys_PT(t::Array{T,1}, params::Array{T,1}) where T<:Real

end

function J_jeffreys_PT(t::Array{T,1}, params::Array{T,1}) where T<:Real
    η₁, k, η₂ = params

    J = (- exp.(-k*t/η₁) .+ 1 )/k + t/η₂
end

function Gp_jeffreys_PT(ω::Array{T,1}, params::Array{T,1}) where T<:Real

end

function Gpp_jeffreys_PT(ω::Array{T,1}, params::Array{T,1}) where T<:Real

end

JeffreysPT() = RheologyModel(G_jeffreys_PT, J_jeffreys_PT, Gp_jeffreys_PT, Gpp_jeffreys_PT, [1.0, 1.0, 1.0], ["model created with default parameters"])
JeffreysPT(params::Array{T, 1}) where T<:Real = RheologyModel(G_jeffreys_PT, J_jeffreys_PT, Gp_jeffreys_PT, Gpp_jeffreys_PT, params, ["model created by user with parameters $params"])
