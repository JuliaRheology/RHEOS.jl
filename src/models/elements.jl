#!/usr/bin/env julia

# Spring-Pot Model
function G_springpot(t::Array{T,1}, params::Array{T,1}) where T<:Real
    cᵦ, β = params

    G = cᵦ*t.^(-β)/gamma(1 - β)
end

function J_springpot(t::Array{T,1}, params::Array{T,1}) where T<:Real
    cᵦ, β = params

    J = (t.^β)/(cᵦ*gamma(1 + β))
end

function Gp_springpot(ω::Array{T,1}, params::Array{T,1}) where T<:Real
    cᵦ, β = params

    Gp = cᵦ*(ω.^β)*cos(π*β/2)
end

function Gpp_springpot(ω::Array{T,1}, params::Array{T,1}) where T<:Real
    cᵦ, β = params

    Gpp = cᵦ*(ω.^β)*sin(π*β/2)
end

SpringPot() = RheologyModel(G_springpot, J_springpot, Gp_springpot, Gpp_springpot, [2.0, 0.5], ["model created with default parameters"])
SpringPot(params::Array{T, 1}) where T<:Real = RheologyModel(G_springpot, J_springpot, Gp_springpot, Gpp_springpot, params, ["model created by user with parameters $params"])

# Spring
function G_spring(t::Array{T,1}, params::Array{T,1}) where T<:Real
    k = params[1]

    G = [k for time in t]
end

function J_spring(t::Array{T,1}, params::Array{T,1}) where T<:Real
    k = params[1]

    J = [1/k for time in t]
end

function Gp_spring(ω::Array{T,1}, params::Array{T,1}) where T<:Real
    k = params[1]

    Gp = [k for frequency in ω]
end

function Gpp_spring(ω::Array{T,1}, params::Array{T,1}) where T<:Real
    k = params[1]

    Gpp = [0 for frequency in ω]
end

Spring() = RheologyModel(G_spring, J_spring, Gp_spring, Gpp_spring, [1.0], ["model created with default parameters"])
Spring(params::Array{T, 1}) where T<:Real = RheologyModel(G_spring, J_spring, Gp_spring, Gpp_spring, params, ["model created by user with parameters $params"])

# Dashpot
function G_dashpot(t::Array{T,1}, params::Array{T,1}) where T<:Real
    # undefined, use null-modulus instead
end

function J_dashpot(t::Array{T,1}, params::Array{T,1}) where T<:Real
    η = params[1]

    J = t/η
end

function Gp_dashpot(ω::Array{T,1}, params::Array{T,1}) where T<:Real
    η = params[1]

    Gp = [0 for frequency in ω]
end

function Gpp_dashpot(ω::Array{T,1}, params::Array{T,1}) where T<:Real
    η = params[1]

    Gpp = η*ω
end

DashPot() = RheologyModel(null_modulus, J_dashpot, Gp_dashpot, Gpp_dashpot, [1.0], ["model created with default parameters"])
DashPot(params::Array{T, 1}) where T<:Real = RheologyModel(null_modulus, J_dashpot, Gp_dashpot, Gpp_dashpot, params, ["model created by user with parameters $params"])
