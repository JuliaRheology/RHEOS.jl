#!/usr/bin/env julia

# Spring-Pot Model
function G_springpot(t::T, params::Vector{T}) where T<:Real
    cᵦ, β = params

    G = cᵦ*t^(-β)/gamma(1 - β)
end
G_springpot(t::Vector{T}, params::Vector{T}) where T<:Real = G_springpot.(t, (params,))

function J_springpot(t::T, params::Vector{T}) where T<:Real
    cᵦ, β = params

    J = (t^β)/(cᵦ*gamma(1 + β))
end
J_springpot(t::Vector{T}, params::Vector{T}) where T<:Real = J_springpot.(t, (params,))

function Gp_springpot(ω::T, params::Vector{T}) where T<:Real
    cᵦ, β = params

    Gp = cᵦ*(ω^β)*cos(π*β/2)
end
Gp_springpot(ω::Vector{T}, params::Vector{T}) where T<:Real = Gp_springpot.(ω, (params,))

function Gpp_springpot(ω::T, params::Vector{T}) where T<:Real
    cᵦ, β = params

    Gpp = cᵦ*(ω^β)*sin(π*β/2)
end
Gp_springpot(ω::Vector{T}, params::Vector{T}) where T<:Real = Gp_springpot.(ω, (params,))

SpringPot() = RheologyModel(G_springpot, J_springpot, Gp_springpot, Gpp_springpot, [2.0, 0.5], ["model created with default parameters"])
SpringPot(params::Vector{T}) where T<:Real = RheologyModel(G_springpot, J_springpot, Gp_springpot, Gpp_springpot, params, ["model created by user with parameters $params"])

# Spring
function G_spring(t::T, params::Vector{T}) where T<:Real
    k, = params

    G = k
end
G_spring(t::Vector{T}, params::Vector{T}) where T<:Real = G_spring.(t, (params,))

function J_spring(t::T, params::Vector{T}) where T<:Real
    k, = params

    J = 1/k
end
J_spring(t::Vector{T}, params::Vector{T}) where T<:Real = J_spring.(t, (params,))

function Gp_spring(ω::T, params::Vector{T}) where T<:Real
    k, = params

    Gp = k
end
Gp_spring(ω::Vector{T}, params::Vector{T}) where T<:Real = Gp_spring.(ω, (params,))

function Gpp_spring(ω::T, params::Vector{T}) where T<:Real
    k, = params

    Gpp = 0.0
end
Gpp_spring(ω::Vector{T}, params::Vector{T}) where T<:Real = Gpp_spring.(ω, (params,))

Spring() = RheologyModel(G_spring, J_spring, Gp_spring, Gpp_spring, [1.0], ["model created with default parameters"])
Spring(params::Vector{T}) where T<:Real = RheologyModel(G_spring, J_spring, Gp_spring, Gpp_spring, params, ["model created by user with parameters $params"])

# Dashpot
function G_dashpot(t::T, params::Vector{T}) where T<:Real
    # undefined, use null-modulus instead
end
G_dashpot(t::Vector{T}, params::Vector{T}) where T<:Real = G_dashpot.(t, (params,))

function J_dashpot(t::T, params::Vector{T}) where T<:Real
    η, = params

    J = t/η
end
J_dashpot(t::Vector{T}, params::Vector{T}) where T<:Real = J_dashpot.(t, (params,))

function Gp_dashpot(ω::T, params::Vector{T}) where T<:Real
    η, = params

    Gp = 0.0
end
Gp_dashpot(ω::Vector{T}, params::Vector{T}) where T<:Real = Gp_dashpot.(ω, (params,))

function Gpp_dashpot(ω::T, params::Vector{T}) where T<:Real
    η, = params

    Gpp = η*ω
end
Gpp_dashpot(ω::Vector{T}, params::Vector{T}) where T<:Real = Gpp_dashpot.(ω, (params,))

DashPot() = RheologyModel(null_modulus, J_dashpot, Gp_dashpot, Gpp_dashpot, [1.0], ["model created with default parameters"])
DashPot(params::Vector{T}) where T<:Real = RheologyModel(null_modulus, J_dashpot, Gp_dashpot, Gpp_dashpot, params, ["model created by user with parameters $params"])
