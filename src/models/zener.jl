#!/usr/bin/env julia

# Fractional Zener Model
function G_fraczener(t::T, params::Vector{T}) where T<:Real
    cₐ, a, cᵦ, β, cᵧ, γ = params

    G = cᵦ*t^(-β)*mittleff(a - β, 1 - β, -cᵦ*t^(a - β)/cₐ) + cᵧ*t^(-γ)/gamma(1 - γ)
end
G_fraczener(t::Vector{T}, params::Vector{T}) where T<:Real = G_fraczener.(t, (params,))

function J_fraczener(t::Vector{T}, params::Vector{T}) where T<:Real
    cₐ, a, cᵦ, β, cᵧ, γ = params

    Ĵ(s) = (1/s)*(cₐ*s^a + cᵦ*s^β)/(cₐ*s^a*cᵦ*s^β + cᵧ*s^γ*(cₐ*s^a + cᵦ*s^β))

    # method gives NaN at time 0.0
    # J = InverseLaplace.ILt(s -> Jbar(s))
    # return [J(t_val) for t_val in t]

    return InverseLaplace.talbotarr(s -> Ĵ(s), t)
end

function Gp_fraczener(ω::T, params::Vector{T}) where T<:Real
    cₐ, a, cᵦ, β, cᵧ, γ = params

    denominator = (cₐ*ω^a)^2 + (cᵦ*ω^β)^2 + 2*(cₐ*ω^a)*(cᵦ*ω^β)*cos((a - β)*π/2)

    numerator = ((cᵦ*ω^β)^2)*(cₐ*ω^a)*cos(a*π/2) + ((cₐ*ω^a)^2)*(cᵦ*ω^β)*cos(β*π/2)

    Gp = numerator/denominator + cᵧ*ω^γ*cos(γ*π/2)
end
Gp_fraczener(ω::Vector{T}, params::Vector{T}) where T<:Real = Gp_fraczener.(ω, (params,))

function Gpp_fraczener(ω::T, params::Vector{T}) where T<:Real
    cₐ, a, cᵦ, β, cᵧ, γ = params

    denominator = (cₐ*ω^a)^2 + (cᵦ*ω^β)^2 + 2*(cₐ*ω^a)*(cᵦ*ω^β)*cos((a - β)*π/2)

    numerator = ((cᵦ*ω^β)^2)*(cₐ*ω^a)*sin(a*π/2) + ((cₐ*ω^a)^2)*(cᵦ*ω^β)*sin(β*π/2)

    Gpp = numerator/denominator + cᵧ*ω^γ*sin(γ*π/2)
end
Gpp_fraczener(ω::Vector{T}, params::Vector{T}) where T<:Real = Gpp_fraczener.(ω, (params,))

FractionalZener() = RheologyModel(G_fraczener, J_fraczener, Gp_fraczener, Gpp_fraczener, [1.0, 0.7, 1.0, 0.5, 1.0, 0.2], ["model created with default parameters"])
FractionalZener(params::Vector{T}) where T<:Real = RheologyModel(G_fraczener, J_fraczener, Gp_fraczener, Gpp_fraczener, params, ["model created by user with parameters $params"])

# Fractional Standard Linear Solid in (Fractional) Maxwell Form
function G_fracsls(t::T, params::Vector{T}) where T<:Real
    cₐ, a, kᵦ, kᵧ = params

    G = kᵦ*mittleff(a, -(kᵦ/cₐ)*t^a) + kᵧ
end
G_fracsls(t::Vector{T}, params::Vector{T}) where T<:Real = G_fracsls.(t, (params,))

function J_fracsls(t::Vector{T}, params::Vector{T}) where T<:Real
    cₐ, a, kᵦ, kᵧ = params

    Ĵ(s) = (1/s)*(cₐ*s^a + kᵦ)/(cₐ*s^a*kᵦ + kᵧ*(cₐ*s^a + kᵦ))

    # method gives NaN at time 0.0
    # J = InverseLaplace.ILt(s -> Ĵ(s))
    # return [J(t_val) for t_val in t]

    return InverseLaplace.talbotarr(s -> Ĵ(s), t)
end

function Gp_fracsls(ω::T, params::Vector{T}) where T<:Real
    cₐ, a, kᵦ, kᵧ = params

    denominator = (cₐ*ω^a)^2 + kᵦ^2 + 2*(cₐ*ω^a)*kᵦ*cos(a*π/2)

    numerator = kᵦ^2*(cₐ*ω^a)*cos(a*π/2) + ((cₐ*ω^a)^2)*kᵦ

    Gp = numerator/denominator + kᵧ
end
Gp_fracsls(ω::Vector{T}, params::Vector{T}) where T<:Real = Gp_fracsls.(ω, (params,))

function Gpp_fracsls(ω::T, params::Vector{T}) where T<:Real
    cₐ, a, kᵦ, kᵧ = params

    denominator = (cₐ*ω^a)^2 + kᵦ^2 + 2*(cₐ*ω^a)*kᵦ*cos(a*π/2)

    numerator = kᵦ^2*(cₐ*ω^a)*sin(a*π/2)

    Gpp = numerator/denominator
end
Gpp_fracsls(ω::Vector{T}, params::Vector{T}) where T<:Real = Gpp_fracsls.(ω, (params,))

FractionalSLS() = RheologyModel(G_fracsls, J_fracsls, Gp_fracsls, Gpp_fracsls, [0.5, 0.5, 1.0, 0.2], ["model created with default parameters"])
FractionalSLS(params::Vector{T}) where T<:Real = RheologyModel(G_fracsls, J_fracsls, Gp_fracsls, Gpp_fracsls, params, ["model created by user with parameters $params"])

# Standard Linear Solid in Maxwell Form
function G_sls(t::T, params::Vector{T}) where T<:Real
    η, kᵦ, kᵧ = params

    G = kᵧ + kᵦ*exp(-t*kᵦ/η)
end
G_sls(t::Vector{T}, params::Vector{T}) where T<:Real = G_sls.(t, (params,))

function J_sls(t::T, params::Vector{T}) where T<:Real
    η, kᵦ, kᵧ = params

    c₀ = 1/kᵧ
    c₁ = kᵦ/(kᵧ*(kᵧ + kᵦ))
    τᵣ = η*(kᵧ + kᵦ)/(kᵧ*kᵦ)

    J = c₀ - c₁*exp(-t/τᵣ)
end
J_sls(t::Vector{T}, params::Vector{T}) where T<:Real = J_sls.(t, (params,))

function Gp_sls(ω::T, params::Vector{T}) where T<:Real
    η, kᵦ, kᵧ = params

    τ = η/kᵦ

    denominator = 1 + τ^2*ω^2
    numerator = ω^2*τ^2*kᵦ

    Gp = numerator/denominator + kᵧ
end
Gp_sls(ω::Vector{T}, params::Vector{T}) where T<:Real = Gp_sls.(ω, (params,))

function Gpp_sls(ω::T, params::Vector{T}) where T<:Real
    η, kᵦ, kᵧ = params

    τ = η/kᵦ

    denominator = 1 + τ^2*ω^2
    numerator = ω*τ*kᵦ

    Gpp = numerator/denominator
end
Gpp_sls(ω::Vector{T}, params::Vector{T}) where T<:Real = Gpp_sls.(ω, (params,))

SLS() = RheologyModel(G_sls, J_sls, Gp_sls, Gpp_sls, [1.0, 0.5, 1.0], ["model created with default parameters"])
SLS(params::Vector{T}) where T<:Real = RheologyModel(G_sls, J_sls, Gp_sls, Gpp_sls, params, ["model created by user with parameters $params"])
