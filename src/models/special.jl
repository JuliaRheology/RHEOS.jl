#!/usr/bin/env julia


FractionalSpecial = RheoModelClass(
        # Model name
        name="fractspecial",
        # Model parameters,
        p = [:η, :cᵦ, :β, :k],
        # Relaxation modulus
        G = quote
              G= k + cᵦ*t^(-β)*mittleff(1 - β, 1 - β, -cᵦ*(t^(1 - β))/η)
              return convert(Float64,G)
            end,
        # Creep modulus
        J = quote
              a = η/cᵦ
              b = k*η/cᵦ
              Ĵ(s) = (1/s^2)*(1+a*s^(1-β))/(η+k/s+b*s^(-β))
              return invLaplace(s -> Ĵ(s), t)
            end,
        # Storage modulus
        Gp = quote
               denominator = (η*ω)^2 + (cᵦ*ω^β)^2
               numerator = ((η*ω)^2)*(cᵦ*ω^β)*cos(β*π/2)
               return numerator/denominator + k
             end,
       # Loss modulus
        Gpp = quote
                denominator = (η*ω)^2 + (cᵦ*ω^β)^2
                numerator = ((cᵦ*ω^β)^2)*(η*ω) + ((η*ω)^2)*(cᵦ*ω^β)*sin(β*π/2)
                return numerator/denominator
              end,
        # Network
        info= "
                      ___
                  _____| |__________╱╲__________
                 |    _|_|          ╲╱          |
             ___ |      η              cᵦ, β    |___
                 |                              |
                 |__________╱╲  ╱╲  ╱╲  ________|
                              ╲╱  ╲╱  ╲╱
                                k
               "
        )

#
#
#
#
#
#
# # Fractional Special Model
# function G_fractspecial(t::T, params::Vector{T}) where T<:Real
#     η, cᵦ, β, k = params
#
#     G = k + cᵦ*t^(-β)*mittleff(1 - β, 1 - β, -cᵦ*(t^(1 - β))/η)
#     convert(typeof(params[1]),G)
# end
# G_fractspecial(t::Vector{T}, params::Vector{T}) where T<:Real = G_fractspecial.(t, (params,))
#
# function J_fractspecial(t::Vector{T}, params::Vector{T}) where T<:Real
#     η, cᵦ, β, k = params
#
#     a = η/cᵦ
#     b = k*η/cᵦ
#
#     Ĵ(s) = (1/s^2)*(1+a*s^(1-β))/(η+k/s+b*s^(-β))
#
#     # method gives NaN at time 0.0
#     # J = InverseLaplace.ILt(s -> Jbar(s))
#     # return [J(t_val) for t_val in t]
#
#     return InverseLaplace.talbotarr(s -> Ĵ(s), t)
# end
#
# function Gp_fractspecial(ω::T, params::Vector{T}) where T<:Real
#     η, cᵦ, β, k = params
#
#     denominator = (η*ω)^2 + (cᵦ*ω^β)^2
#     numerator = ((η*ω)^2)*(cᵦ*ω^β)*cos(β*π/2)
#
#     Gp = numerator/denominator + k
# end
# Gp_fractspecial(ω::Vector{T}, params::Vector{T}) where T<:Real = Gp_fractspecial.(ω, (params,))
#
# function Gpp_fractspecial(ω::T, params::Vector{T}) where T<:Real
#     η, cᵦ, β, k = params
#
#     denominator = (η*ω)^2 + (cᵦ*ω^β)^2
#     numerator = ((cᵦ*ω^β)^2)*(η*ω) + ((η*ω)^2)*(cᵦ*ω^β)*sin(β*π/2)
#
#     Gpp = numerator/denominator
# end
# Gpp_fractspecial(ω::Vector{T}, params::Vector{T}) where T<:Real = Gpp_fractspecial.(ω, (params,))
#
# """
#     FractionalSpecial([params::Vector{T}]) where T<:Real
#
# A springpot and dashpot in series, both in parallel with a spring.
# """
# FractionalSpecial(params::Vector{T}) where T<:Real = RheologyModel(G_fractspecial, J_fractspecial, Gp_fractspecial, Gpp_fractspecial, params, ["model created by user with parameters $params"])
# FractionalSpecial() = RheologyModel(G_fractspecial, J_fractspecial, Gp_fractspecial, Gpp_fractspecial, [1.0, 1.0, 0.2, 1.0], ["model created with default parameters"])
