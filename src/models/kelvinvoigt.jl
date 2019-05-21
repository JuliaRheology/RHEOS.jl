#!/usr/bin/env julia

FractionalKelvinVoigt =  RheoModelClass(
        # Model name
        name="fractKV",
        # Model parameters,
        p = [:cₐ, :a, :cᵦ, :β],
        # Relaxation modulus
        G = quote
              cₐ*t^(-a)/gamma(1 - a) + cᵦ*t^(-β)/gamma(1 - β)
            end,
        # Creep modulus
        J = quote
              (t^(a)/cₐ)*mittleff(a - β, 1 + a, -cᵦ*t^(a - β)/cₐ)
            end,
        # Storage modulus
        Gp = quote
                # cosine floating point error work-around
                if a!=1.0 && β!=1.0
                    cₐ*ω^a*cos(a*π/2) + cᵦ*ω^β*cos(β*π/2)
                elseif a==1.0 && β!=1.0
                    cᵦ*ω^β*cos(β*π/2)
                elseif a==1.0 && β==1.0
                    0.0
                end
             end,
        # Loss modulus
        Gpp = quote
                cₐ*ω^a*sin(a*π/2) + cᵦ*ω^β*sin(β*π/2)
              end,
        # Constraints
        constraint = quote
                 all([  (a<1) & (a>0)
                         (β<1) & (β>0)
                          -a+β < 0])
                end,
        # Network
        info= "
                ________ ╱╲ ________
               |         ╲╱  cₐ, a  |
           ____|                    |____
               |                    |
               |________ ╱╲ ________|
                return         ╲╱  cᵦ, β
                "
        )


FractionalKVspring =  RheoModelClass(
        # Model name
        name="fractKVspring",
        # Model parameters,
        p = [:cₐ, :a, :k],
        # Relaxation modulus
        G = quote
              cₐ*t^(-a)/gamma(1 - a) + k
            end,
        # Creep modulus
        J = quote
              (t^(a)/cₐ)*mittleff(a, 1 + a, -k*t^a/cₐ)
            end,
        # Storage modulus
        Gp = quote
                if a!=1.0
                    cₐ*ω^a*cos(a*π/2) + k
                else
                    k
                end
             end,
        # Loss modulus
        Gpp = quote
                cₐ*ω^a*sin(a*π/2)
              end,
        # Constraints
        constraint = quote
                 (a<1) & (a>0)
                end,
        # Network
        info= "
                ________ ╱╲ ________
               |         ╲╱  cₐ, a  |
           ____|                    |____
               |                    |
               |____╱╲  ╱╲  ╱╲  ____|
                      ╲╱  ╲╱  ╲╱
                                k
                "
        )

FractionalKVdashpot =  RheoModelClass(
        # Model name
        name="fractKVdashpot",
        # Model parameters,
        p = [:η, :cᵦ, :β],
        # Relaxation modulus
        G = quote
              t!=0.0 ? cᵦ*t^(-β)/gamma(1 - β) : Inf
            end,
        # Creep modulus
        J = quote
              (t/η)*mittleff(1 - β, 1 + 1, -cᵦ*t^(1.0 - β)/η)
            end,
        # Storage modulus
        Gp = quote
                # cosine floating point error work-around
                if β==1.0
                    0.0
                else
                    cᵦ*ω^β*cos(β*π/2)
                end
             end,
        # Loss modulus
        Gpp = quote
                η*ω + cᵦ*ω^β*sin(β*π/2)
              end,
        # Constraints
        constraint = quote
                 (β<1) & (β>0)
                end,
        # Network
        info= "
                        ___
                _________| |________
               |        _|_| η      |
           ____|                    |____
               |                    |
               |________ ╱╲ ________|
                         ╲╱  cᵦ, β
                "
                )

KelvinVoigt =  RheoModelClass(
        # Model name
        name="KV",
        # Model parameters,
        p = [:η, :k],
        # Relaxation modulus
        G = quote
              t!=0.0 ? k : Inf
            end,
        # Creep modulus
        J = quote
              (1 - exp(-k*t/η))/k
            end,
        # Storage modulus
        Gp = quote
                k
             end,
        # Loss modulus
        Gpp = quote
                η*ω
              end,
        # Network
        info= "
                        ___
                _________| |________
               |        _|_| η      |
           ____|                    |____
               |                    |
               |____╱╲  ╱╲  ╱╲  ____|
                      ╲╱  ╲╱  ╲╱
                                k
                "
                )

# # Fraction Kelvin-Voigt model
# function G_fractKV(t::T, params::Vector{T}) where T<:Real
#     cₐ, a, cᵦ, β = params
#
#     G = cₐ*t^(-a)/gamma(1 - a) + cᵦ*t^(-β)/gamma(1 - β)
# end
# G_fractKV(t::Vector{T}, params::Vector{T}) where T<:Real = G_fractKV.(t, (params,))
#
# function J_fractKV(t::T, params::Vector{T}) where T<:Real
#     cₐ, a, cᵦ, β = params
#
#     J = (t^(a)/cₐ)*mittleff(a - β, 1 + a, -cᵦ*t^(a - β)/cₐ)
# end
# J_fractKV(t::Vector{T}, params::Vector{T}) where T<:Real = J_fractKV.(t, (params,))
#
# function Gp_fractKV(ω::T, params::Vector{T}) where T<:Real
#     cₐ, a, cᵦ, β = params
#
#     # cosine floating point error work-around
#     if a!=1.0 && β!=1.0
#         return cₐ*ω^a*cos(a*π/2) + cᵦ*ω^β*cos(β*π/2)
#     elseif a==1.0 && β!=1.0
#         return cᵦ*ω^β*cos(β*π/2)
#     elseif a==1.0 && β==1.0
#         return 0.0
#     end
# end
# Gp_fractKV(ω::Vector{T}, params::Vector{T}) where T<:Real = Gp_fractKV.(ω, (params,))
#
# function Gpp_fractKV(ω::T, params::Vector{T}) where T<:Real
#     cₐ, a, cᵦ, β = params
#
#     Gpp = cₐ*ω^a*sin(a*π/2) + cᵦ*ω^β*sin(β*π/2)
# end
# Gpp_fractKV(ω::Vector{T}, params::Vector{T}) where T<:Real = Gpp_fractKV.(ω, (params,))
#
# """
#     FractionalKelvinVoigt([params::Vector{T}]) where T<:Real
#
# Two springpots in parallel.
# """
# FractionalKelvinVoigt(params::Vector{T}) where T<:Real = RheologyModel(G_fractKV, J_fractKV, Gp_fractKV, Gpp_fractKV, params, ["model created by user with parameters $params"])
# FractionalKelvinVoigt() = RheologyModel(G_fractKV, J_fractKV, Gp_fractKV, Gpp_fractKV, [2.0, 0.2, 1.0, 0.5], ["model created with default parameters"])

# # Fraction Kelvin-Voigt model with lower spring-pot specialized to a spring
# function G_fractKVspring(t::T, params::Vector{T}) where T<:Real
#     cₐ, a, k = params
#
#     G = cₐ*t^(-a)/gamma(1 - a) + k
# end
# G_fractKVspring(t::Vector{T}, params::Vector{T}) where T<:Real = G_fractKVspring.(t, (params,))
#
# function J_fractKVspring(t::T, params::Vector{T}) where T<:Real
#     cₐ, a, k = params
#
#     J = (t^(a)/cₐ)*mittleff(a, 1 + a, -k*t^a/cₐ)
# end
# J_fractKVspring(t::Vector{T}, params::Vector{T}) where T<:Real = J_fractKVspring.(t, (params,))
#
# function Gp_fractKVspring(ω::T, params::Vector{T}) where T<:Real
#     cₐ, a, k = params
#
#     # cosine floating point error work-around
#     if a!=1.0
#         return cₐ*ω^a*cos(a*π/2) + k
#     else
#         return k
#     end
# end
# Gp_fractKVspring(ω::Vector{T}, params::Vector{T}) where T<:Real = Gp_fractKVspring.(ω, (params,))
#
# function Gpp_fractKVspring(ω::T, params::Vector{T}) where T<:Real
#     cₐ, a, k = params
#
#     Gpp = cₐ*ω^a*sin(a*π/2)
# end
# Gpp_fractKVspring(ω::Vector{T}, params::Vector{T}) where T<:Real = Gpp_fractKVspring.(ω, (params,))
#
# """
#     FractionalKVspring([params::Vector{T}]) where T<:Real
#
# A springpot and spring in parallel.
# """
# FractionalKVspring(params::Vector{T}) where T<:Real = RheologyModel(G_fractKVspring, J_fractKVspring, Gp_fractKVspring, Gpp_fractKVspring, params, ["model created by user with parameters $params"])
# FractionalKVspring() = RheologyModel(G_fractKVspring, J_fractKVspring, Gp_fractKVspring, Gpp_fractKVspring, [2.0, 0.2, 1.0], ["model created with default parameters"])

# # Fraction Kelvin-Voigt model with upper spring-pot specialized to a dash-pot
# function G_fractKVdashpot(t::T, params::Vector{T}) where T<:Real
#     η, cᵦ, β = params
#
#     # G = η*t.^(-1)/gamma(1.0 - 1.0) + cᵦ*t.^(-β)/gamma(1 - β)
#
#     # dashpot (first term) denominator is undefined as 0*∞, physical intuition
#     # suggest an instantaneous singularity as an dashpot/spring model. As this
#     # is already captured by the second power law term, just keep that one in.
#     G = t!=0.0 ? cᵦ*t^(-β)/gamma(1 - β) : Inf
# end
# G_fractKVdashpot(t::Vector{T}, params::Vector{T}) where T<:Real = G_fractKVdashpot.(t, (params,))
#
# function J_fractKVdashpot(t::T, params::Vector{T}) where T<:Real
#     η, cᵦ, β = params
#
#     J = (t/η)*mittleff(1 - β, 1 + 1, -cᵦ*t^(1.0 - β)/η)
# end
# J_fractKVdashpot(t::Vector{T}, params::Vector{T}) where T<:Real = J_fractKVdashpot.(t, (params,))
#
# function Gp_fractKVdashpot(ω::T, params::Vector{T}) where T<:Real
#     η, cᵦ, β = params
#
#     # cosine floating point error work-around
#     if β==1.0
#         return 0.0
#     else
#         return cᵦ*ω^β*cos(β*π/2)
#     end
# end
# Gp_fractKVdashpot(ω::Vector{T}, params::Vector{T}) where T<:Real = Gp_fractKVdashpot.(ω, (params,))
#
# function Gpp_fractKVdashpot(ω::T, params::Vector{T}) where T<:Real
#     η, cᵦ, β = params
#
#     Gpp = η*ω + cᵦ*ω^β*sin(β*π/2)
# end
# Gpp_fractKVdashpot(ω::Vector{T}, params::Vector{T}) where T<:Real = Gpp_fractKVdashpot.(ω, (params,))
#
# """
#     FractionalKVdashpot([params::Vector{T}]) where T<:Real
#
# A springpot and dashpot in parallel.
# """
# FractionalKVdashpot(params::Vector{T}) where T<:Real = RheologyModel(G_fractKVdashpot, J_fractKVdashpot, Gp_fractKVdashpot, Gpp_fractKVdashpot, params, ["model created by user with parameters $params"])
# FractionalKVdashpot() = RheologyModel(G_fractKVdashpot, J_fractKVdashpot, Gp_fractKVdashpot, Gpp_fractKVdashpot, [2.0, 0.2, 1.0], ["model created with default parameters"])
#
# # Standard Kelvin-Voigt model
# function G_KV(t::T, params::Vector{T}) where T<:Real
#     η, k = params
#
#     # G = k # + η*δ(t) in Findley et al. (1989)
#     G = t!=0.0 ? k : Inf
# end
# G_KV(t::Vector{T}, params::Vector{T}) where T<:Real = G_KV.(t, (params,))
#
# function J_KV(t::T, params::Vector{T}) where T<:Real
#     η, k = params
#
#     J = (1 - exp(-k*t/η))/k
# end
# J_KV(t::Vector{T}, params::Vector{T}) where T<:Real = J_KV.(t, (params,))
#
# function Gp_KV(ω::T, params::Vector{T}) where T<:Real
#     η, k = params
#
#     Gp = k
# end
# Gp_KV(ω::Vector{T}, params::Vector{T}) where T<:Real = Gp_KV.(ω, (params,))
#
# function Gpp_KV(ω::T, params::Vector{T}) where T<:Real
#     η, k = params
#
#     Gpp = η*ω
# end
# Gpp_KV(ω::Vector{T}, params::Vector{T}) where T<:Real = Gpp_KV.(ω, (params,))
#
# """
#     KelvinVoigt([params::Vector{T}]) where T<:Real
#
# A spring and dashpot in parallel.
# """
# KelvinVoigt(params::Vector{T}) where T<:Real = RheologyModel(G_KV, J_KV, Gp_KV, Gpp_KV, params, ["model created by user with parameters $params"])
# KelvinVoigt() = RheologyModel(G_KV, J_KV, Gp_KV, Gpp_KV, [2.0, 1.0], ["model created with default parameters"])
