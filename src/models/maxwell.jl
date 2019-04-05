#!/usr/bin/env julia

FractionalMaxwell = RheoModelClass(
          # Model name
          name="fractmaxwell",
          # Model parameters,
          p = [:cₐ, :a, :cᵦ, :β],
          # Relaxation modulus
          G = quote
                 return cᵦ*t^(-β)*mittleff(a - β, 1 - β, -cᵦ*t^(a - β)/cₐ)
              end,
          # Creep modulus
          J = quote
                return t^(a)/(cₐ*gamma(1 + a)) + t^(β)/(cᵦ*gamma(1 + β))
              end,
          # Storage modulus
          Gp = quote
                 denominator = (cₐ*ω^a)^2 + (cᵦ*ω^β)^2 + 2*(cₐ*ω^a)*(cᵦ*ω^β)*cos((a-β)*π/2)
                 numerator = ((cᵦ*ω^β)^2)*(cₐ*ω^a)*cos(a*π/2) + ((cₐ*ω^a)^2)*(cᵦ*ω^β)*cos(β*π/2)
                 return numerator/denominator
               end,
         # Loss modulus
          Gpp = quote
                  denominator = (cₐ*ω^a)^2 + (cᵦ*ω^β)^2 + 2*(cₐ*ω^a)*(cᵦ*ω^β)*cos((a-β)*π/2)
                  numerator = ((cᵦ*ω^β)^2)*(cₐ*ω^a)*sin(a*π/2) + ((cₐ*ω^a)^2)*(cᵦ*ω^β)*sin(β*π/2)
                  return numerator/denominator
                end,
          # Constraints
          Ineq = quote
                   return [(a<1) & (a>0)
                           (β<1) & (β>0)
                            -a+β < 0]
                  end,
          # Network
          info= "
             ___╱╲__________╱╲____
                ╲╱          ╲╱
                  cₐ,a         cᵦ, β
                 "
          )

FractionalMaxwellSpring = RheoModelClass(
        # Model name
        name="fractmaxwell_spring",
        # Model parameters,
        p = [:cₐ, :a, :k],
        # Relaxation modulus
        G = quote
                # special case when α==0.0, 2 springs in series
                a==0.0 && return 1.0/(1.0/cₐ + 1.0/k)
                # normal case
                return k*mittleff(a, -k*t^a/cₐ)
            end,
        # Creep modulus
        J = quote
              return t^(a)/(cₐ*gamma(1 + a)) + 1/k
            end,
        # Storage modulus
        Gp = quote
                denominator = (cₐ*ω^a)^2 + k^2 + 2*(cₐ*ω^a)*k*cos(a*π/2)
                numerator = k^2*(cₐ*ω^a)*cos(a*π/2) + (cₐ*ω^a)^2*k
                return numerator/denominator
             end,
       # Loss modulus
        Gpp = quote
                denominator = (cₐ*ω^a)^2 + k^2 + 2*(cₐ*ω^a)*k*cos(a*π/2)
                numerator = k^2*(cₐ*ω^a)*sin(a*π/2)
                return numerator/denominator
              end,
        # Constraints
        Ineq = quote
                 return (a<1) & (a>0)
                end,
        # Network
        info= "
           ___╱╲_________╱╲  ╱╲  ╱╲  ________
              ╲╱           ╲╱  ╲╱  ╲╱
                cₐ,a               k
               "
        )


FractionalMaxwellDashpot = RheoModelClass(
          # Model name
          name="fractmaxwell_dashpot",
          # Model parameters,
          p = [:η, :cᵦ, :β],
          # Relaxation modulus
          G = quote
                return cᵦ*t^(-β)*mittleff(1 - β, 1 - β, -cᵦ*t^(1 - β)/η)
              end,
          # Creep modulus
          J = quote
                return t/η + t^β/(cᵦ*gamma(1 + β))
              end,
          # Storage modulus
          Gp = quote
                  denominator = (η*ω)^2 + (cᵦ*ω^β)^2 + 2*(η*ω)*(cᵦ*ω^β)*cos((1-β)*π/2)
                  numerator = ((η*ω)^2)*(cᵦ*ω^β)*cos(β*π/2)
                  return numerator/denominator
               end,
         # Loss modulus
          Gpp = quote
                  denominator = (η*ω)^2 + (cᵦ*ω^β)^2 + 2*(η*ω)*(cᵦ*ω^β)*cos((1-β)*π/2)
                  numerator = ((cᵦ*ω^β)^2)*(η*ω) + ((η*ω)^2)*(cᵦ*ω^β)*sin(β*π/2)
                  return numerator/denominator
                end,
          # Constraints
          Ineq = quote
                   return (β<1) & (β>0)
                  end,
          # Network
          info= "
                  ___
              _____| |_________╱╲____
                  _|_|         ╲╱
                    η            cᵦ, β
                 "
          )


Maxwell = RheoModelClass(
        # Model name
        name="maxwell",
        # Model parameters,
        p = [:η, :k],
        # Relaxation modulus
        G = quote
              return k*exp(-k*t/η)
            end,
        # Creep modulus
        J = quote
              return t/η + 1/k
            end,
        # Storage modulus
        Gp = quote
                denominator = 1 + η^2*ω^2/k^2
                numerator = η^2*ω^2/k
                return numerator/denominator
             end,
       # Loss modulus
        Gpp = quote
                denominator = 1 + η^2*ω^2/k^2
                numerator = η*ω
                return numerator/denominator
              end,
        # Network
        info= "
                ___
            _____| |________╱╲  ╱╲  ╱╲  ___
                _|_|          ╲╱  ╲╱  ╲╱
                  η                  k
               "
        )


# # Fractional Maxwell Model
# function G_fractmaxwell(t::T, params::Vector{T}) where T<:Real
#     cₐ, a, cᵦ, β = params
#
#     G = cᵦ*t^(-β)*mittleff(a - β, 1 - β, -cᵦ*t^(a - β)/cₐ)
# end
# G_fractmaxwell(t::Vector{T}, params::Vector{T}) where T<:Real = G_fractmaxwell.(t, (params,))
#
# function J_fractmaxwell(t::T, params::Vector{T}) where T<:Real
#     cₐ, a, cᵦ, β = params
#
#     J = t^(a)/(cₐ*gamma(1 + a)) + t^(β)/(cᵦ*gamma(1 + β))
# end
# J_fractmaxwell(t::Vector{T}, params::Vector{T}) where T<:Real = J_fractmaxwell.(t, (params,))
#
# function Gp_fractmaxwell(ω::T, params::Vector{T}) where T<:Real
#     cₐ, a, cᵦ, β = params
#
#     denominator = (cₐ*ω^a)^2 + (cᵦ*ω^β)^2 + 2*(cₐ*ω^a)*(cᵦ*ω^β)*cos((a-β)*π/2)
#
#     numerator = ((cᵦ*ω^β)^2)*(cₐ*ω^a)*cos(a*π/2) + ((cₐ*ω^a)^2)*(cᵦ*ω^β)*cos(β*π/2)
#
#     Gp = numerator/denominator
# end
# Gp_fractmaxwell(ω::Vector{T}, params::Vector{T}) where T<:Real = Gp_fractmaxwell.(ω, (params,))
#
# function Gpp_fractmaxwell(ω::T, params::Vector{T}) where T<:Real
#     cₐ, a, cᵦ, β = params
#
#     denominator = (cₐ*ω^a)^2 + (cᵦ*ω^β)^2 + 2*(cₐ*ω^a)*(cᵦ*ω^β)*cos((a-β)*π/2)
#
#     numerator = ((cᵦ*ω^β)^2)*(cₐ*ω^a)*sin(a*π/2) + ((cₐ*ω^a)^2)*(cᵦ*ω^β)*sin(β*π/2)
#
#     Gpp = numerator/denominator
# end
# Gpp_fractmaxwell(ω::Vector{T}, params::Vector{T}) where T<:Real = Gpp_fractmaxwell.(ω, (params,))
#
# """
#     FractionalMaxwell([params::Vector{T}]) where T<:Real
#
# Two springpots in series.
# """
# FractionalMaxwell(params::Vector{T}) where T<:Real = RheologyModel(G_fractmaxwell, J_fractmaxwell, Gp_fractmaxwell, Gpp_fractmaxwell, params, ["model created by user with parameters $params"])
# FractionalMaxwell() = RheologyModel(G_fractmaxwell, J_fractmaxwell, Gp_fractmaxwell, Gpp_fractmaxwell, [2.0, 0.2, 1.0, 0.5], ["model created with default parameters"])

# # Fractional Maxwell Model (β=0, second spring-pot specialized to spring)
# function G_fractmaxwell_spring(t::T, params::Vector{T}) where T<:Real
#     cₐ, a, k = params
#
#     # special case when α==0.0, 2 springs in series
#     a==0.0 && return 1.0/(1.0/cₐ + 1.0/k)
#
#     # normal case
#     G = k*mittleff(a, -k*t^a/cₐ)
#     return G
# end
# G_fractmaxwell_spring(t::Vector{T}, params::Vector{T}) where T<:Real = G_fractmaxwell_spring.(t, (params,))
#
# function J_fractmaxwell_spring(t::T, params::Vector{T}) where T<:Real
#     cₐ, a, k = params
#
#     J = t^(a)/(cₐ*gamma(1 + a)) + 1/k
# end
# J_fractmaxwell_spring(t::Vector{T}, params::Vector{T}) where T<:Real = J_fractmaxwell_spring.(t, (params,))
#
# function Gp_fractmaxwell_spring(ω::T, params::Vector{T}) where T<:Real
#     cₐ, a, k = params
#
#     denominator = (cₐ*ω^a)^2 + k^2 + 2*(cₐ*ω^a)*k*cos(a*π/2)
#
#     numerator = k^2*(cₐ*ω^a)*cos(a*π/2) + (cₐ*ω^a)^2*k
#
#     Gp = numerator/denominator
# end
# Gp_fractmaxwell_spring(ω::Vector{T}, params::Vector{T}) where T<:Real = Gp_fractmaxwell_spring.(ω, (params,))
#
# function Gpp_fractmaxwell_spring(ω::T, params::Vector{T}) where T<:Real
#     cₐ, a, k = params
#
#     denominator = (cₐ*ω^a)^2 + k^2 + 2*(cₐ*ω^a)*k*cos(a*π/2)
#
#     numerator = k^2*(cₐ*ω^a)*sin(a*π/2)
#
#     Gpp = numerator/denominator
# end
# Gpp_fractmaxwell_spring(ω::Vector{T}, params::Vector{T}) where T<:Real = Gpp_fractmaxwell_spring.(ω, (params,))
#
# """
#     FractionalMaxwellSpring([params::Vector{T}]) where T<:Real
#
# A springpot and spring in series.
# """
# FractionalMaxwellSpring(params::Vector{T}) where T<:Real = RheologyModel(G_fractmaxwell_spring, J_fractmaxwell_spring, Gp_fractmaxwell_spring, Gpp_fractmaxwell_spring, params, ["model created by user with parameters $params"])
# FractionalMaxwellSpring() = RheologyModel(G_fractmaxwell_spring, J_fractmaxwell_spring, Gp_fractmaxwell_spring, Gpp_fractmaxwell_spring, [2.0, 0.2, 1.0], ["model created with default parameters"])
#
# # Fractional Maxwell Model (α=0, first spring-pot specialized to dash-pot)
# function G_fractmaxwell_dashpot(t::T, params::Vector{T}) where T<:Real
#     η, cᵦ, β = params
#
#     G = cᵦ*t^(-β)*mittleff(1 - β, 1 - β, -cᵦ*t^(1 - β)/η)
# end
# G_fractmaxwell_dashpot(t::Vector{T}, params::Vector{T}) where T<:Real = G_fractmaxwell_dashpot.(t, (params,))
#
# function J_fractmaxwell_dashpot(t::T, params::Vector{T}) where T<:Real
#     η, cᵦ, β = params
#
#     J = t/η + t^β/(cᵦ*gamma(1 + β))
# end
# J_fractmaxwell_dashpot(t::Vector{T}, params::Vector{T}) where T<:Real = J_fractmaxwell_dashpot.(t, (params,))
#
# function Gp_fractmaxwell_dashpot(ω::T, params::Vector{T}) where T<:Real
#     η, cᵦ, β = params
#
#     denominator = (η*ω)^2 + (cᵦ*ω^β)^2 + 2*(η*ω)*(cᵦ*ω^β)*cos((1-β)*π/2)
#
#     numerator = ((η*ω)^2)*(cᵦ*ω^β)*cos(β*π/2)
#
#     Gp = numerator/denominator
# end
# Gp_fractmaxwell_dashpot(ω::Vector{T}, params::Vector{T}) where T<:Real = Gp_fractmaxwell_dashpot.(ω, (params,))
#
# function Gpp_fractmaxwell_dashpot(ω::T, params::Vector{T}) where T<:Real
#     η, cᵦ, β = params
#
#     denominator = (η*ω)^2 + (cᵦ*ω^β)^2 + 2*(η*ω)*(cᵦ*ω^β)*cos((1-β)*π/2)
#
#     numerator = ((cᵦ*ω^β)^2)*(η*ω) + ((η*ω)^2)*(cᵦ*ω^β)*sin(β*π/2)
#
#     Gpp = numerator/denominator
# end
# Gpp_fractmaxwell_dashpot(ω::Vector{T}, params::Vector{T}) where T<:Real = Gpp_fractmaxwell_dashpot.(ω, (params,))
#
# """
#     FractionalMaxwellDashpot([params::Vector{T}]) where T<:Real
#
# A springpot and dashpot in series.
# """
# FractionalMaxwellDashpot(params::Vector{T}) where T<:Real = RheologyModel(G_fractmaxwell_dashpot, J_fractmaxwell_dashpot, Gp_fractmaxwell_dashpot, Gpp_fractmaxwell_dashpot, params, ["model created by user with parameters $params"])
# FractionalMaxwellDashpot() = RheologyModel(G_fractmaxwell_dashpot, J_fractmaxwell_dashpot, Gp_fractmaxwell_dashpot, Gpp_fractmaxwell_dashpot, [2.0, 1.0, 0.5], ["model created with default parameters"])
#
# # Maxwell Model (From Findley, Lai, Onaran for comparison/debug)
# function G_maxwell(t::T, params::Vector{T}) where T<:Real
#     η, k = params
#
#     G = k*exp(-k*t/η)
# end
# G_maxwell(t::Vector{T}, params::Vector{T}) where T<:Real = G_maxwell.(t, (params,))
#
# function J_maxwell(t::T, params::Vector{T}) where T<:Real
#     η, k = params
#
#     J = t/η + 1/k
# end
# J_maxwell(t::Vector{T}, params::Vector{T}) where T<:Real = J_maxwell.(t, (params,))
#
# function Gp_maxwell(ω::T, params::Vector{T}) where T<:Real
#     η, k = params
#
#     denominator = 1 + η^2*ω^2/k^2
#
#     numerator = η^2*ω^2/k
#
#     Gp = numerator/denominator
# end
# Gp_maxwell(ω::Vector{T}, params::Vector{T}) where T<:Real = Gpp_fractmaxwell_dashpot.(ω, (params,))
#
# function Gpp_maxwell(ω::T, params::Vector{T}) where T<:Real
#     η, k = params
#
#     denominator = 1 + η^2*ω^2/k^2
#
#     numerator = η*ω
#
#     Gpp = numerator/denominator
# end
# Gpp_maxwell(ω::Vector{T}, params::Vector{T}) where T<:Real = Gpp_maxwell.(ω, (params,))
#
# """
#     Maxwell([params::Vector{T}]) where T<:Real
#
# A spring and dashpot in series.
# """
# Maxwell(params::Vector{T}) where T<:Real = RheologyModel(G_maxwell, J_maxwell, Gp_maxwell, Gpp_maxwell, params, ["model created by user with parameters $params"])
# Maxwell() = RheologyModel(G_maxwell, J_maxwell, Gp_maxwell, Gpp_maxwell, [2.0, 1.0], ["model created with default parameters"])
