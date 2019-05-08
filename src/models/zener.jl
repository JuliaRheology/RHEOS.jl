#!/usr/bin/env julia

FractionalZener = RheoModelClass(
          # Model name
          name="fraczener",
          # Model parameters,
          p = [ :cₐ, :a, :cᵦ, :β, :cᵧ, :γ],
          # Relaxation modulus
          G = quote
                 return cᵦ*t^(-β)*mittleff(a - β, 1 - β, -cᵦ*t^(a - β)/cₐ) + cᵧ*t^(-γ)/gamma(1 - γ)
              end,
          # Creep modulus
          J = quote
                  Ĵ(s) = (1/s)*(cₐ*s^a + cᵦ*s^β)/(cₐ*s^a*cᵦ*s^β + cᵧ*s^γ*(cₐ*s^a + cᵦ*s^β))
                  return InverseLaplace.talbotarr(s -> Ĵ(s), t)
              end,
          # Storage modulus
          Gp = quote
                   denominator = (cₐ*ω^a)^2 + (cᵦ*ω^β)^2 + 2*(cₐ*ω^a)*(cᵦ*ω^β)*cos((a - β)*π/2)
                   numerator = ((cᵦ*ω^β)^2)*(cₐ*ω^a)*cos(a*π/2) + ((cₐ*ω^a)^2)*(cᵦ*ω^β)*cos(β*π/2)
                   return numerator/denominator + cᵧ*ω^γ*cos(γ*π/2)
               end,
         # Loss modulus
          Gpp = quote
                   denominator = (cₐ*ω^a)^2 + (cᵦ*ω^β)^2 + 2*(cₐ*ω^a)*(cᵦ*ω^β)*cos((a - β)*π/2)
                   numerator = ((cᵦ*ω^β)^2)*(cₐ*ω^a)*sin(a*π/2) + ((cₐ*ω^a)^2)*(cᵦ*ω^β)*sin(β*π/2)
                   return numerator/denominator + cᵧ*ω^γ*sin(γ*π/2)
                end,
          # Constraints
          Ineq = quote
                   return [(a<1) & (a>0)
                           (β<1) & (β>0)
                            -a+β < 0]
                  end,
          # Network
          info= "

                  ______╱╲__________╱╲______
                 |      ╲╱          ╲╱      |
          _______|      cₐ,a         cᵦ, β  |_______
                 |                          |
                 |____________╱╲____________|
                              ╲╱
                              cᵧ, γ
                     "
          )


FractionalSLS = RheoModelClass(
        # Model name
        name="fracsls",
        # Model parameters,
        p = [ :cₐ, :a, :kᵦ, :kᵧ],
        # Relaxation modulus
        G = quote
               return kᵦ*mittleff(a, -(kᵦ/cₐ)*t^a) + kᵧ
            end,
        # Creep modulus
        J = quote
                Ĵ(s) = (1/s)*(cₐ*s^a + kᵦ)/(cₐ*s^a*kᵦ + kᵧ*(cₐ*s^a + kᵦ))
                return InverseLaplace.talbotarr(s -> Ĵ(s), t)
            end,
        # Storage modulus
        Gp = quote
                denominator = (cₐ*ω^a)^2 + kᵦ^2 + 2*(cₐ*ω^a)*kᵦ*cos(a*π/2)
                numerator = kᵦ^2*(cₐ*ω^a)*cos(a*π/2) + ((cₐ*ω^a)^2)*kᵦ
                return numerator/denominator + kᵧ
             end,
       # Loss modulus
        Gpp = quote
                denominator = (cₐ*ω^a)^2 + kᵦ^2 + 2*(cₐ*ω^a)*kᵦ*cos(a*π/2)
                numerator = kᵦ^2*(cₐ*ω^a)*sin(a*π/2)
                return numerator/denominator
              end,
        # Constraints
        Ineq = quote
                 return (a<1) & (a>0)
                end,
        # Network
        info= "

                _____╱╲_______╱╲  ╱╲  ╱╲  _____
               |     ╲╱         ╲╱  ╲╱  ╲╱     |
        _______|       cₐ,a              kᵦ    |_______
               |                               |
               |__________╱╲  ╱╲  ╱╲  _________|
                            ╲╱  ╲╱  ╲╱
                                 kᵧ
                   "
        )


SLS = RheoModelClass(
          # Model name
          name="SLS",
          # Model parameters,
          p = [ :η, :kᵦ, :kᵧ],
          # Relaxation modulus
          G = quote
                 return kᵧ + kᵦ*exp(-t*kᵦ/η)
              end,
          # Creep modulus
          J = quote
                  c₀ = 1/kᵧ
                  c₁ = kᵦ/(kᵧ*(kᵧ + kᵦ))
                  τᵣ = η*(kᵧ + kᵦ)/(kᵧ*kᵦ)
                  return c₀ - c₁*exp(-t/τᵣ)
              end,
          # Storage modulus
          Gp = quote
                  τ = η/kᵦ
                  denominator = 1 + τ^2*ω^2
                  numerator = ω^2*τ^2*kᵦ
                  return numerator/denominator + kᵧ
               end,
         # Loss modulus
          Gpp = quote
                  τ = η/kᵦ
                  denominator = 1 + τ^2*ω^2
                  numerator = ω*τ*kᵦ
                  return numerator/denominator
                end,
          # Network
          info= "
                      ___
                  _____| |________╱╲  ╱╲  ╱╲  ___
                 |    _|_|          ╲╱  ╲╱  ╲╱   |
          _______|      η                  kᵦ    |_______
                 |                               |
                 |__________╱╲  ╱╲  ╱╲  _________|
                              ╲╱  ╲╱  ╲╱
                                   kᵧ
                     "
          )

FractionalJeffreys = RheoModelClass(
                    # Model name
                    name="fjeff",
                    # Model parameters,
                    p = [ :ηₐ, :cᵦ, :β, :ηᵧ],
                    # Relaxation modulus
                    G = quote
                            diracterm = t!=0.0 ? 0.0 : Inf
                            return cᵦ*t^(-β)*mittleff(1-β, 1-β, -cᵦ*t^(1-β)/ηₐ) + ηᵧ*diracterm
                        end,
                    # Creep modulus
                    J = quote
                            Ĵ(s) = (1/s)*(ηₐ*s + cᵦ*s^β)/(ηₐ*s*cᵦ*s^β + ηᵧ*s*(ηₐ*s + cᵦ*s^β))
                            return InverseLaplace.talbotarr(s -> Ĵ(s), t)
                        end,
                    # Storage modulus
                    Gp = quote
                            denominator = (ηₐ*ω)^2 + (cᵦ*ω^β)^2 + 2*(ηₐ*ω)*(cᵦ*ω^β)*cos((1 - β)*π/2)
                            numerator = ((ηₐ*ω)^2)*(cᵦ*ω^β)*cos(β*π/2)
                            # cosine floating point error work-around
                            if β==1.0
                                return 0.0
                            else
                                return numerator/denominator
                            end
                         end,
                   # Loss modulus
                    Gpp = quote
                            denominator = (ηₐ*ω)^2 + (cᵦ*ω^β)^2 + 2*(ηₐ*ω)*(cᵦ*ω^β)*cos((1 - β)*π/2)
                            numerator = ((cᵦ*ω^β)^2)*(ηₐ*ω) + ((ηₐ*ω)^2)*(cᵦ*ω^β)*sin(β*π/2)
                            return numerator/denominator + ηᵧ*ω
                          end,
                    # Constraints
                    Ineq = quote
                             return (β<1) & (β>0)
                            end,
                    # Network
                    info= "

                                    ___
                            _________| |_________╱╲________
                           |        _|_|         ╲╱        |
                    _______|          ηₐ            cᵦ, β  |_______
                           |              ___              |
                           |_______________| |_____________|
                                          _|_|
                                             ηᵧ
                               "
                    )


Jeffreys = RheoModelClass(
                # Model name
                name="jeffreys",
                # Model parameters,
                p = [ :ηₐ, :k, :ηᵧ],
                # Relaxation modulus
                G = quote
                        diracterm = t!=0.0 ? 0.0 : Inf
                        return k*exp(-k*t/ηₐ) + ηᵧ*diracterm
                    end,
                # Creep modulus
                J = quote
                        Ĵ(s) = (1/s)*(ηₐ*s + k)/(ηₐ*s*k + ηᵧ*s*(ηₐ*s + k))
                        return InverseLaplace.talbotarr(s -> Ĵ(s), t)
                    end,
                # Storage modulus
                Gp = quote
                        denominator = (ηₐ*ω)^2 + k^2
                        numerator = ((ηₐ*ω)^2)*k
                        return numerator/denominator
                     end,
               # Loss modulus
                Gpp = quote
                        denominator = (ηₐ*ω)^2 + k^2
                        numerator = (k^2)*(ηₐ*ω)
                        return numerator/denominator + ηᵧ*ω
                      end,
                # Network
                info= "

                              ___
                        _______| |_______╱╲  ╱╲  ╱╲  ___
                       |      _|_|         ╲╱  ╲╱  ╲╱  |
                _______|          ηₐ            k      |_______
                       |              ___              |
                       |_______________| |_____________|
                                      _|_|
                                         ηᵧ
                           "
                )







          # # Fractional Jeffrey's Model in Zener Form
          # function G_fjeff(t::T, params::Vector{T}) where T<:Real
          #     ηₐ, cᵦ, β, ηᵧ = params
          #
          #     diracterm = t!=0.0 ? 0.0 : Inf
          #
          #     G = cᵦ*t^(-β)*mittleff(1-β, 1-β, -cᵦ*t^(1-β)/ηₐ) + ηᵧ*diracterm
          # end
          # G_fjeff(t::Vector{T}, params::Vector{T}) where T<:Real = G_fjeff.(t, (params,))
          #
          # function J_fjeff(t::Vector{T}, params::Vector{T}) where T<:Real
          #     ηₐ, cᵦ, β, ηᵧ = params
          #
          #     Ĵ(s) = (1/s)*(ηₐ*s + cᵦ*s^β)/(ηₐ*s*cᵦ*s^β + ηᵧ*s*(ηₐ*s + cᵦ*s^β))
          #
          #     # # method gives NaN at time 0.0
          #     # J = InverseLaplace.ILt(s -> Ĵ(s))
          #     # return [J(t_val) for t_val in t]
          #
          #     return InverseLaplace.talbotarr(s -> Ĵ(s), t)
          # end
          #
          # function Gp_fjeff(ω::T, params::Vector{T}) where T<:Real
          #     ηₐ, cᵦ, β, ηᵧ = params
          #
          #     denominator = (ηₐ*ω)^2 + (cᵦ*ω^β)^2 + 2*(ηₐ*ω)*(cᵦ*ω^β)*cos((1 - β)*π/2)
          #
          #     numerator = ((ηₐ*ω)^2)*(cᵦ*ω^β)*cos(β*π/2)
          #
          #     # cosine floating point error work-around
          #     if β==1.0
          #         return 0.0
          #     else
          #         return numerator/denominator
          #     end
          # end
          # Gp_fjeff(ω::Vector{T}, params::Vector{T}) where T<:Real = Gp_fjeff.(ω, (params,))
          #
          # function Gpp_fjeff(ω::T, params::Vector{T}) where T<:Real
          #     ηₐ, cᵦ, β, ηᵧ = params
          #
          #     denominator = (ηₐ*ω)^2 + (cᵦ*ω^β)^2 + 2*(ηₐ*ω)*(cᵦ*ω^β)*cos((1 - β)*π/2)
          #
          #     numerator = ((cᵦ*ω^β)^2)*(ηₐ*ω) + ((ηₐ*ω)^2)*(cᵦ*ω^β)*sin(β*π/2)
          #
          #     Gpp = numerator/denominator + ηᵧ*ω
          # end
          # Gpp_fjeff(ω::Vector{T}, params::Vector{T}) where T<:Real = Gpp_fjeff.(ω, (params,))
          #
          # """
          #     FractionalJeffreys(params::Vector{T}) where T<:Real
          # """
          # FractionalJeffreys(params::Vector{T}) where T<:Real = RheologyModel(G_fjeff, J_fjeff, Gp_fjeff, Gpp_fjeff, params, ["model created by user with parameters $params"])
          # FractionalJeffreys() = RheologyModel(G_fjeff, J_fjeff, Gp_fjeff, Gpp_fjeff, [1.0, 1.0, 0.5, 1.0], ["model created with default parameters"])
          #
          # # Jeffrey's Model in Zener Form
          # function G_jeffreys(t::T, params::Vector{T}) where T<:Real
          #     ηₐ, k, ηᵧ = params
          #
          #     diracterm = t!=0.0 ? 0.0 : Inf
          #
          #     G = k*exp(-k*t/ηₐ) + ηᵧ*diracterm
          #
          # end
          # G_jeffreys(t::Vector{T}, params::Vector{T}) where T<:Real = G_jeffreys.(t, (params,))
          #
          # function J_jeffreys(t::Vector{T}, params::Vector{T}) where T<:Real
          #     ηₐ, k, ηᵧ = params
          #
          #     Ĵ(s) = (1/s)*(ηₐ*s + k)/(ηₐ*s*k + ηᵧ*s*(ηₐ*s + k))
          #
          #     # # method gives NaN at time 0.0
          #     # J = InverseLaplace.ILt(s -> Ĵ(s))
          #     # return [J(t_val) for t_val in t]
          #
          #     return InverseLaplace.talbotarr(s -> Ĵ(s), t)
          # end
          #
          # function Gp_jeffreys(ω::T, params::Vector{T}) where T<:Real
          #     ηₐ, k, ηᵧ = params
          #
          #     denominator = (ηₐ*ω)^2 + k^2
          #
          #     numerator = ((ηₐ*ω)^2)*k
          #
          #     Gp = numerator/denominator
          # end
          # Gp_jeffreys(ω::Vector{T}, params::Vector{T}) where T<:Real = Gp_jeffreys.(ω, (params,))
          #
          # function Gpp_jeffreys(ω::T, params::Vector{T}) where T<:Real
          #     ηₐ, k, ηᵧ = params
          #
          #     denominator = (ηₐ*ω)^2 + k^2
          #
          #     numerator = (k^2)*(ηₐ*ω)
          #
          #     Gpp = numerator/denominator + ηᵧ*ω
          # end
          # Gpp_jeffreys(ω::Vector{T}, params::Vector{T}) where T<:Real = Gpp_jeffreys.(ω, (params,))
          #
          # """
          #     Jeffreys(params::Vector{T}) where T<:Real
          # """
          # Jeffreys(params::Vector{T}) where T<:Real = RheologyModel(G_jeffreys, J_jeffreys, Gp_jeffreys, Gpp_jeffreys, params, ["model created by user with parameters $params"])
          # Jeffreys() = RheologyModel(G_jeffreys, J_jeffreys, Gp_jeffreys, Gpp_jeffreys, [1.0, 1.0, 1.0], ["model created with default parameters"])


# # Fractional Zener Model
# function G_fraczener(t::T, params::Vector{T}) where T<:Real
#     cₐ, a, cᵦ, β, cᵧ, γ = params
#
#     G = cᵦ*t^(-β)*mittleff(a - β, 1 - β, -cᵦ*t^(a - β)/cₐ) + cᵧ*t^(-γ)/gamma(1 - γ)
# end
# G_fraczener(t::Vector{T}, params::Vector{T}) where T<:Real = G_fraczener.(t, (params,))
#
# function J_fraczener(t::Vector{T}, params::Vector{T}) where T<:Real
#     cₐ, a, cᵦ, β, cᵧ, γ = params
#
#     Ĵ(s) = (1/s)*(cₐ*s^a + cᵦ*s^β)/(cₐ*s^a*cᵦ*s^β + cᵧ*s^γ*(cₐ*s^a + cᵦ*s^β))
#
#     # # method gives NaN at time 0.0
#     # J = InverseLaplace.ILt(s -> Ĵ(s))
#     # return [J(t_val) for t_val in t]
#
#     return InverseLaplace.talbotarr(s -> Ĵ(s), t)
# end
#
# function Gp_fraczener(ω::T, params::Vector{T}) where T<:Real
#     cₐ, a, cᵦ, β, cᵧ, γ = params
#
#     denominator = (cₐ*ω^a)^2 + (cᵦ*ω^β)^2 + 2*(cₐ*ω^a)*(cᵦ*ω^β)*cos((a - β)*π/2)
#
#     numerator = ((cᵦ*ω^β)^2)*(cₐ*ω^a)*cos(a*π/2) + ((cₐ*ω^a)^2)*(cᵦ*ω^β)*cos(β*π/2)
#
#     Gp = numerator/denominator + cᵧ*ω^γ*cos(γ*π/2)
# end
# Gp_fraczener(ω::Vector{T}, params::Vector{T}) where T<:Real = Gp_fraczener.(ω, (params,))
#
# function Gpp_fraczener(ω::T, params::Vector{T}) where T<:Real
#     cₐ, a, cᵦ, β, cᵧ, γ = params
#
#     denominator = (cₐ*ω^a)^2 + (cᵦ*ω^β)^2 + 2*(cₐ*ω^a)*(cᵦ*ω^β)*cos((a - β)*π/2)
#
#     numerator = ((cᵦ*ω^β)^2)*(cₐ*ω^a)*sin(a*π/2) + ((cₐ*ω^a)^2)*(cᵦ*ω^β)*sin(β*π/2)
#
#     Gpp = numerator/denominator + cᵧ*ω^γ*sin(γ*π/2)
# end
# Gpp_fraczener(ω::Vector{T}, params::Vector{T}) where T<:Real = Gpp_fraczener.(ω, (params,))
#
# """
#     FractionalZener([params::Vector{T}]) where T<:Real
#
# 2 springpots in series, both in parallel with a springpot.
# """
# FractionalZener(params::Vector{T}) where T<:Real = RheologyModel(G_fraczener, J_fraczener, Gp_fraczener, Gpp_fraczener, params, ["model created by user with parameters $params"])
# FractionalZener() = RheologyModel(G_fraczener, J_fraczener, Gp_fraczener, Gpp_fraczener, [1.0, 0.7, 1.0, 0.5, 1.0, 0.2], ["model created with default parameters"])

# # Fractional Standard Linear Solid in (Fractional) Maxwell Form
# function G_fracsls(t::T, params::Vector{T}) where T<:Real
#     cₐ, a, kᵦ, kᵧ = params
#
#     G = kᵦ*mittleff(a, -(kᵦ/cₐ)*t^a) + kᵧ
# end
# G_fracsls(t::Vector{T}, params::Vector{T}) where T<:Real = G_fracsls.(t, (params,))
#
# function J_fracsls(t::Vector{T}, params::Vector{T}) where T<:Real
#     cₐ, a, kᵦ, kᵧ = params
#
#     Ĵ(s) = (1/s)*(cₐ*s^a + kᵦ)/(cₐ*s^a*kᵦ + kᵧ*(cₐ*s^a + kᵦ))
#
#     # method gives NaN at time 0.0
#     # J = InverseLaplace.ILt(s -> Ĵ(s))
#     # return [J(t_val) for t_val in t]
#
#     return InverseLaplace.talbotarr(s -> Ĵ(s), t)
# end
#
# function Gp_fracsls(ω::T, params::Vector{T}) where T<:Real
#     cₐ, a, kᵦ, kᵧ = params
#
#     denominator = (cₐ*ω^a)^2 + kᵦ^2 + 2*(cₐ*ω^a)*kᵦ*cos(a*π/2)
#
#     numerator = kᵦ^2*(cₐ*ω^a)*cos(a*π/2) + ((cₐ*ω^a)^2)*kᵦ
#
#     Gp = numerator/denominator + kᵧ
# end
# Gp_fracsls(ω::Vector{T}, params::Vector{T}) where T<:Real = Gp_fracsls.(ω, (params,))
#
# function Gpp_fracsls(ω::T, params::Vector{T}) where T<:Real
#     cₐ, a, kᵦ, kᵧ = params
#
#     denominator = (cₐ*ω^a)^2 + kᵦ^2 + 2*(cₐ*ω^a)*kᵦ*cos(a*π/2)
#
#     numerator = kᵦ^2*(cₐ*ω^a)*sin(a*π/2)
#
#     Gpp = numerator/denominator
# end
# Gpp_fracsls(ω::Vector{T}, params::Vector{T}) where T<:Real = Gpp_fracsls.(ω, (params,))
#
# """
#     FractionalSLS([params::Vector{T}]) where T<:Real
#
# A springpot and spring in series, both in parallel with a spring.
# """
# FractionalSLS(params::Vector{T}) where T<:Real = RheologyModel(G_fracsls, J_fracsls, Gp_fracsls, Gpp_fracsls, params, ["model created by user with parameters $params"])
# FractionalSLS() = RheologyModel(G_fracsls, J_fracsls, Gp_fracsls, Gpp_fracsls, [0.5, 0.5, 1.0, 0.2], ["model created with default parameters"])

# # Standard Linear Solid in Maxwell Form
# function G_sls(t, params)
#     η, kᵦ, kᵧ = params
#
#     G = kᵧ + kᵦ*exp(-t*kᵦ/η)
# end
# G_sls(t::Vector{T}, params::Vector{T}) where {T<:Real} = G_sls.(t, (params,))
#
# function J_sls(t::T, params::Vector{T}) where T<:Real
#     η, kᵦ, kᵧ = params
#
#     c₀ = 1/kᵧ
#     c₁ = kᵦ/(kᵧ*(kᵧ + kᵦ))
#     τᵣ = η*(kᵧ + kᵦ)/(kᵧ*kᵦ)
#
#     J = c₀ - c₁*exp(-t/τᵣ)
# end
# J_sls(t::Vector{T}, params::Vector{T}) where T<:Real = J_sls.(t, (params,))
#
# function Gp_sls(ω::T, params::Vector{T}) where T<:Real
#     η, kᵦ, kᵧ = params
#
#     τ = η/kᵦ
#
#     denominator = 1 + τ^2*ω^2
#     numerator = ω^2*τ^2*kᵦ
#
#     Gp = numerator/denominator + kᵧ
# end
# Gp_sls(ω::Vector{T}, params::Vector{T}) where T<:Real = Gp_sls.(ω, (params,))
#
# function Gpp_sls(ω::T, params::Vector{T}) where T<:Real
#     η, kᵦ, kᵧ = params
#
#     τ = η/kᵦ
#
#     denominator = 1 + τ^2*ω^2
#     numerator = ω*τ*kᵦ
#
#     Gpp = numerator/denominator
# end
# Gpp_sls(ω::Vector{T}, params::Vector{T}) where T<:Real = Gpp_sls.(ω, (params,))
#
# """
#     SLS(params::Vector{T}) where T<:Real
#
# A spring and dashpot in series, both in parallel with a spring.
# """
# SLS(params::Vector{T}) where T<:Real = RheologyModel(G_sls, J_sls, Gp_sls, Gpp_sls, params, ["model created by user with parameters $params"])
# SLS() = RheologyModel(G_sls, J_sls, Gp_sls, Gpp_sls, [1.0, 0.5, 1.0], ["model created with default parameters"])
