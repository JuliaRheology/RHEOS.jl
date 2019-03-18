# 2 time scale Generalized Maxwell-Wiechert type models

SLS2 = fun_gen(
        # Model name
        name="SLS2",
        # Model parameters,
        p = [:G₀, :G₁, :η₁, :G₂, :η₂],
        # Relaxation modulus
        G = quote
              return G₀ .+ G₁*exp.(-t*G₁/η₁) .+ G₂*exp.(-t*G₂/η₂)
            end,
        # Creep modulus
        J = quote
              Jbar(s) = (1/s^2)*(G₀/s+G₁*1/(s+tau1)+G₂*1/(s+tau2))
              return invLaplace(s -> Jbar(s), t)
            end,
        # Storage modulus
        Gp = quote
               return [-1]
             end,
        # Loss modulus
        Gpp = quote
                return [-1]
              end,
        # Network
        info= "
                                       ___
                 ___╱╲  ╱╲  ╱╲  ________| |_____
                 |    ╲╱  ╲╱  ╲╱       _|_|     |
                 |      G₂              η₂      |
                 |                     ___      |
              ___|__╱╲  ╱╲  ╱╲  ________| |_____|___
                 |    ╲╱  ╲╱  ╲╱       _|_|     |
                 |      G₁              η₁      |
                 |                              |
                 |__________╱╲  ╱╲  ╱╲  ________|
                              ╲╱  ╲╱  ╲╱
                                G₀
               "
        )

#
# """
# SLS2(params::Vector{T}) where T<:Real
#
# A spring in parallel to two Maxwell branches (G' and G'' temporary missing).
# """
# SLS2 = RheoModelClass(G_SLS2, J_SLS2, Gp_SLS2, Gpp_SLS2, params_SLS2, info_SLS2)
#

# function G_SLS2(t::Union{Array{RheoFloat,1},RheoFloat}, params::Array{RheoFloat,1})
#
#    G₀, G₁, η₁, G₂, η₂ = params
#
#    return G₀ .+ G₁*exp.(-t*G₁/η₁) .+ G₂*exp.(-t*G₂/η₂)
#
# end



# function J_SLS2(t::Union{Array{RheoFloat,1},RheoFloat}, params::Array{RheoFloat,1})
#     G₀, G₁, η₁, G₂, η₂ = params
#
#     tau1 = G₁/η₁
#     tau2 = G₂/η₂
#
#     Jbar(s) = (1/s^2)*(G₀/s+G₁*1/(s+tau1)+G₂*1/(s+tau2))
#
#     return invLaplace(s -> Jbar(s), t)
# end


#
# Gp_SLS2 = null_modulus
#
# Gpp_SLS2 = null_modulus
