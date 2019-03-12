# 2 time scale Generalized Maxwell-Wiechert type models

params_SLS2 = [:G₀, :G₁, :η₁, :G₂, :η₂]

function info_SLS2()
    Base.print("
                               ___
         ___╱╲  ╱╲  ╱╲  ________| |_____
         |    ╲╱  ╲╱  ╲╱       _|_|     |
       __|      k1              η1      |____
         |             ___              |
         |______________| |_____________|
                       _|_|
                        η2

       ")
end

function G_SLS2(t::Union{Array{T,1},T}, params::Array{T,1}) where T<:Real

   G₀, G₁, η₁, G₂, η₂ = params

   G = G₀ .+ G₁*exp.(-t*G₁/η₁) .+ G₂*exp.(-t*G₂/η₂)

end

function J_SLS2(t::Union{Array{T,1},T}, params::Array{T,1}) where T<:Real
    G₀, G₁, η₁, G₂, η₂ = params

    tau1 = G₁/η₁
    tau2 = G₂/η₂

    Jbar(s) = (1/s^2)*(G₀/s+G₁*1/(s+tau1)+G₂*1/(s+tau2))

    return InverseLaplace.talbotarr(s -> Jbar(s), t)
end


Gp_SLS2 = null_modulus

Gpp_SLS2 = null_modulus


"""
SLS2(params::Vector{T}) where T<:Real

A spring in parallel to two Maxwell branches (G' and G'' temporary missing).
"""
#SLS2(params::Array{T, 1}) where T<:Real = RheologyModel(G_2SLS, J_2SLS, null_modulus, null_modulus, params, ["model created with default parameters"])
#SLS2() = RheologyModel(G_2SLS, J_2SLS, null_modulus, null_modulus, [1.0, 1.0, 1.0, 1.0, 1.0], ["model created with default parameters"])

#SLS2(params::Array{T, 1}) where T<:Real = RheologyModel(G_2SLS, J_2SLS, null_modulus, null_modulus, params, ["model created with default parameters"])
SLS2 = RheoModelClass(G_SLS2, J_SLS2, Gp_SLS2, Gpp_SLS2, params_SLS2, info_SLS2)
