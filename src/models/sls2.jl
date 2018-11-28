# 2 time scale Generalized Maxwell-Wiechert type models
function G_2SLS(t::Array{T,1}, params::Array{T,1}) where T<:Real

   G₀, G₁, η₁, G₂, η₂ = params

   G = G₀ .+ G₁*exp.(-t*G₁/η₁) .+ G₂*exp.(-t*G₂/η₂)

end

function J_2SLS(t::Array{T,1}, params::Array{T,1}) where T<:Real
    G₀, G₁, η₁, G₂, η₂ = params

    tau1 = G₁/η₁
    tau2 = G₂/η₂

    Jbar(s) = (1/s^2)*(G₀/s+G₁*1/(s+tau1)+G₂*1/(s+tau2))

    return InverseLaplace.talbotarr(s -> Jbar(s), t)
end

"""
    SLS2(params::Vector{T}) where T<:Real

A spring in parallel to two Maxwell branches (G' and G'' temporary missing).
"""
SLS2() = RheologyModel(G_2SLS, J_2SLS, null_modulus, null_modulus, [1.0, 1.0, 1.0, 1.0, 1.0], ["model created with default parameters"])
SLS2(params::Array{T, 1}) where T<:Real = RheologyModel(G_2SLS, J_2SLS, null_modulus, null_modulus, params, ["model created with default parameters"])
