#!/usr/bin/env julia

# Fraction Kelvin-Voigt model
function G_fractKV(t::Vector{T}, params::Vector{T}) where T<:Real
    cₐ, a, cᵦ, β = params

    G = cₐ*t.^(-a)/gamma(1 - a) + cᵦ*t.^(-β)/gamma(1 - β) 
end

function J_fractKV(t::Vector{T}, params::Vector{T}) where T<:Real
    cₐ, a, cᵦ, β = params

    J = (t.^(a)/cₐ).*mittleff.(a - β, 1 + a, -cᵦ*t.^(a - β)/cₐ)
end

function Gp_fractKV(ω::Vector{T}, params::Vector{T}) where T<:Real
    cₐ, a, cᵦ, β = params

    # cosine floating point error work-around
    if a!=1.0 && β!=1.0
        return cₐ*ω.^a*cos(a*π/2) + cᵦ*ω.^β*cos(β*π/2)
    elseif a==1.0 && β!=1.0
        return cᵦ*ω.^β*cos(β*π/2)
    elseif a==1.0 && β==1.0
        return [0.0 for i in ω]
    end
end

function Gpp_fractKV(ω::Vector{T}, params::Vector{T}) where T<:Real
    cₐ, a, cᵦ, β = params

    Gpp = cₐ*ω.^a*sin(a*π/2) + cᵦ*ω.^β*sin(β*π/2)
end

FractionalKelvinVoigt() = RheologyModel(G_fractKV, J_fractKV, Gp_fractKV, Gpp_fractKV, [2.0, 0.2, 1.0, 0.5], ["model created with default parameters"])
FractionalKelvinVoigt(params::Array{T, 1}) where T<:Real = RheologyModel(G_fractKV, J_fractKV, Gp_fractKV, Gpp_fractKV, params, ["model created by user with parameters $params"])

# Fraction Kelvin-Voigt model with lower spring-pot specialized to a spring
function G_fractKVspring(t::Vector{T}, params::Vector{T}) where T<:Real
    cₐ, a, k = params

    G = cₐ*t.^(-a)/gamma(1 - a) + k
end

function J_fractKVspring(t::Vector{T}, params::Vector{T}) where T<:Real
    cₐ, a, k = params

    J = (t.^(a)/cₐ).*mittleff.(a, 1 + a, -k*t.^a/cₐ)
end

function Gp_fractKVspring(ω::Vector{T}, params::Vector{T}) where T<:Real
    cₐ, a, k = params
    
    # cosine floating point error work-around
    if a!=1.0
        return cₐ*ω.^a*cos(a*π/2) + k
    else
        return [k for i in ω]
    end
end

function Gpp_fractKVspring(ω::Vector{T}, params::Vector{T}) where T<:Real
    cₐ, a, k = params

    Gpp = cₐ*ω.^a*sin(a*π/2)
end

FractionalKVspring() = RheologyModel(G_fractKVspring, J_fractKVspring, Gp_fractKVspring, Gpp_fractKVspring, [2.0, 0.2, 1.0], ["model created with default parameters"])
FractionalKVspring(params::Array{T, 1}) where T<:Real = RheologyModel(G_fractKVspring, J_fractKVspring, Gp_fractKVspring, Gpp_fractKVspring, params, ["model created by user with parameters $params"])

# Fraction Kelvin-Voigt model with upper spring-pot specialized to a dash-pot
function G_fractKVdashpot(t::Vector{T}, params::Vector{T}) where T<:Real
    η, cᵦ, β = params

    # G = η*t.^(-1)/gamma(1.0 - 1.0) + cᵦ*t.^(-β)/gamma(1 - β)

    # dashpot (first term) denominator is undefined as 0*∞, physical intuition
    # suggest an instantaneous singularity as an dashpot/spring model. As this
    # is already captured by the second power law term, just keep that one in.
    G = [i!=0.0 ? cᵦ*i^(-β)/gamma(1 - β) : Inf for i in t]
end

function J_fractKVdashpot(t::Vector{T}, params::Vector{T}) where T<:Real
    η, cᵦ, β = params

    J = (t/η).*mittleff.(1 - β, 1 + 1, -cᵦ*t.^(1.0 - β)/η)
end

function Gp_fractKVdashpot(ω::Vector{T}, params::Vector{T}) where T<:Real
    η, cᵦ, β = params

    # cosine floating point error work-around
    if β==1.0
        return [0.0 for i in ω]
    else
        return cᵦ*ω.^β*cos(β*π/2)
    end
end

function Gpp_fractKVdashpot(ω::Vector{T}, params::Vector{T}) where T<:Real
    η, cᵦ, β = params

    Gpp = η*ω + cᵦ*ω.^β*sin(β*π/2)
end

FractionalKVdashpot() = RheologyModel(G_fractKVdashpot, J_fractKVdashpot, Gp_fractKVdashpot, Gpp_fractKVdashpot, [2.0, 0.2, 1.0], ["model created with default parameters"])
FractionalKVdashpot(params::Array{T, 1}) where T<:Real = RheologyModel(G_fractKVdashpot, J_fractKVdashpot, Gp_fractKVdashpot, Gpp_fractKVdashpot, params, ["model created by user with parameters $params"])

# Standard Kelvin-Voigt model
function G_KV(t::Vector{T}, params::Vector{T}) where T<:Real
    η, k = params

    # G = k # + η*δ(t) in Findley et al. (1989)
    G = [i!=0.0 ? k : Inf for i in t]
end

function J_KV(t::Vector{T}, params::Vector{T}) where T<:Real
    η, k = params

    J = (1 - exp.(-k*t./η))/k
end

function Gp_KV(ω::Vector{T}, params::Vector{T}) where T<:Real
    η, k = params

    Gp = [k for i in ω]
end

function Gpp_KV(ω::Vector{T}, params::Vector{T}) where T<:Real
    η, k = params

    Gpp = η*ω
end

KelvinVoigt() = RheologyModel(G_KV, J_KV, Gp_KV, Gpp_KV, [2.0, 1.0], ["model created with default parameters"])
KelvinVoigt(params::Array{T, 1}) where T<:Real = RheologyModel(G_KV, J_KV, Gp_KV, Gpp_KV, params, ["model created by user with parameters $params"])
