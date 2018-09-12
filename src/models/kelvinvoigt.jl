#!/usr/bin/env julia

# Fraction Kelvin-Voigt model
function G_fractKV(t::Array{T,1}, params::Array{T,1}) where T<:Real
    cₐ, a, cᵦ, β = params

    G = cₐ*t.^(-a)/gamma(1 - a) + cᵦ*t.^(-β)/gamma(1 - β) 
end

function J_fractKV(t::Array{T,1}, params::Array{T,1}) where T<:Real
    cₐ, a, cᵦ, β = params

    J = cₐ*t.^(a)*mittleff.(a - β, 1 + a, -cᵦ*t.^(a - β)/cₐ)
end

function Gp_fractKV(ω::Array{T,1}, params::Array{T,1}) where T<:Real
    cₐ, a, cᵦ, β = params

    Gp = cₐ*ω.^a*cos(a*π/2) + cᵦ*ω.^β*cos(β*π/2)
end

function Gpp_fractKV(ω::Array{T,1}, params::Array{T,1}) where T<:Real
    cₐ, a, cᵦ, β = params

    Gpp = cₐ*ω.^a*sin(a*π/2) + cᵦ*ω.^β*sin(β*π/2)
end

FractionalKelvinVoigt() = RheologyModel(G_fractKV, J_fractKV, Gp_fractKV, Gpp_fractKV, [2.0, 0.2, 1.0, 0.5], ["model created with default parameters"])
FractionalKelvinVoigt(params::Array{T, 1}) where T<:Real = RheologyModel(G_fractKV, J_fractKV, Gp_fractKV, Gpp_fractKV, params, ["model created by user with parameters $params"])

# Fraction Kelvin-Voigt model with lower spring-pot specialized to a spring
function G_fractKVspring(t::Array{T,1}, params::Array{T,1}) where T<:Real
    cₐ, a, k = params

    G = cₐ*t.^(-a)/gamma(1 - a) + k
end

function J_fractKVspring(t::Array{T,1}, params::Array{T,1}) where T<:Real
    cₐ, a, k = params

    J = cₐ*t.^(a)*mittleff.(a - β, 1 + a, -cᵦ*t.^(a - β)/cₐ)
end

function Gp_fractKVspring(ω::Array{T,1}, params::Array{T,1}) where T<:Real
    cₐ, a, k = params

    Gp = cₐ*ω.^a*cos(a*π/2) + k
end

function Gpp_fractKVspring(ω::Array{T,1}, params::Array{T,1}) where T<:Real
    cₐ, a, k = params

    Gpp = cₐ*ω.^a*sin(a*π/2)
end

FractionalKVspring() = RheologyModel(G_fractKVspring, J_fractKVspring, Gp_fractKVspring, Gpp_fractKVspring, [2.0, 0.2, 1.0], ["model created with default parameters"])
FractionalKVspring(params::Array{T, 1}) where T<:Real = RheologyModel(G_fractKVspring, J_fractKVspring, Gp_fractKVspring, Gpp_fractKVspring, params, ["model created by user with parameters $params"])

# Fraction Kelvin-Voigt model with upper spring-pot specialized to a dash-pot
function G_fractKVdashpot(t::Array{T,1}, params::Array{T,1}) where T<:Real
    η, cᵦ, β = params

    # warn("Using G_fractKVdashpot, a single spring-pot or fractional KV with a~0.99 may be more effective.")

    # note, dashpot does not contribute anything after t=0.0 for step input
    # and first is singularity so should just use a spring-pot
    G = η*t.^(-1)/gamma(1.0 - 1.0) + cᵦ*t.^(-β)/gamma(1 - β)
end

function J_fractKVdashpot(t::Array{T,1}, params::Array{T,1}) where T<:Real
    η, cᵦ, β = params

    J = η*t.*mittleff.(1 - β, 1 + 1, -cᵦ*t.^(1.0 - β)/η)
end

function Gp_fractKVdashpot(ω::Array{T,1}, params::Array{T,1}) where T<:Real
    η, cᵦ, β = params

    Gp = cᵦ*ω.^β*cos(β*π/2)
end

function Gpp_fractKVdashpot(ω::Array{T,1}, params::Array{T,1}) where T<:Real
    η, cᵦ, β = params

    Gpp = η*ω + cᵦ*ω.^β*sin(β*π/2)
end

FractionalKVdashpot() = RheologyModel(G_fractKVdashpot, J_fractKVdashpot, Gp_fractKVdashpot, Gpp_fractKVdashpot, [2.0, 0.2, 1.0], ["model created with default parameters"])
FractionalKVdashpot(params::Array{T, 1}) where T<:Real = RheologyModel(G_fractKVdashpot, J_fractKVdashpot, Gp_fractKVdashpot, Gpp_fractKVdashpot, params, ["model created by user with parameters $params"])

# Standard Kelvin-Voigt model
function G_KV(t::Array{T,1}, params::Array{T,1}) where T<:Real
    η, k = params

    warn("Unimplemented singularity in Kelvin-Voigt Relaxation Modulus")

    G = k # + η*δ(t) in Findley et al. (1989)
end

function J_KV(t::Array{T,1}, params::Array{T,1}) where T<:Real
    η, k = params

    J = (1 - exp.(-k*t./η))/k
end

function Gp_KV(ω::Array{T,1}, params::Array{T,1}) where T<:Real
    η, k = params

    Gp = k
end

function Gpp_KV(ω::Array{T,1}, params::Array{T,1}) where T<:Real
    η, k = params

    Gpp = η*ω
end

KelvinVoigt() = RheologyModel(G_KV, J_KV, Gp_KV, Gpp_KV, [2.0, 1.0], ["model created with default parameters"])
KelvinVoigt(params::Array{T, 1}) where T<:Real = RheologyModel(G_KV, J_KV, Gp_KV, Gpp_KV, params, ["model created by user with parameters $params"])
