#!/usr/bin/env julia

# Standard Linear Solid (Maxwell Form)
function G_SLS(t::Array{Real,1}, params::Array{Real,1})
    k₀, k₁, η₁ = params

    G = k₀ + k₁*exp.(-t*k₁/η₁)
end

function J_SLS(t::Array{Real,1}, params::Array{Real,1})
    c₀ = 1/k
    c₁ = cᵦ/(k*(k + cᵦ))
    τᵣ = β*(k + cᵦ)/(k*cᵦ)

    J = c₀ - c₁*exp.(-t/τᵣ)
end

SLS() = RheologyModel(G_SLS, J_SLS, [1.0, 0.5, 1.0], ["model created with default parameters"])
SLS(params::Array{Real, 1}) = RheologyModel(G_SLS, J_SLS, params, ["model created by user with parameters $params"])

# Spring-Pot Model
function G_springpot(t::Array{Real,1}, params::Array{Real,1})
    cᵦ, β = params

    G = cᵦ*t.^(-β)/gamma(1 - β)
end

function J_springpot(t::Array{Real,1}, params::Array{Real,1})
    cᵦ, β = params

    J = (t.^β)/(cᵦ*gamma(1 + β))
end

SpringPot() = RheologyModel(G_springpot, J_springpot, [2.0, 0.5], ["model created with default parameters"])
SpringPot(params::Array{Real, 1}) = RheologyModel(G_springpot, J_springpot, params, ["model created by user with parameters $params"])

# Fractional Maxwell Model
function G_fractmaxwell(t::Array{Real,1}, params::Array{Real,1})
    cₐ, a, cᵦ, β = params

    G = cᵦ*t.^(-β)*mittleff.(a - β, 1 - β, -cᵦ*t.^(a - β)/cₐ)
end

function J_fractmaxwell(t::Array{Real,1}, params::Array{Real,1})
    cₐ, a, cᵦ, β = params

    J = t.^(a)/(cₐ*gamma(1 + a)) + t.^(β)/(cᵦ*gamma(1 + β))
end

FractionalMaxwell() = RheologyModel(G_fractmaxwell, J_fractmaxwell, [2.0, 0.2, 1.0, 0.5], ["model created with default parameters"])
FractionalMaxwell(params::Array{Real, 1}) = RheologyModel(G_fractmaxwell, J_fractmaxwell, params, ["model created by user with parameters $params"])

# Fraction Kelvin-Voigt model
function G_fractKV(t::Array{Real,1}, params::Array{Real,1})
    cₐ, a, cᵦ, β = params

    G = cₐ*t.^(-a)/gamma(1 - a) + cᵦ*t.^(-β)/gamma(1 - β) 
end

function J_fractKV(t::Array{Real,1}, params::Array{Real,1})
    cₐ, a, cᵦ, β = params

    J = cₐ*t.^(a)*mittleff.(a - β, 1 + a, -cᵦ*t.^(a - β)/cₐ)
end

FractionalKelvinVoigt() = RheologyModel(G_fractKV, J_fractKV, [2.0, 0.2, 1.0, 0.5], ["model created with default parameters"])
FractionalKelvinVoigt(params::Array{Real, 1}) = RheologyModel(G_fractKV, J_fractKV, params, ["model created by user with parameters $params"])

# Fractional Special Model
function G_fractspecial(t::Array{Real,1}, params::Array{Real,1})::Array{Real,1}
    k, cᵦ, β, η = params

    G = k + cᵦ*t.^(-β).*mittleff.(1 - β, 1 - β, -cᵦ*(t.^(1 - β))/η)
end

function J_fractspecial_slow(t::Array{Real,1}, params::Array{Real,1})::Array{Real,1}
    k, cᵦ, β, η = params

    a = η/cᵦ
    b = k*η/cᵦ

    J = InverseLaplace.ILt( s -> 1/s^2 * ((1+a*s^(1-β)) / (η+k/s+b*s^(-β)) ))

    return [J(t_val) for t_val in t]
end

function J_fractspecial_fast(t::Array{Real,1}, params::Array{Real,1})::Array{Real,1}
    k, cᵦ, β, η = params

    a = η/cᵦ
    b = k*η/cᵦ

    J = InverseLaplace.talbotarr( s -> 1/s^2 * ((1+a*s^(1-β)) / (η+k/s+b*s^(-β)) ), t)
end
                
FractionalSpecial() = RheologyModel(G_fractspecial, J_fractspecial_fast, [1.0, 1.0, 0.2, 1.0], ["model created with default parameters"])
FractionalSpecial(params::Array{Real, 1}) = RheologyModel(G_fractspecial, J_fractspecial_fast, params, ["model created by user with parameters $params"])

FractionalSpecial_Precision() = RheologyModel(G_fractspecial, J_fractspecial_slow, [1.0, 1.0, 0.2, 1.0], ["model created with default parameters"])
FractionalSpecial_Precision(params::Array{Real, 1}) = RheologyModel(G_fractspecial, J_fractspecial_slow, params, ["model created by user with parameters $params"])

# function singularitytest(modulus::Function)

#     (symbol, params) = modeldatabase(modulus)

#     startval = modulus([0.0], params)[1]

#     if startval == NaN || startval == Inf
#         return true
#     else
#         return false 
#     end

# end
