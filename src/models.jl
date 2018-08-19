#!/usr/bin/env julia

# Standard Linear Solid (Maxwell Form)
# Parameters in order: k₀, k₁, η₁
function G_SLS(t::Array{Real,1}, params::Array{Real,1})
    G = params[1] + params[2]*exp.(-t*params[2]/params[3])
end

function J_SLS(t::Array{Real,1}, params::Array{Real,1})
    c₀ = 1/params[1]
    c₁ = params[2]/(params[1]*(params[1] + params[2]))
    τᵣ = params[3]*(params[1] + params[2])/(params[1]*params[2])

    J = c₀ - c₁*exp.(-t/τᵣ)
end

SLS() = RheologyModel(G_SLS, J_SLS, [1.0, 0.5, 1.0], ["model created with default parameters"])
SLS(params::Array{Real, 1}) = RheologyModel(G_SLS, J_SLS, params, ["model created by user with parameters $params"])

# Spring-Pot Model
# Parameters in order: cᵦ, β
function G_springpot(t::Array{Real,1}, params::Array{Real,1})::Array{Real,1}
    G = params[1]*t.^(-params[2])/gamma(1-params[2])
end

function J_springpot(t::Array{Real,1}, params::Array{Real,1})::Array{Real,1}
    J = (t.^params[2])/(params[1]*gamma(1+params[2]))
end

SpringPot() = RheologyModel(G_springpot, J_springpot, [2.0, 0.5], ["model created with default parameters"])
SpringPot(params::Array{Real, 1}) = RheologyModel(G_springpot, J_springpot, params, ["model created by user with parameters $params"])

# Fractional Maxwell Model
# Parameters in order: cₐ, a, cᵦ, β
function G_fractmaxwell(t::Array{Real,1}, params::Array{Real,1})::Array{Real,1}

    G = 2*params[1]*mittleff(params[3], -params[2]*(t.^params[3]))
end


# Fractional Kelvin-Voigt Model
function G_fractKV(t::Array{Real,1}, params::Array{Real,1})::Array{Real,1}

    G = 2*params[1]*(1+(params[2]*t.^(-params[3]))/gamma(1-params[3]))
end





"""
    G_fractzener(t::Array{Real,1}, params::Array{Real,1})

Fraction Zener model as defined by Mainardi & Spada (2011) parameters by index
order: μ*, r, cᵦ, ν
"""
function G_fractzener(t::Array{Real,1}, params::Array{Real,1})::Array{Real,1}

    G = 2*params[1]*(1 + params[2]*mittleff(params[4], -params[3]*(t.^params[4])))
end

"""
    G_fractspecial(t::Array{Real,1}, params::Array{Real,1})

Special fractional model defined by A. Bonfanti (2017) defined as a 1 spring in
parallel with a (spring-pot and dash-pot in series) parameters by index order:
K, Cᵦ, β, η
"""
function G_fractspecial(t::Array{Real,1}, params::Array{Real,1})::Array{Real,1}

    G = params[1] + params[2]*(t.^(-params[3])).*mittleff(1-params[3],1-params[3],-params[2]*(t.^(1-params[3]))/params[4])
end



"""
    J_fractKV(t::Array{Real,1}, params::Array{Real,1})

Fractional Kelvin-Voigt model as defined by Mainardi & Spada (2011) parameters
by index order: μ, ν, τ
"""
function J_fractKV(t::Array{Real,1}, params::Array{Real,1})::Array{Real,1}

    J = (1-mittleff(params[2], -((t/params[3]).^params[2])))/(2*params[1])
end

"""
    J_fractmaxwell(t::Array{Real,1}, params::Array{Real,1})

Fractional Maxwell model as defined by Mainardi & Spada (2011) parameters by
index order: μ, ν, τ
"""
function J_fractmaxwell(t::Array{Real,1}, params::Array{Real,1})::Array{Real,1}

    J = (1 + ((t/params[3]).^params[2])/gamma(1+params[2]))/(2*params[1])
end

"""
    J_fractzener(t::Array{Real,1}, params::Array{Real,1})

Fraction Zener model as defined by Mainardi & Spada (2011) parameters by index
order: μ₁, μ₂, ν, τ₂
"""
function J_fractzener(t::Array{Real,1}, params::Array{Real,1})::Array{Real,1}

    J = (1/params[1] + (1-mittleff(params[3], -(t/params[4]).^params[3]))/params[2])/2
end

"""
    J_fractspecial_slow(t::Array{Real,1}, params::Array{Real,1})

Special fractional model defined by A. Bonfanti (2017) defined as a 1 spring in
parallel with a (spring-pot and dash-pot in series) parameters by index order:
K, Cᵦ, β, η.
Using default Talbot() method for the Laplace transform.
"""
function J_fractspecial_slow(t::Array{Real,1}, params::Array{Real,1})::Array{Real,1}
    a = params[4]/params[2]
    b = params[1]*params[4]/params[2]

    J = InverseLaplace.ILt( s -> 1/s^2 * ((1+a*s^(1-params[3])) / (params[4]+params[1]/s+b*s^(-params[3])) ))

    return [J(t_val) for t_val in t]
end

"""
    J_fractspecial_fast(t::Array{Real,1}, params::Array{Real,1})

Special fractional model defined by A. Bonfanti (2017) defined as a 1 spring in
parallel with a (spring-pot and dash-pot in series) parameters by index order:
K, Cᵦ, β, η.
Using the 'fixed' Talbot method for the Laplace transform (less accurate).
"""
function J_fractspecial_fast(t::Array{Real,1}, params::Array{Real,1})::Array{Real,1}
    a = params[4]/params[2]
    b = params[1]*params[4]/params[2]

    J = InverseLaplace.talbotarr( s -> 1/s^2 * ((1+a*s^(1-params[3])) / (params[4]+params[1]/s+b*s^(-params[3])) ), t)
end

####################
#~ Model Database ~#
####################

"""
    modeldatabase(modulus::Function)

For a specified modulus function, yields controlled variable, singularity presence
and suggested initial conditions as a tuple.

Additional models can be added at any time by the user by defining the modulus as
a function in "RHEOS/src/base_models.jl" and appending it the dictionary in this
function along with the its controlled variable (e.g. σ for creep moduli), singularity
presence and suggested initial conditions.
"""
function modeldatabase(modulus::Function)

    database = Dict(# creep moduli
                    J_springpot => ("σ", [2.0, 0.5]),
                    G_springpot => ("ϵ",),
                    J_fractKV => ("σ", [0.5, 0.65, 1.5]),
                    J_fractmaxwell => ("σ", [1.0, 0.35, 10.5]),
                    J_fractzener => ("σ", [1.0, 1.0, 0.7, 1.0]),
                    # Laplaced - takes relaxation modulus parameters but returns creep modulus
                    J_fractspecial_slow => ("σ", [1.0, 1.0, 0.2, 1.0]),
                    J_fractspecial_fast => ("σ", [1.0, 1.0, 0.2, 1.0]),


                    G_fractKV => ("ϵ",[0.25, 1.0, 0.25]),
                    G_fractmaxwell => ("ϵ", [0.5, 2.0, 0.52]),
                    G_fractzener => ("ϵ", [0.5, 1.0, 1.0, 0.5]),
                    G_fractspecial => ("ϵ",[1.0, 1.0, 0.2, 1.0]) 
                    
                    )

        database[modulus]

end

function singularitytest(modulus::Function)

    (symbol, params) = modeldatabase(modulus)

    startval = modulus([0.0], params)[1]

    if startval == NaN || startval == Inf
        return true
    else
        return false 
    end

end
