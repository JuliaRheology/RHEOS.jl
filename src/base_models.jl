#!/usr/bin/env julia

#########################
#~ G Relaxation Moduli ~#
#########################

"""
    G_SLS(t::Array{Float64,1}, params::Array{Float64,1})

Standard Linear Solid / Zener Model with 1 time scale. Parameters by index order:
Gₑ, G₁, τ₁
"""
function G_SLS(t::Array{Float64,1}, params::Array{Float64,1})::Array{Float64,1}

    G = params[1] + params[2]*exp.(-t/params[3])
end

"""
    G_SLS2(t::Array{Float64,1}, params::Array{Float64,1})

Standard Linear Solid / Zener Model with 2 time scales, also known as Wiechert.
model. Parameters by index order: Gₑ, G₁, τ₁, G₂, τ₂
"""
function G_SLS2(t::Array{Float64,1}, params::Array{Float64,1})::Array{Float64,1}

    G = params[1] + params[2]*exp.(-t/params[3]) + params[4]*exp.(-t/params[5])
end

"""
    G_burgers(t::Array{Float64,1}, params::Array{Float64,1})

Burgers model. 2 time scales, relaxes to 0. Parameters by index order: G₁, τ₁, G₂, τ₂
"""
function G_burgers(t::Array{Float64,1}, params::Array{Float64,1})::Array{Float64,1}

    G = params[1]*exp.(-t/params[2]) + params[3]*exp.(-t/params[4])
end

"""
    G_springpot(t::Array{Float64,1}, params::Array{Float64,1})

Single Spring-pot model. Parameters by index order: cᵦ, ν
"""
function G_springpot(t::Array{Float64,1}, params::Array{Float64,1})::Array{Float64,1}

    G = params[1]*t.^(-params[2])/gamma(1-params[2])
end

"""
    G_fractKV(t::Array{Float64,1}, params::Array{Float64,1})

Fractional Kelvin-Voigt model as defined by Mainardi & Spada (2011) parameters
by index order: μ, cᵦ, ν
"""
function G_fractKV(t::Array{Float64,1}, params::Array{Float64,1})::Array{Float64,1}

    G = 2*params[1]*(1+(params[2]*t.^(-params[3]))/gamma(1-params[3]))
end

"""
    G_fractmaxwell(t::Array{Float64,1}, params::Array{Float64,1})

Fractional Maxwell model as defined by Mainardi & Spada (2011) parameters by
index order: μ, cᵦ, ν
"""
function G_fractmaxwell(t::Array{Float64,1}, params::Array{Float64,1})::Array{Float64,1}

    G = 2*params[1]*mittleff(params[3], -params[2]*(t.^params[3]))
end

"""
    G_fractzener(t::Array{Float64,1}, params::Array{Float64,1})

Fraction Zener model as defined by Mainardi & Spada (2011) parameters by index
order: μ*, r, cᵦ, ν
"""
function G_fractzener(t::Array{Float64,1}, params::Array{Float64,1})::Array{Float64,1}

    G = 2*params[1]*(1 + params[2]*mittleff(params[4], -params[3]*(t.^params[4])))
end

"""
    G_fractspecial(t::Array{Float64,1}, params::Array{Float64,1})

Special fractional model defined by A. Bonfanti (2017) defined as a 1 spring in
parallel with a (spring-pot and dash-pot in series) parameters by index order:
K, Cᵦ, β, η
"""
function G_fractspecial(t::Array{Float64,1}, params::Array{Float64,1})::Array{Float64,1}

    G = params[1] + params[2]*(t.^(-params[3])).*mittleff(1-params[3],1-params[3],-params[2]*(t.^(1-params[3]))/params[4])
end

####################
#~ J Creep Moduli ~#
####################

"""
    J_SLS(t::Array{Float64,1}, params::Array{Float64,1})

Standard Linear Solid / Zener Model with 1 time scale parameters by index order:
Jₑ, J₁, τ₁
"""
function J_SLS(t::Array{Float64,1}, params::Array{Float64,1})::Array{Float64,1}

    J = params[1] - params[2]*exp.(-t/params[3])
end

"""
    J_SLS2(t::Array{Float64,1}, params::Array{Float64,1})

Standard Linear Solid / Zener Model with 2 time scales, also known as Wiechert
model parameters by index order: Jₑ, J₁, τ₁, J₂, τ₂
"""
function J_SLS2(t::Array{Float64,1}, params::Array{Float64,1})::Array{Float64,1}

    J = params[1] - params[2]*exp.(-t/params[3]) - params[4]*exp.(-t/params[5])
end

"""
    J_burgers(t::Array{Float64,1}, params::Array{Float64,1})

Burgers model. Parameters by index order: Jₐ, Jᵦ, J₁, τ₁
"""
function J_burgers(t::Array{Float64,1}, params::Array{Float64,1})::Array{Float64,1}

    J = params[1] + params[2]*t - params[3]*exp.(-t/params[4])
end

"""
    J_springpot(t::Array{Float64,1}, params::Array{Float64,1})

Single Spring-pot model parameters by index order: Cᵦ, β
"""
function J_springpot(t::Array{Float64,1}, params::Array{Float64,1})::Array{Float64,1}

    J = (t.^params[2])/(params[1]*gamma(1+params[2]))
end

"""
    J_fractKV(t::Array{Float64,1}, params::Array{Float64,1})

Fractional Kelvin-Voigt model as defined by Mainardi & Spada (2011) parameters
by index order: μ, ν, τ
"""
function J_fractKV(t::Array{Float64,1}, params::Array{Float64,1})::Array{Float64,1}

    J = (1-mittleff(params[2], -((t/params[3]).^params[2])))/(2*params[1])
end

"""
    J_fractmaxwell(t::Array{Float64,1}, params::Array{Float64,1})

Fractional Maxwell model as defined by Mainardi & Spada (2011) parameters by
index order: μ, ν, τ
"""
function J_fractmaxwell(t::Array{Float64,1}, params::Array{Float64,1})::Array{Float64,1}

    J = (1 + ((t/params[3]).^params[2])/gamma(1+params[2]))/(2*params[1])
end

"""
    J_fractzener(t::Array{Float64,1}, params::Array{Float64,1})

Fraction Zener model as defined by Mainardi & Spada (2011) parameters by index
order: μ₁, μ₂, ν, τ₂
"""
function J_fractzener(t::Array{Float64,1}, params::Array{Float64,1})::Array{Float64,1}

    J = (1/params[1] + (1-mittleff(params[3], -(t/params[4]).^params[3]))/params[2])/2
end

"""
    J_fractspecial_slow(t::Array{Float64,1}, params::Array{Float64,1})

Special fractional model defined by A. Bonfanti (2017) defined as a 1 spring in
parallel with a (spring-pot and dash-pot in series) parameters by index order:
K, Cᵦ, β, η.
Using default Talbot() method for the Laplace transform.
"""
function J_fractspecial_slow(t::Array{Float64,1}, params::Array{Float64,1})::Array{Float64,1}
    a = params[4]/params[2]
    b = params[1]*params[4]/params[2]

    J = InverseLaplace.ILt( s -> 1/s^2 * ((1+a*s^(1-params[3])) / (params[4]+params[1]/s+b*s^(-params[3])) ))

    return [J(t_val) for t_val in t]
end

"""
    J_fractspecial_fast(t::Array{Float64,1}, params::Array{Float64,1})

Special fractional model defined by A. Bonfanti (2017) defined as a 1 spring in
parallel with a (spring-pot and dash-pot in series) parameters by index order:
K, Cᵦ, β, η.
Using the 'fixed' Talbot method for the Laplace transform (less accurate).
"""
function J_fractspecial_fast(t::Array{Float64,1}, params::Array{Float64,1})::Array{Float64,1}
    a = params[4]/params[2]
    b = params[1]*params[4]/params[2]

    J = InverseLaplace.talbotarr( s -> 1/s^2 * ((1+a*s^(1-params[3])) / (params[4]+params[1]/s+b*s^(-params[3])) ), t)
end

