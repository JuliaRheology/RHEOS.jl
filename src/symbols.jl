#!/usr/bin/env julia

# This file contains data and functions helping with the translation of symbols for the use of RHEOS without unicode symbols.


# translation table for symbols. left = symbol use by user, right = symbol used by model declaraction or RHEOS structs.
const symbol_convertion_table=(
    time = :t,
    t_col = :t,
    strain = :ϵ,
    epsilon = :ϵ,
    ϵ_col = :ϵ,
    ε = :ϵ,
    stress = :σ,
    sigma = :σ,
    σ_col = :σ,
    omega = :ω,
    ω_col = :ω,
    Gp_col = :Gp,
    Gpp_col = :Gpp,
    eta = :η,
    alpha = :α,
    beta = :β,
    gamma = :γ,
    c_alpha = :cₐ,
    c_beta = :cᵦ,
    c_gamma = :cᵧ,
    k_1 = :k₁,
    k_2 = :k₂,
    k_a = :kₐ,
    k_beta = :kᵦ,
    k_gamma = :kᵧ,
    eta_1 = :η₁,
    eta_2 = :η₂,
    eta_a = :ηₐ,
    eta_beta = :ηᵦ,
    eta_gamma = :ηᵧ,
    tau = :τ,
    G_phi = :Gᵩ,
    G_0 = :G₀,
    G_1 = :G₁,
    G_2 = :G₂
    )


function symbol_to_unicode(s::Symbol)    
    s in keys(symbol_convertion_table) ? symbol_convertion_table[s] : s 
end

function symbol_to_unicode(nt::NamedTuple)    
    NamedTuple{Tuple([ symbol_to_unicode(s) for s in keys(nt) ])}( values(nt) )
end

# Edge case useful for function taking option tuple parameters, e.g. modelfit
function symbol_to_unicode(nt::Nothing)    
    nothing 
end


function unicode_to_text(s::Symbol)
    r=findfirst(e->e==s, symbol_convertion_table)    
    isnothing(r) ? s : r
end

function unicode_to_text(nt::NamedTuple)
    NamedTuple{Tuple([ unicode_to_text(s) for s in keys(nt) ])}( values(nt) )
end
