#!/usr/bin/env julia

function forcing(t::Array{Float64,1}, t₀::Float64, t₁::Float64, τ::Float64)::Array{Float64,1}
	#forcing function, can be used to prescribe an
	#approximate step force or a step strain

	L = 1000.0 # amplitude
	k = 10/τ # convert from time scale to inverse

	L./(1+exp.(-k*(t - t₀))) - L./(1+exp.(-k*(t - t₁)))
end

function forcingdot(t::Array{Float64,1}, t₀::Float64, t₁::Float64, τ::Float64)::Array{Float64,1}

	deriv(forcing(t, t₀, t₁, τ), t)
end
