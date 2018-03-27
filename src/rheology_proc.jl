#!/usr/bin/env julia

"""
    modelfit!(self::RheologyType, model::String, params_init::Array{Float64,1}, low_bounds::Array{Float64,1}, hi_bounds::Array{Float64,1})

Fit RheologyType struct to model and store fitted parameters in self.fittedmodels.

# Arguments

- `self`: RheologyType struct containing all data
- `model`: E.g. "SLS", "springpot", "burgers" etc. See models.jl for full list
- `params_init`: Initial parameters to use in fit
- `low_bounds`: Lower bounds for parameters
- `hi_bounds`: Higher bounds for parameters
"""
function modelfit!(self::RheologyType,
                  model::String,
                  params_init::Array{Float64,1},
                  low_bounds::Array{Float64,1},
                  hi_bounds::Array{Float64,1})

    # generate time series difference array (for convolution)
    dt_series = deriv(self.t)

    # get appropriate RheologyModel struct
    modulus = moduli(model, self.test_type)

    # start fit
    (minf, minx, ret) = leastsquares_init(params_init, low_bounds, hi_bounds,
                                          modulus, self.t, dt_series, self.dcontrolled,
                                          self.measured; sampling = self.sampling,
                                          insight = self.insight)

    # store fit results in RheologyType struct's fittedmodels dictionary
    self.fittedmodels[model] = minx
end

"""
    modelcomplete!(self::RheologyData, model::String, params::Array{Float64,1})

Given partial data (just t and ϵ for "strlx" test or t and σ for "creep"),
model and parameters, find missing data (σ for "strlx" and ϵ for "creep").
Currently only works for RheologyData and not more general RheologyType.
"""
function modelcomplete!(self::RheologyData, model::String, params::Array{Float64,1})

    # generate time series difference array (for convolution)
    dt_series = deriv(self.t)

    # get modulus function
    modulus = moduli(model, self.test_type)

    # get convolution
    if self.sampling == "constant"
        convolved = boltzconvolve(modulus, self.t, dt_series, params, self.dcontrolled)
    elseif self.sampling == "variable"
        convolved = boltzintegral(modulus, self.t, params, self.dcontrolled)
    end

    # store output in RheologyData object
    self.measured = convolved
    
    if self.test_type == "strlx"
        self.σ = convolved
    elseif self.test_type == "creep"
        self.ϵ = convolved
    end
end
#= think about extending what is stored in fittedmodels dict, could be
params, convolved with those params, cost and more. Then modelcomplete!
should place this convolved there, not in main self =#
