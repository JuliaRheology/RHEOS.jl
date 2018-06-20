#!/usr/bin/env julia

"""
    modelfit!(self::RheologyData, model::String, params_init::Array{Float64,1}, low_bounds::Array{Float64,1}, hi_bounds::Array{Float64,1})

Fit RheologyData struct to model and store fitted parameters in self.fittedmodels.

# Arguments

- `self`: RheologyData struct containing all data
- `model`: E.g. "SLS", "springpot", "burgers" etc. See models.jl for full list
- `params_init`: Initial parameters to use in fit
- `low_bounds`: Lower bounds for parameters
- `hi_bounds`: Higher bounds for parameters
"""
function modelfit!(self::RheologyData,
                  model::String,
                  params_init::Array{Float64,1},
                  low_bounds::Array{Float64,1},
                  hi_bounds::Array{Float64,1})

    # generate time series difference array (for convolution)
    dt_series = deriv(self.t)

    # get appropriate RheologyModel struct
    modulus = moduli(model, self.test_type)
    # start fit
    tic()
    (minf, minx, ret) = leastsquares_init(params_init, low_bounds, hi_bounds,
                                          modulus, self.t, dt_series, self.dcontrolled,
                                          self.measured; sampling = self.sampling,
                                          insight = self.insight)
    timetaken = toc()
    # store fit results in RheologyData struct's fittedmodels dictionary
    self.fittedmodels[model] = (minf, minx, ret, timetaken)
end

"""
    modelcomplete!(self::RheologyData, model::String, params::Array{Float64,1})

Given partial data (just t and ϵ for "strlx" test or t and σ for "creep"),
model and parameters, find missing data (σ for "strlx" and ϵ for "creep").
Currently only works for RheologyData and not more general RheologyData.
"""
function modelcomplete!(self::RheologyData, modelname::String, params::Array{Float64,1})

    # generate time series difference array (for convolution)
    dt_series = deriv(self.t)

    # get model
    model = moduli(modelname, self.test_type)

    # get convolution
    if !model.singularity && self.sampling == "constant"
        convolved = boltzconvolve_nonsing(model, self.t, deriv(self.t), params, self.dcontrolled)

    elseif model.singularity && self.sampling == "constant"
        convolved = boltzconvolve_sing(model, self.t, deriv(self.t), params, self.dcontrolled)

    elseif !model.singularity && self.sampling == "variable"
        convolved = boltzintegral_nonsing(model, self.t, params, self.dcontrolled)

    elseif model.singularity && self.sampling == "variable"
        convolved = boltzintegral_sing(model, self.t, params, self.dcontrolled)

    end

    # store output in RheologyData object
    self.measured = convolved

    # store operation
    self.appliedops = vcat(self.appliedops, "modelname: $modelname, params: $params")

end
