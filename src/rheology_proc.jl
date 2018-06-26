#!/usr/bin/env julia

"""
    modelfit(self::RheologyData, modulus::Function[, p0::Array{Float64,1}, lo::Array{Float64,1}, hi::Array{Float64,1}; verbose::Bool = false])

Fit RheologyData struct to model and store fitted parameters in self.fittedmodels.

# Arguments

- `self`: RheologyData struct containing all data
- `modulus`: E.g. G_SLS, J_springpot etc. See base_models.jl for full list
- `p0`: Initial parameters to use in fit
- `lo`: Lower bounds for parameters
- `hi`: Higher bounds for parameters
- `verbose`: If true, prints parameters on each optimisation iteration
"""
function modelfit(self::RheologyData,
                  modulus::Function;
                  p0::Array{Float64,1} = [-1.0],
                  lo::Array{Float64,1} = [-1.0],
                  hi::Array{Float64,1} = [-1.0],
                  verbose:Bool = false)::RheologyModel

    # get modulus info from database
    (controlledvar, p0_default) = modeldatabase(modulus)

    # get singularity presence
    sing = singularitytest(modulus)

    # use default p0 if not provided
    if quasinull(p0)
        p0 = p0_default
    end 

    # generate time series difference array (for convolution)
    dt_series = deriv(self.t)

    # get derivative of controlled variable and measured variable
    if controlledvar == "σ"
        dcontrolled = deriv(self.σ, self.t)
        measured = self.ϵ
    elseif controlledvar == "ϵ"
        dcontrolled = deriv(self.ϵ, self.t)
        measured = self.σ
    end

    # start fit
    tic()
    (minf, minx, ret) = leastsquares_init(p0,
                                          lo, 
                                          hi,
                                          modulus, 
                                          self.t, 
                                          dt_series, 
                                          dcontrolled,
                                          measured; 
                                          insight = verbose,
                                          sampling = self.sampling, 
                                          singularity = sing)
    timetaken = toc()

    # store fit results in RheologyData struct's fittedmodels dictionary
    (minf, minx, ret, timetaken)
    
end

"""
    modelpredict(self::RheologyData, model::String, params::Array{Float64,1})

Given partial data (just t and ϵ for "strlx" test or t and σ for "creep"),
model and parameters, find missing data (σ for "strlx" and ϵ for "creep").
Currently only works for RheologyData and not more general RheologyData.
"""
function modepredict(self::RheologyData, modelname::String, params::Array{Float64,1})::RheologyData

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
