#!/usr/bin/env julia

"""
    modelfit(data::RheologyData, modulus::Function[, p0::Array{Float64,1}, lo::Array{Float64,1}, hi::Array{Float64,1}; verbose::Bool = false])

Fit RheologyData struct to model and return a fitted model as a RheologyModel object.

# Arguments

- `data`: RheologyData struct containing all data
- `modulus`: E.g. G_SLS, J_springpot etc. See base_models.jl for full list
- `p0`: Initial parameters to use in fit
- `lo`: Lower bounds for parameters
- `hi`: Higher bounds for parameters
- `verbose`: If true, prints parameters on each optimisation iteration
"""
function modelfit(data::RheologyData,
                  modulus::Function;
                  p0::Array{Float64,1} = [-1.0],
                  lo::Array{Float64,1} = [-1.0],
                  hi::Array{Float64,1} = [-1.0],
                  verbose::Bool = false)::RheologyModel

    # get modulus info from database
    (controlledvar, p0_default) = modeldatabase(modulus)

    # get singularity presence
    sing = singularitytest(modulus)

    # use default p0 if not provided
    if quasinull(p0)
        p0 = p0_default
    end 

    # generate time series difference array (for convolution)
    dt_series = deriv(data.t)

    # get derivative of controlled variable and measured variable
    if controlledvar == "σ"
        dcontrolled = deriv(data.σ, data.t)
        measured = data.ϵ
    elseif controlledvar == "ϵ"
        dcontrolled = deriv(data.ϵ, data.t)
        measured = data.σ
    end

    # start fit
    tic()
    (minf, minx, ret) = leastsquares_init(p0,
                                          lo, 
                                          hi,
                                          modulus, 
                                          data.t, 
                                          dt_series, 
                                          dcontrolled,
                                          measured; 
                                          insight = verbose,
                                          sampling = data.sampling, 
                                          singularity = sing)
    timetaken = toq()

    #
    modulusname = string(:modulus)

    log = vcat(data.log, "Fitted modulus($modulusname) in $timetaken, finished due to $ret with parans $minx")

    RheologyModel(modulus, minx, log)

end

"""
    modelpredict(data::RheologyData, model::RheologyModel)

Given data and model and parameters, predict new dataset based on both.
"""
function modelpredict(data::RheologyData, model::RheologyModel)::RheologyData

    # get modulus info from database
    (controlledvar, p0_default) = modeldatabase(model.modulus)

    # get singularity presence
    sing = singularitytest(model.modulus)

    if controlledvar == "σ"
        dcontrolled = deriv(data.σ, data.t)
    elseif controlledvar == "ϵ"
        dcontrolled = deriv(data.ϵ, data.t)
    end

    # generate time series difference array (for convolution)
    dt_series = deriv(data.t)

    # get convolution
    if !sing && data.sampling == "constant"
        convolved = boltzconvolve_nonsing(model.modulus, data.t, deriv(data.t), model.parameters, dcontrolled)
    elseif sing && data.sampling == "constant"
        convolved = boltzconvolve_sing(model.modulus, data.t, deriv(data.t), model.parameters, dcontrolled)

    elseif !sing && data.sampling == "variable"
        convolved = boltzintegral_nonsing(model.modulus, data.t, model.parameters, dcontrolled)

    elseif sing && data.sampling == "variable"
        convolved = boltzintegral_sing(model.modulus, data.t, model.parameters, dcontrolled)

    end

    if controlledvar == "σ"
        σ = data.σ
        ϵ = convolved
    elseif controlledvar == "ϵ"
        σ = convolved
        ϵ = data.ϵ
    end   

    # store operation
    log = vcat(data.log, "Predicted data from model:", model.log)

    RheologyData(σ, ϵ, data.t, data.sampling, log)

end
