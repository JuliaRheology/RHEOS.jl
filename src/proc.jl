#!/usr/bin/env julia

"""
    boltzintegral(model::Function, time_series::Array{Float64,1}, params::Array{Float64,1}, prescribed_dot::Array{Float64,1})

Calculate Boltzmann Superposition integral using direct integration method.

This is much slower and slightly less accurate (depending on sample resolution)
than the convolution method. However, it works for variable sample rate.

The singularity keyword argument should be set to 'true' if the viscoelastic
modulus used has a singularity at t = 0.0. If this is the case, returns an array
of length(time_series) - 1. Should be compared with [2:end] of reference array
when fitting.

# Arguments

- `model`: RheologyModel struct containing viscoelastic modulus function
- `time_series`: The array of times
- `params`: Parameters passed to viscoelastic modulus
- `prescribed_dot`: Derivative of (usually prescribed) variable inside the integration kernel
"""
function boltzintegral(model::RheologyModel, time_series::Array{Float64,1}, params::Array{Float64,1},
                    prescribed_dot::Array{Float64,1})::Array{Float64,1}

    # check whether using model with t -> 0, model -> ∞
    if model.singularity
        # singularity, so disregard 1st element of length
        I = zeros(length(time_series)-1)
        # only traverse after time t=0, first element
        for (i,v) in enumerate(time_series[2:end])
            # generate integral for each time step
            @inbounds τ = time_series[1:i]
            Modulus_arg = v - τ
            Modulusᵢ = model.form(Modulus_arg, params)
            @inbounds df_dtᵢ = prescribed_dot[1:i]
            intergrand = Modulusᵢ.*df_dtᵢ
            @inbounds I[i] = trapz(intergrand, τ)
        end
        return I

    else
        I = zeros(length(time_series))
        for (i,v) in enumerate(time_series[1:end])
            # generate integral for each time step
            @inbounds τ = time_series[1:i]
            Modulus_arg = v - τ
            Modulusᵢ = model.form(Modulus_arg, params)
            @inbounds df_dtᵢ = prescribed_dot[1:i]
            intergrand = Modulusᵢ.*df_dtᵢ
            @inbounds I[i] = trapz(intergrand, τ)
        end
        return I
    end
end

"""
    boltzconvolve(model::Function, time_series::Array{Float64,1}, dt_series::Array{Float64,1}, params::Array{Float64,1}, prescribed_dot::Array{Float64,1})

Calculate Boltzmann Superposition integral using convolution method.

This is much faster and slightly more accurate (depending on sample resolution)
than the integral method. However, it works for constant sample rate.

The singularity keyword argument should be set to 'true' if the viscoelastic
modulus used has a singularity at t = 0.0. If this is the case, returns an array
of length(time_series) - 1. Should be compared with [2:end] of reference array
when fitting.

# Arguments

- `model`: RheologyModel struct containing viscoelastic modulus function
- `time_series`: The array of times
- `dt_series`: Array of gradient of time series, can be found using deriv(time_series)
- `params`: Parameters passed to viscoelastic modulus
- `prescribed_dot`: Derivative of (usually prescribed) variable inside the integration kernel
"""
function boltzconvolve(model::RheologyModel, time_series::Array{Float64,1}, dt_series::Array{Float64,1},
                        params::Array{Float64,1}, prescribed_dot::Array{Float64,1})::Array{Float64,1}

    # check whether using model with t -> 0, model -> ∞
    if model.singularity
        # if so, convolved length will be original_length-1
        len = length(dt_series)-1
        Modulus = model.form(time_series, params)
        # fast convolution, ignoring initial singularity
        β = FastConv.convn(Modulus[2:end], prescribed_dot[1:end])
        # pick out relevant elements (1st half) and multiply by dt
        β = β[1:len].*dt_series[1:len]

    else
        Modulus = model.form(time_series, params)
        β = FastConv.convn(Modulus, prescribed_dot)
        # pick out relevant elements (1st half) and multiply by dt
        β = β[1:length(dt_series)].*dt_series
    end

    return β
end

"""
    objectivefunc(params::Array{Float64,1}, grad::Array{Float64,1}, model, time_series, dt_series, prescribed_dot, measured; _singularity = false, _sampling = "constant", _insight = false)

Generate the sum-of-squares of the difference between `measured` and the Boltzmann
convolution integral of the selected `model` and `prescribed_dot`.

# Arguments

- `params`: Array of parameters sent to `model`
- `grad`: Gradient argument used by NLOpt
- `model`: RheologyModel struct containing viscoelastic modulus function
- `time_series`: Array of time data
- `dt_series`: Array of time data differences,  equal to deriv(time_series)
- `prescribed_dot`: Convolved with `model` over `time_series`, usually dϵ/dt for stress relaxation and dσ/dt for creep
- `measured`: Data for comparison against, usually σ for stress relaxation and ϵ for creep
- `_sampling`: Declare whether sample rate is `constant` or `variable` so that convolution or integration is used respectively
- `_insight`: Declare whether insight info should be shown when this function is called, true or false
"""
function objectivefunc(params::Array{Float64,1}, grad::Array{Float64,1},
                            model, time_series, dt_series, prescribed_dot,
                            measured; _insight = false, _sampling = "constant")

    # convolution if constant, brute force integration if variable
    if _sampling == "constant"
        convolved = boltzconvolve(model, time_series, dt_series, params, prescribed_dot)
    elseif _sampling == "variable"
        convolved = boltzintegral(model, time_series, params, prescribed_dot)
    end

    # don't use first element if singularity exists in model
    if model.singularity
        cost = sum(0.5*(measured[2:end] - convolved).^2)
    else
        cost = sum(0.5*(measured - convolved).^2)
    end

    # if insight=true then output parameters used, add auto-updating plot?
    if _insight
        println("Current Parameters: ", params)
    end

    cost
end

"""
    leastsquares_init(params_init, low_bounds, hi_bounds, model, time_series, dt_series, prescribed_dot, measured; sampling = "constant", insight = false)

Initialise then begin a least squares fitting of the supplied data.

# Arguments

- `params_init`: Initial parameters to be used (starting guess)
- `low_bounds`: Lower bounds for parameters
- `hi_bounds`: Higher bounds for parameters
- `model`: RheologyModel struct containing viscoelastic modulus function
- `time_series`: Array of time data
- `dt_series`: Array of time data differences,  equal to deriv(time_series)
- `prescribed_dot`: Convolved with `model` over `time_series`, usually dϵ/dt for stress relaxation and dσ/dt for creep
- `measured`: Data for comparison against, usually σ for stress relaxation and ϵ for creep
- `sampling`: Declare whether sample rate is `constant` or `variable` so that convolution or integration is used respectively
- `insight`: Declare whether insight info should be shown when this function is called, true or false
"""
function leastsquares_init(params_init, low_bounds, hi_bounds,
                            model, time_series, dt_series, prescribed_dot,
                            measured; sampling = "constant",
                            insight = false)

    # initialise NLOpt.Opt object with :LN_SBPLX Subplex algorithm
    opt = Opt(:LN_SBPLX, length(params_init))
    # set lower bounds and upper bounds
    lower_bounds!(opt, low_bounds)
    upper_bounds!(opt, hi_bounds)
    # set relative tolerance
    xtol_rel!(opt,1e-4)

    # set Opt object as a minimisation objective. Use a closure for additional
    # arguments sent to object objectivefunc
    min_objective!(opt, (params, grad) -> objectivefunc(params, grad, model,
                                                        time_series, dt_series,
                                                        prescribed_dot, measured;
                                                        _insight = insight,
                                                        _sampling = sampling))

    # minimise objective func, minx are the parameters resulting in minimum
    (minf,minx,ret) = optimize(opt, params_init)

    # return all
    (minf,minx,ret)
end
