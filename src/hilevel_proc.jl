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
                  hi_bounds::Array{Float64,1};)

    # initialise NLOpt.Opt object with :LN_SBPLX Subplex algorithm
    opt = Opt(:LN_SBPLX, length(params_init))
    # set lower bounds and upper bounds
    lower_bounds!(opt, low_bounds)
    upper_bounds!(opt, hi_bounds)
    # set relative tolerance
    xtol_rel!(opt,1e-4)

    # generate time series difference array (for convolution)
    dt_series = deriv(self.t)

    # measured and prescribed_dot depend on test type
    if self.test_type == "strlx"
        measured = self.σ
        prescribed_dot = self.dϵ
    elseif self.test_type == "creep"
        measured = self.ϵ
        prescribed_dot = self.dσ
    end

    # get modulus function
    modulus = moduli(model, self.test_type)

    # set Opt object as a minimisation objective. Use a closure for additional
    # arguments sent to object objectivefunc
    min_objective!(opt, (params, grad) -> objectivefunc(params, grad, modulus,
                                                        self.t, dt_series,
                                                        prescribed_dot, measured;
                                                        _insight = self.insight,
                                                        _sampling = self.sampling))

    # minimise objective func, minx are the parameters resulting in minimum
    (minf, minx, ret) = optimize(opt, params_init)

    # store fit results in RheologyData struct's fittedmodels dictionary
    self.fittedmodels[model] = minx
end

"""
    modelcomplete!(self::RheologyData, model::String, params::Array{Float64,1})

Given partial data (just t and ϵ for "strlx" test or t and σ for "creep"),
model and parameters, find missing data (σ for "strlx" and ϵ for "creep").
"""
function modelcomplete!(self::RheologyData, model::String, params::Array{Float64,1})

    # generate time series difference array (for convolution)
    dt_series = deriv(self.t)

    # measured and prescribed_dot depend on test type
    if self.test_type == "strlx"
        prescribed_dot = self.dϵ
    elseif self.test_type == "creep"
        prescribed_dot = self.dσ
    end

    # get modulus function
    modulus = moduli(model, self.test_type)

    # get convolution
    if self.sampling == "constant"
        convolved = boltzconvolve(modulus, self.t, dt_series, params, prescribed_dot)
    elseif self.sampling == "variable"
        convolved = boltzintegral(modulus, self.t, params, prescribed_dot)
    end

    # store output in RheologyData object
    if self.test_type == "strlx"
        self.σ = convolved
    elseif self.test_type == "creep"
        self.ϵ = convolved
    end
end
#= think about extending what is stored in fittedmodels dict, could be
params, convolved with those params, cost and more. Then modelcomplete!
should place this convolved there, not in main self =#
