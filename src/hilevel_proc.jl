#!/usr/bin/env julia

"""
    modelfit!(self::RheologyData, test_type::String, model::String, params_init::Array{Float64,1}, low_bounds::Array{Float64,1}, hi_bounds::Array{Float64,1}; singularity = false)

Fit RheologyData struct to model and store fitted parameters in self.fittedmodels.

# Arguments

- `self`: RheologyData struct containing all data
- `test_type`: Either "creep" or "strlx"
- `model`: E.g. "SLS", "springpot", "burgers" etc. See models.jl for full list
- `params_init`: Initial parameters to use in fit
- `low_bounds`: Lower bounds for parameters
- `hi_bounds`: Higher bounds for parameters
- `singularity`: Declare whether there is a singularity in viscoelastic modulus at t=0, true or false
"""
function modelfit!(self::RheologyData,
                  test_type::String,
                  model::String,
                  params_init::Array{Float64,1},
                  low_bounds::Array{Float64,1},
                  hi_bounds::Array{Float64,1};
                  singularity = false)

    # initialise NLOpt.Opt object with :LN_SBPLX Subplex algorithm
    opt = Opt(:LN_SBPLX, length(params_init))
    # set lower bounds and upper bounds
    lower_bounds!(opt, low_bounds)
    upper_bounds!(opt, hi_bounds)
    # set relative tolerance
    xtol_rel!(opt,1e-4)

    # generate time series difference array (for convolution)
    dt_series = deriv(self.tᵦ)

    # measured and prescribed_dot depend on test type
    if test_type == "strlx"
        measured = self.σᵦ
        prescribed_dot = self.dϵᵦ
    elseif test_type == "creep"
        measured = self.ϵᵦ
        prescribed_dot = self.dσᵦ
    end

    # get modulus function
    modulus = moduli(model, test_type)

    # set Opt object as a minimisation objective. Use a closure for additional
    # arguments sent to object objectivefunc
    min_objective!(opt, (params, grad) -> objectivefunc(params, grad, modulus,
                                                        self.tᵦ, dt_series,
                                                        prescribed_dot, measured;
                                                        _singularity = singularity,
                                                        _insight = self.insight,
                                                        _sampling = self.sampling))

    # minimise objective func, minx are the parameters resulting in minimum
    (minf, minx, ret) = optimize(opt, params_init)

    # store fit results in RheologyData struct
    self.fittedmodels[model] = minx
end
