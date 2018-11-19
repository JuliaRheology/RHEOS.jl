#!/usr/bin/env julia

#############################
#~ Preprocessing Functions ~#
#############################

"""
    var_resample(self::RheologyData, refvar::Symbol, pcntdownsample::Float64; mapback::Bool = false)

Convert a fixed sample rate array to a variable sample rate, with sampling points
added according to a relative change in chosen variable `refvar`, 1st derivative
of `refvar` and 2nd derivative of `refvar` (WRT time). Usually chosen as the
measured variable, so :σ for a stress relaxation test and :ϵ for a creep test.

Currently only variable downsampling supported. pcntdown sample is approximate,
works well in some cases and very poorly in others. If required, compare resampled
length vs original length after processing has finished. If data is noisy, may
benefit from sending smoothed signal to this algorithm and either using mapback
function or interpolating onto unsmoothed data.

See help docstring for `var_resample` for more details on algorithm implementation.
"""
function var_resample(self::RheologyData, refvar::Symbol, pcntdownsample::Float64; _mapback::Bool = false)

    # enforce minimum period of original period/10
    _minperiod = (self.t[2] - self.t[1])/10.0

    # get time resampled with respect to refvar
    (tᵦ, dummy)  = var_resample(self.t, getfield(self, refvar), pcntdownsample, _minperiod)

    local σ::Array{Float64,1}
    local ϵ::Array{Float64,1}
    local t::Array{Float64,1}

    if _mapback
        # get mapped indices wrt original
        mapback_indices = mapback(tᵦ, self.t)

        σ = self.σ[mapback_indices]
        ϵ = self.ϵ[mapback_indices]
        t = self.t[mapback_indices]
        
    elseif !_mapback
        # interpolate with respect to t
        σ_interp = Interpolations.interpolate((self.t,), self.σ, Interpolations.Gridded(Interpolations.Linear()))
        ϵ_interp = Interpolations.interpolate((self.t,), self.ϵ, Interpolations.Gridded(Interpolations.Linear()))
        t_interp = Interpolations.interpolate((self.t,), self.t, Interpolations.Gridded(Interpolations.Linear()))

        # resample all using new timepoints tᵦ
        σ = σ_interp[tᵦ]
        ϵ = ϵ_interp[tᵦ]
        t = t_interp[tᵦ]
    end

    # change sampling type to variable
    sampling = "variable"
    
    # add record of operation applied
    log = vcat(self.log, "var_resample - refvar: $refvar, pcntdownsample: $pcntdownsample, mapback: $_mapback")

    self_new = RheologyData(σ, ϵ, t, sampling, log)

end

"""
    downsample(self::RheologyData, time_boundaries::Array{Float64,1}, elperiods::Array{Int64,1})

High-level RheologyData interface to downsample in base.jl. Boundaries are floating point times which are then converted to the closest elements.
"""
function downsample(self::RheologyData, time_boundaries::Array{Float64,1}, elperiods::Array{Int64,1})

    # convert boundaries from times to element indicies
    boundaries = closestindices(self.t, time_boundaries)

    # get downsampled indices
    indices = downsample(boundaries::Array{Int64,1}, elperiods::Array{Int64,1})

    # downsample data
    σ = self.σ[indices]
    ϵ = self.ϵ[indices]
    t = self.t[indices]

    # change to variable sampling rate if more than one section, if not then as original
    local sampling::String
    if length(elperiods) > 1
        sampling = "variable"
    else
        sampling = self.sampling
    end

    # add record of operation applied
    log = vcat(self.log, "downsample - boundaries: $boundaries, elperiods: $elperiods")

    self_new = RheologyData(σ, ϵ, t, sampling, log)

end

"""
    fixed_resample(self::RheologyData, time_boundaries::Array{Float64,1}, elperiods::Array{Int64,1}, direction::Array{String,1})

High-level RheologyData interface to fixed_resample in base.jl
"""
function fixed_resample(self::RheologyData, time_boundaries::Array{Float64,1}, elperiods::Array{Int64,1}, direction::Array{String,1})

    # convert boundaries from times to element indicies
    boundaries = closestindices(self.t, time_boundaries)

    # resample all data
    (t, σ) = fixed_resample(self.t, self.σ, boundaries, elperiods, direction)
    (t, ϵ) = fixed_resample(self.t, self.ϵ, boundaries, elperiods, direction)

    # change to variable sampling rate if more than one section
    local sampling::String
    if length(elperiods) > 1
        sampling = "variable"
    else
        sampling = self.sampling
    end

    # add record of operation applied
    log = vcat(self.log, "fixed_resample - boundaries: $boundaries, elperiods: $elperiods, direction: $direction")

    self_new = RheologyData(σ, ϵ, t, sampling, log)

end

"""
    smooth(self::RheologyData, τ::Float64)

Smooth data using a Gaussian Kernel to time scale τ (approximately half power).

Smooths both σ and ϵ.
"""
function smooth(self::RheologyData, τ::Float64)

    @assert self.sampling=="constant" "Sample rate must be constant for Gaussian smoothing kernel to function properly"

    t = self.t

    # get constant sample period
    dt = t[2] - t[1]

    # sample rate
    samplerate = 1.0/dt;

    σ = smoothgauss(self.σ, τ, samplerate)
    ϵ = smoothgauss(self.ϵ, τ, samplerate)

    # add record of operation applied
    log = vcat(self.log, "smooth - τ: $τ")

    self_new = RheologyData(σ, ϵ, t, self.sampling, log)

end

"""
    function mapbackdata(self_new::RheologyData, self_original::RheologyData)

Map back elements (WRT closest time elements) of all data from self_new to
self_original. See mapback help docstring for more info on how algorithm works.
"""
function mapbackdata(self_new::RheologyData, self_original::RheologyData)

    # get mapped back indices
    indices = mapback(self_new.t, self_original.t)

    #  set variables
    σ = self_original.σ[indices]
    ϵ = self_original.ϵ[indices]
    t = self_original.t[indices]

    # add record of operation applied
    log = vcat(self_new.log, "mapped back") # mapped back to what? Include self_origina.log as well? Need to change Array type if so.

    # to add, check if sample rate is now variable or constant (unlikely but could fall back to constant?)

    self_new = RheologyData(σ, ϵ, t, self_new.sampling, log)

end

function zerotime(self::RheologyData)

    return RheologyData(self.σ, self.ϵ, self.t .- minimum(self.t), self.sampling, vcat(self.log, ["Normalized time to start at 0.0"]))

end

##########################
#~ Processing Functions ~#
##########################

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
                  model::RheologyModel,
                  modtouse::Symbol;
                  p0::Array{Float64,1} = [-1.0],
                  lo::Array{Float64,1} = [-1.0],
                  hi::Array{Float64,1} = [-1.0],
                  verbose::Bool = false,
                  rel_tol = 1e-4,
                  diff_method="BD")::RheologyModel

    # get modulus function
    modulus = getfield(model, modtouse)

    # use default p0 if not provided
    if quasinull(p0)
        p0 = model.parameters
    end 

    # get singularity presence
    sing = singularitytest(modulus, p0)

    # get time step (only needed for convolution, which requires constant so t[2]-t[1] is sufficient)
    dt = data.t[2] - data.t[1]

    # TEMP - CHECK WITH ALE AND ALEXANDRE BUT IS DEFINITELY NECESSARY
    # time must start at 0 for convolution to work properly!
    t_zeroed = data.t .- minimum(data.t)

    # use correct method for derivative
    if diff_method=="BD"
        deriv = derivBD
    elseif diff_method=="CD"
        deriv = derivCD
    end

    # get derivative of controlled variable and measured variable
    if modtouse == :J
        dcontrolled = deriv(data.σ, data.t)
        measured = data.ϵ
    elseif modtouse == :G
        dcontrolled = deriv(data.ϵ, data.t)
        measured = data.σ
    end

    # fit
    (minf, minx, ret), timetaken, bytes, gctime, memalloc = @timed leastsquares_init(p0,
                                                                                    lo, 
                                                                                    hi,
                                                                                    modulus, 
                                                                                    t_zeroed, 
                                                                                    dt, 
                                                                                    dcontrolled,
                                                                                    measured; 
                                                                                    insight = verbose,
                                                                                    sampling = data.sampling, 
                                                                                    singularity = sing,
                                                                                    _rel_tol = rel_tol)

    modulusname = string(modulus)

    log = vcat(data.log, "Fitted $modulusname, Time: $timetaken s, Why: $ret, Parameters: $minx, Error: $minf")

    RheologyModel(model.G, model.J, model.Gp, model.Gpp, minx, log)

end

"""
    modelpredict(data::RheologyData, model::RheologyModel)

Given data and model and parameters, predict new dataset based on both.
"""
function modelpredict(data::RheologyData, model::RheologyModel, modtouse::Symbol; diff_method="BD")::RheologyData

    # get modulus
    modulus = getfield(model, modtouse)

    # get singularity presence
    sing = singularitytest(modulus, model.parameters)

    # use correct method for derivative
    if diff_method=="BD"
        deriv = derivBD
    elseif diff_method=="CD"
        deriv = derivCD
    end

    if modtouse == :J
        dcontrolled = deriv(data.σ, data.t)
    elseif modtouse == :G
        dcontrolled = deriv(data.ϵ, data.t)
    end

    # get time step (only needed for convolution, which requires constant so t[2]-t[1] is sufficient)
    dt = data.t[2] - data.t[1]

    # TEMP - CHECK WITH ALE AND ALEXANDRE BUT IS DEFINITELY NECESSARY
    # time must start at 0 for convolution to work properly!
    t_zeroed = data.t .- minimum(data.t)

    # get convolution
    if !sing && data.sampling == "constant"
        convolved = boltzconvolve_nonsing(modulus, t_zeroed, dt, model.parameters, dcontrolled)

    elseif sing && data.sampling == "constant"
        # convolved = boltzconvolve_sing(modulus, t_zeroed, dt, model.parameters, dcontrolled)
        t_zeroed[1] = 0.0 + (t_zeroed[2] - t_zeroed[1])/10.0
        # convolved = boltzconvolve_sing(modulus, t_zeroed, dt, model.parameters, dcontrolled)
        convolved = boltzconvolve(modulus, t_zeroed, dt, model.parameters, dcontrolled)

    elseif !sing && data.sampling == "variable"
        convolved = boltzintegral_nonsing(modulus, t_zeroed, model.parameters, dcontrolled)

    elseif sing && data.sampling == "variable"
        convolved = boltzintegral_sing(modulus, t_zeroed, model.parameters, dcontrolled)

    end

    # code commented out has off by one error
    # need to find a better way to handle singularities
    # if sing
    #     if modtouse == :J
    #         σ = data.σ[2:end]
    #         ϵ = convolved
    #     elseif modtouse == :G
    #         σ = convolved
    #         ϵ = data.ϵ[2:end]
    #     end   
    #     t = data.t[2:end]

    if sing
        if modtouse == :J
            σ = data.σ
            ϵ = convolved
        elseif modtouse == :G
            σ = convolved
            ϵ = data.ϵ
        end 
        t = data.t

    elseif !sing
        if modtouse == :J
            σ = data.σ
            ϵ = convolved
        elseif modtouse == :G
            σ = convolved
            ϵ = data.ϵ
        end 
        t = data.t

    end
        
    # store operation
    log = vcat(data.log, "Predicted data from model:", model.log)

    RheologyData(σ, ϵ, t, data.sampling, log)

end

function modelstepfit(data::RheologyData,
                  model::RheologyModel,
                  modtouse::Symbol;
                  p0::Array{Float64,1} = [-1.0],
                  lo::Array{Float64,1} = [-1.0],
                  hi::Array{Float64,1} = [-1.0],
                  verbose::Bool = false,
                  rel_tol = 1e-4)::RheologyModel

    # get modulus function
    modulus = getfield(model, modtouse)

    # use default p0 if not provided
    if quasinull(p0)
        p0 = model.parameters
    end 

    # get singularity presence
    sing = singularitytest(modulus, p0; t1 = data.t[1])
    
    # get controlled variable (just take first element as it's a step) and measured variable
    if modtouse == :J
        controlled = data.σ[1]
        measured = data.ϵ
    elseif modtouse == :G
        controlled = data.ϵ[1]
        measured = data.σ
    end

    # start fit
    (minf, minx, ret), timetaken, bytes, gctime, memalloc = @timed leastsquares_stepinit(p0,
                                                                                        lo, 
                                                                                        hi,
                                                                                        modulus, 
                                                                                        data.t, 
                                                                                        controlled,
                                                                                        measured; 
                                                                                        insight = verbose,
                                                                                        singularity = sing,
                                                                                        _rel_tol = rel_tol)

    modulusname = string(modulus)

    log = vcat(data.log, "Fitted $modulusname, Time: $timetaken s, Why: $ret, Parameters: $minx, Error: $minf")

    RheologyModel(model.G, model.J, model.Gp, model.Gpp, minx, log)

end

function modelsteppredict(data::RheologyData, model::RheologyModel, modtouse::Symbol; step_on::Real = 0.0)::RheologyData

    # get modulus
    modulus = getfield(model, modtouse)

    # check singularity presence at time closest to step
    stepon_el = closestindex(data.t, step_on)
    sing = singularitytest(modulus, model.parameters; t1 = (data.t[stepon_el] - step_on))

    if modtouse == :J
        controlled = data.σ[1]
    elseif modtouse == :G
        controlled = data.ϵ[1]
    end

    # get predicted
    if !sing
        predicted = zeros(length(data.t))
        predicted[stepon_el:end] = controlled*modulus(data.t[stepon_el:end] .- step_on, model.parameters)

    elseif sing
        predicted = zeros(length(data.t))
        predicted[(stepon_el + 1):end] = controlled*modulus(data.t[(stepon_el + 1):end] .- step_on, model.parameters)

    end

    if modtouse == :J
        σ = data.σ
        ϵ = predicted
    elseif modtouse == :G
        σ = predicted
        ϵ = data.ϵ
    end 
        
    # store operation
    log = vcat(data.log, "Predicted data from model:", model.log)

    RheologyData(σ, ϵ, data.t, data.sampling, log)

end

function obj_dynamic(params::Vector{T},
                     grad::Vector{T},
                     ω::Vector{T},
                     dataGp::Vector{T},
                     dataGpp::Vector{T},
                     modelGp::Function,
                     modelGpp::Function;
                     _insight::Bool = false) where T<:Real

    if _insight
        println("Current Parameters: ", params)
    end

    costGp = sum(0.5*(dataGp - modelGp(ω, params)).^2)
    costGpp = sum(0.5*(dataGpp - modelGpp(ω, params)).^2)

    cost = costGp + costGpp

end

function obj_dynamic_linear(params::Vector{T},
                            grad::Vector{T},
                            ω::Vector{T},
                            dataGp::Vector{T},
                            dataGpp::Vector{T},
                            modelGp::Function,
                            modelGpp::Function,
                            meanGp::T,
                            meanGpp::T;
                            _insight::Bool = false) where T<:Real

    if _insight
        println("Current Parameters: ", params)
    end

    costGp = sum(0.5*(dataGp/meanGp - modelGp(ω, params)/meanGp).^2)
    costGpp = sum(0.5*(dataGpp/meanGpp - modelGpp(ω, params)/meanGpp).^2)

    cost = costGp + costGpp

end

function obj_dynamic_log(params::Vector{T},
                     grad::Vector{T},
                     ω::Vector{T},
                     dataGp::Vector{T},
                     dataGpp::Vector{T},
                     modelGp::Function,
                     modelGpp::Function;
                     _insight::Bool = false) where T<:Real

    if _insight
        println("Current Parameters: ", params)
    end

    costGp = sum(0.5*(log.(dataGp) - log.(modelGp(ω, params))).^2)
    costGpp = sum(0.5*(log.(dataGpp) - log.(modelGpp(ω, params))).^2)

    cost = costGp + costGpp

end

function obj_dynamic_global(params::Vector{T},
                     grad::Vector{T},
                     ω::Vector{T},
                     dataGp::Vector{T},
                     dataGpp::Vector{T},
                     modelGp::Function,
                     modelGpp::Function;
                     _insight::Bool = false) where T<:Real

    if _insight
        println("Current Parameters: ", params)
    end

    costGp = sum(0.5*(((dataGp - modelGp(ω, params))./dataGp).^2))
    costGpp = sum(0.5*(((dataGpp - modelGpp(ω, params))./dataGpp).^2))

    cost = costGp + costGpp

end

function obj_dynamic_manual(params::Vector{T},
                            grad::Vector{T},
                            ω::Vector{T},
                            dataGp::Vector{T},
                            dataGpp::Vector{T},
                            modelGp::Function,
                            modelGpp::Function,
                            weights::Vector{T};
                            _insight::Bool = false) where T<:Real

    if _insight
        println("Current Parameters: ", params)
    end

    costGp = sum(0.5*(dataGp - modelGp(ω, params)).^2)
    costGpp = sum(0.5*(dataGpp - modelGpp(ω, params)).^2)

    cost = weights[1]*costGp + weights[2]*costGpp
end

function dynamicmodelfit(data::RheologyDynamic,
                model::RheologyModel;
                p0::Vector{T} = [-1.0],
                lo::Vector{T} = [-1.0],
                hi::Vector{T} = [-1.0],
                verbose::Bool = false,
                rel_tol::T = 1e-4,
                weights::Union{String, Vector{T}}="log") where T<:Real

    # get initial paramaters
    if quasinull(p0)
        p0 = model.parameters
    end 

    # initialise NLOpt.Opt object with :LN_SBPLX Subplex algorithm
    opt = Opt(:LN_SBPLX, length(p0))

    # apply parameter boundaries if prescribed
    if !quasinull(lo)
        lower_bounds!(opt, lo)
    end

    if !quasinull(hi)
        upper_bounds!(opt, hi)
    end

    # set relative tolerance
    xtol_rel!(opt, rel_tol)

    # set objective/cost function
    if weights=="none"
        min_objective!(opt, (params, grad) -> obj_dynamic(params, grad, data.ω, data.Gp, data.Gpp, model.Gp, model.Gpp; _insight = verbose))

    elseif weights=="linear"
        min_objective!(opt, (params, grad) -> obj_dynamic_linear(params, grad, data.ω, data.Gp, data.Gpp, model.Gp, model.Gpp, mean(data.Gp), mean(data.Gpp); _insight = verbose))

    elseif weights=="log"
        min_objective!(opt, (params, grad) -> obj_dynamic_log(params, grad, data.ω, data.Gp, data.Gpp, model.Gp, model.Gpp; _insight = verbose))

    elseif weights=="global"
        min_objective!(opt, (params, grad) -> obj_dynamic_global(params, grad, data.ω, data.Gp, data.Gpp, model.Gp, model.Gpp; _insight = verbose))

    elseif typeof(weights)==Vector{T} && length(weights)==2
        min_objective!(opt, (params, grad) -> obj_dynamic_manual(params, grad, data.ω, data.Gp, data.Gpp, model.Gp, model.Gpp, weights; _insight = verbose))

    end

    # timed fitting
    # (minf, minx, ret), timetaken, bytes, gctime, memalloc = @timed NLopt.optimize(opt, p0)
    (minf, minx, ret) = NLopt.optimize(opt, p0)
    timetaken = -1.0

    println(ret)

    # log fit details
    modelname = string(model)
    log = vcat(data.log, "Fitted Gp, Gpp of $modelname, Time: $timetaken s, Why: $ret, Parameters: $minx, Error: $minf")

    RheologyModel(model.G, model.J, model.Gp, model.Gpp, minx, log)

end

function dynamicmodelpredict(data::RheologyDynamic, model::RheologyModel)

    # get results
    predGp = model.Gp(data.ω, model.parameters)
    predGpp = model.Gpp(data.ω, model.parameters)

    # store operation
    log = vcat(data.log, "Predicted data from model:", model.log)

    RheologyDynamic(predGp, predGpp, data.ω, log)

end