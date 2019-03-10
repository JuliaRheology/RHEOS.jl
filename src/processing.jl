#!/usr/bin/env julia

#############################
#~ Preprocessing Functions ~#
#############################

export cutting

"""
    fixedresample(self::RheoTimeData, elperiods::Union{Vector{K},K}; time_boundaries::Vector{T}= [-1])

Resample data with new sample rate(s).

Fixedresample can downsample or upsample data. If the number of elperiods is negative it is going to reduce the number of samples,
viceversa if it is positive. If time boundaries are not specified, resampling is applied to the whole set of data.
"""
function fixedresample(self::RheoTimeData, elperiods::Union{Vector{K},K}; time_boundaries::Vector{T}= [-1]) where {K<:Integer,T<:Real}

    @assert (count(iszero,elperiods)==0) "Number of elements cannot be zero"
    # convert boundaries from times to element indicies
    if time_boundaries ==[-1]
        boundaries = [1,length(self.t)];
    else
        boundaries = closestindices(self.t, time_boundaries)
    end

    check = RheoTimeDataType(self)
    if (Int(check) == 1)
        (time, epsilon) = fixed_resample(self.t, self.ϵ, boundaries, elperiods)
        sigma = [];
    elseif (Int(check) == 2)
        (time, sigma) = fixed_resample(self.t, self.σ, boundaries, elperiods)
        epsilon = [];
    elseif (Int(check) == 3)
        (time, sigma) = fixed_resample(self.t, self.σ, boundaries, elperiods)
        (time, epsilon) = fixed_resample(self.t, self.ϵ, boundaries, elperiods)
    end

    # add record of operation applied
    log = vcat(self.log, "fixed_resample - boundaries: $boundaries, elperiods: $elperiods")

    self_new = RheoTimeData(sigma, epsilon, time, log)

end

function cutting(self::RheoTimeData, time_on::T1,time_off::T2) where {T1<:Number, T2<:Number}

    @assert isreal(time_on) && isreal(time_off) "Boundaries cannot be complex numbers"
    
    boundary_on = closestindex(self.t, time_on);
    boundary_off = closestindex(self.t, time_off);
    time = self.t[boundary_on:boundary_off]

    check = RheoTimeDataType(self)
    if (Int(check) == 1)
        epsilon = self.ϵ[boundary_on:boundary_off]
        sigma = [];
    elseif (Int(check) == 2)
        sigma = self.σ[boundary_on:boundary_off]
        epsilon = [];
    elseif (Int(check) == 3)
        epsilon = self.ϵ[boundary_on:boundary_off]
        sigma = self.σ[boundary_on:boundary_off]
    end
    log = vcat(self.log, "Data from $time_on to $time_off extracted.")

    return RheoTimeData(sigma,epsilon,time,log)

end

"""
    smooth(self::RheologyData, τ::Real; pad::String="replicate")

Smooth data using a Gaussian Kernel to time scale τ (approximately half power).

Smooths both σ and ϵ. Essentially a low pass filter with frequencies of 1/τ being cut to approximately
half power. For other pad types available see ImageFiltering documentation.
"""
function smooth(self::RheoTimeData, τ::Real; pad::String="reflect")

    @eval import ImageFiltering: imfilter, Kernel

    # get sample-rate and Gaussian kernel (std. dev)
    samplerate = 1.0/getsampleperiod(self.t)
    Σ = getsigma(τ, samplerate)

    # smooth signal and return
    check = RheoTimeDataType(self)
    if (Int(check) == 1)
        epsilon = Base.invokelatest(imfilter, self.ϵ, Base.invokelatest(Kernel.reflect, Base.invokelatest(Kernel.gaussian, (Σ,))), pad)
        sigma = [];
    elseif (Int(check) == 2)
        sigma = Base.invokelatest(imfilter, self.σ, Base.invokelatest(Kernel.reflect, Base.invokelatest(Kernel.gaussian, (Σ,))), pad)
        epsilon = [];
    elseif (Int(check) == 3)
        sigma = Base.invokelatest(imfilter, self.σ, Base.invokelatest(Kernel.reflect, Base.invokelatest(Kernel.gaussian, (Σ,))), pad)
        epsilon = Base.invokelatest(imfilter, self.ϵ, Base.invokelatest(Kernel.reflect, Base.invokelatest(Kernel.gaussian, (Σ,))), pad)
    end

    # add record of operation applied
    log = vcat(self.log, "smooth - τ: $τ")

    self_new = RheoTimeData(sigma, epsilon, self.t, log)

end

##########################
#~ Processing Functions ~#
##########################

"""
    modelfit(data::RheologyData, model::RheologyModel, modtouse::Symbol; p0::Vector{T} = [-1.0], lo::Vector{T} = [-1.0], hi::Vector{T} = [-1.0], verbose::Bool = false, rel_tol = 1e-4, diff_method="BD") where T<:Real

Fit RheologyData struct to model and return a fitted model as a RheologyModel object.

# Arguments

- `data`: RheologyData struct containing all data
- `model`: RheologyModel containing moduli and default (initial) parameters
- `modtouse`: :G for relaxation modulus, :J for creep modulus
- `p0`: Initial parameters to use in fit (uses 'model' parameters if none given)
- `lo`: Lower bounds for parameters
- `hi`: Upper bounds for parameters
- `verbose`: If true, prints parameters on each optimisation iteration
- `rel_tol`: Relative tolerance of optimization, see NLOpt docs for more details
- `diff_method`: Set finite difference formula to use for derivative, currently "BD" or "CD"
"""
function modelfit(data::RheoTimeData,
                  model::RheologyModel,
                  modloading::Union{LoadingType,Integer};
                  p0::Array{T1,1} = [-1.0],
                  lo::Array{T2,1} = [-1.0],
                  hi::Array{T3,1} = [-1.0],
                  verbose::Bool = false,
                  rel_tol = 1e-4,
                  diff_method="BD") where {T1<:Real, T2<:Real, T3<:Real}

    p0 = convert(Array{RheoFloat,1},p0)
    lo = convert(Array{RheoFloat,1},lo)
    hi = convert(Array{RheoFloat,1},hi)
    rel_tol = convert(RheoFloat,rel_tol)

    check = RheoTimeDataType(data)
    @assert (Int(check) == 3) "Both stress and strain are required"

    # use correct method for derivative
    if diff_method=="BD"
        deriv = derivBD
    elseif diff_method=="CD"
        deriv = derivCD
    end

    # get modulus function and derivative
    if Int(modloading) == 2
        dcontrolled = deriv(data.σ, data.t)
        measured = data.ϵ
        modtouse = :J
    elseif Int(modloading) == 1
        dcontrolled = deriv(data.ϵ, data.t)
        measured = data.σ
        modtouse = :G
    end

    modulus = getfield(model, modtouse)

    # use default p0 if not provided
    if quasinull(p0)
        p0 = convert(Array{RheoFloat,1},model.parameters)
    end

    # get singularity presence
    sing = singularitytest(modulus, p0)

    # get time step (only needed for convolution, which requires constant dt so t[2]-t[1] is sufficient)
    dt = data.t[2] - data.t[1]

    # TEMP - CHECK WITH ALE AND ALEXANDRE BUT IS DEFINITELY NECESSARY
    # time must start at 0 for convolution to work properly!
    t_zeroed = data.t .- minimum(data.t)


    # fit
    sampling_check = constantcheck(data.t)

    (minf, minx, ret), timetaken, bytes, gctime, memalloc = @timed leastsquares_init(p0,
                                                                                    lo,
                                                                                    hi,
                                                                                    modulus,
                                                                                    t_zeroed,
                                                                                    dt,
                                                                                    dcontrolled,
                                                                                    measured;
                                                                                    insight = verbose,
                                                                                    sampling = sampling_check,
                                                                                    singularity = sing,
                                                                                    _rel_tol = rel_tol)

    modulusname = string(modulus)

    log = vcat(data.log, "Fitted $modulusname, Modulus used: $modtouse, Time: $timetaken s, Why: $ret, Parameters: $minx, Error: $minf")

    RheologyModel(model.G, model.J, model.Gp, model.Gpp, minx, log)

end

"""
    modelpredict(data::RheologyData, model::RheologyModel, modtouse::Symbol; diff_method="BD")

Given data and model, return new dataset based on model parameters and using the
modulus specified by 'modtouse'; either creep modulus (:J, only returned strain is new) or
relaxation modulus (:G, only returned stress is new). 'diff_method' sets finite difference for
calculating the derivative used in the hereditary integral and can be either backwards difference
("BD") or central difference ("CD").
"""
function modelpredict(data::RheoTimeData,model::RheologyModel; diff_method="BD")

    # use correct method for derivative
    if diff_method=="BD"
        deriv = derivBD
    elseif diff_method=="CD"
        deriv = derivCD
    end

    check = RheoTimeDataType(data)
    @assert (Int(check) ==1)||(Int(check) ==2) "Both stress and strain are already defined"

    if (Int(check) == 1)
        modtouse = :G;
        dcontrolled = deriv(data.ϵ, data.t)
    elseif (Int(check) == 2)
        modtouse = :J;
        dcontrolled = deriv(data.σ, data.t)
    end

    # get modulus
    modulus = getfield(model, modtouse)

    # get singularity presence
    sing = singularitytest(modulus, model.parameters)

    # get time step (only needed for convolution, which requires constant so t[2]-t[1] is sufficient)
    dt = data.t[2] - data.t[1]

    # TEMP - CHECK WITH ALE AND ALEXANDRE BUT IS DEFINITELY NECESSARY
    # time must start at 0 for convolution to work properly!
    t_zeroed = data.t .- minimum(data.t)

    # get convolution
    if !sing && constantcheck(data.t)
        convolved = boltzconvolve_nonsing(modulus, t_zeroed, dt, model.parameters, dcontrolled)

    elseif sing && constantcheck(data.t)
        # convolved = boltzconvolve_sing(modulus, t_zeroed, dt, model.parameters, dcontrolled)
        t_zeroed[1] = 0.0 + (t_zeroed[2] - t_zeroed[1])/10.0
        # convolved = boltzconvolve_sing(modulus, t_zeroed, dt, model.parameters, dcontrolled)
        convolved = boltzconvolve(modulus, t_zeroed, dt, model.parameters, dcontrolled)

    elseif !sing && !constantcheck(data.t)
        convolved = boltzintegral_nonsing(modulus, t_zeroed, model.parameters, dcontrolled)

    elseif sing && !constantcheck(data.t)
        convolved = boltzintegral_sing(modulus, t_zeroed, model.parameters, dcontrolled)

    end

    if modtouse == :J
        sigma = data.σ
        epsilon = convolved
        pred_mod = model.J
    elseif modtouse == :G
        sigma = convolved
        epsilon = data.ϵ
        pred_mod = model.G
    end
    time = data.t

    modparam = model.parameters
    # store operation
    log = vcat( data.log, "Predicted data using: $pred_mod, Parameters: $modparam")


    return RheoTimeData(sigma,epsilon,time, log)

end

"""
    modelstepfit(data::RheologyData, model::RheologyModel, modtouse::Symbol; p0::Vector{T} = [-1.0], lo::Vector{T} = [-1.0], hi::Vector{T} = [-1.0], verbose::Bool = false, rel_tol = 1e-4) where T<:Real

Same as 'modelfit' except assumes a step loading. If this assumption is appropriate for the data
then fitting can be sped up greatly by use of this function. If modtouse is :G, relaxation modulus,
then the first element of the strain is assumed to be the amplitude of the step. If modtouse is :j,
creep modulus, then the first element of the stress is assumed to be the amplitude of the step.

# Arguments

- `data`: RheologyData struct containing all data
- `model`: RheologyModel containing moduli and default (initial) parameters
- `modtouse`: :G for relaxation modulus, :J for creep modulus
- `p0`: Initial parameters to use in fit (uses 'model' parameters if none given)
- `lo`: Lower bounds for parameters
- `hi`: Upper bounds for parameters
- `verbose`: If true, prints parameters on each optimisation iteration
- `rel_tol`: Relative tolerance of optimization, see NLOpt docs for more details
"""
function modelstepfit(data::RheoTimeData,
                  model::RheologyModel,
                  modloading::Union{LoadingType,Integer};
                  step = nothing,
                  p0::Array{T1,1} = [-1.0],
                  lo::Array{T2,1} = [-1.0],
                  hi::Array{T3,1} = [-1.0],
                  verbose::Bool = false,
                  rel_tol = 1e-4,
                  diff_method="BD") where {T1<:Real, T2<:Real, T3<:Real}

    p0 = convert(Vector{RheoFloat},p0)
    lo = convert(Vector{RheoFloat},lo)
    hi = convert(Vector{RheoFloat},hi)
    # get modulus function

    if diff_method=="BD"
        deriv = derivBD
    elseif diff_method=="CD"
        deriv = derivCD
    end

    if Int(modloading) == 2
        dcontrolled = deriv(data.σ, data.t)
        measured = data.ϵ
        modtouse = :J
    elseif Int(modloading) == 1
        dcontrolled = deriv(data.ϵ, data.t)
        measured = data.σ
        modtouse = :G
    end
    # get modulus function and derivative
    if (step == nothing)
        check = RheoTimeDataType(data)
        @assert (Int(check) == 3) "Both stress and strain are required"
        if Int(modloading) == 2
            controlled = data.σ[convert(Integer,round(length(data.σ)\2))]
            measured = data.ϵ
            modtouse = :J
        elseif Int(modloading) == 1
            controlled = data.ϵ[convert(Integer,round(length(data.ϵ)/2))]
            measured = data.σ
            modtouse = :G
        end
    elseif (step != nothing)
        check = RheoTimeDataType(data)
        if Int(modloading) == 2
            @assert (Int(check) == 3) || (Int(check) == 1) "Strain required"
            modtouse = :J;
            controlled = convert(RheoFloat,step);
            measured = data.ϵ
        elseif Int(modloading) == 1
            @assert (Int(check) == 3) || (Int(check) == 2) "Stress required"
            measured = data.σ
            modtouse = :G;
            controlled =convert(RheoFloat, step);
        end
    end


    print(controlled)
    modulus = getfield(model, modtouse)

    # use default p0 if not provided
    if quasinull(p0)
        p0 = convert(Array{RheoFloat,1},model.parameters)
    end

    # TEMP - CHECK WITH ALE AND ALEXANDRE BUT IS DEFINITELY NECESSARY
    # time must start at 0 for convolution to work properly!
    t_zeroed = data.t .- minimum(data.t)

    # get singularity presence
    sing = singularitytest(modulus, p0; t1 = t_zeroed[1])

    # start fit
    (minf, minx, ret), timetaken, bytes, gctime, memalloc = @timed leastsquares_stepinit(p0,
                                                                                        lo,
                                                                                        hi,
                                                                                        modulus,
                                                                                        t_zeroed,
                                                                                        controlled,
                                                                                        measured;
                                                                                        insight = verbose,
                                                                                        singularity = sing,
                                                                                        _rel_tol = rel_tol)

    modulusname = string(modulus)


    log = vcat(data.log, "Fitted $modulusname, Step value = $controlled, Time: $timetaken s, Why: $ret, Parameters: $minx, Error: $minf")

    RheologyModel(model.G, model.J, model.Gp, model.Gpp, minx, log)

end

"""
    modelsteppredict(data::RheologyData, model::RheologyModel, modtouse::Symbol; step_on::Real = 0.0)

Same as modelpredict but assumes a step loading with step starting at 'step_on'. Singularities are bypassed
by adding 1 to the index of the singular element.
"""
function modelsteppredict(data, model; modtouse::Symbol=:Nothing, step_on::Real = 0.0, diff_method = "BD")

    step_on = convert(RheoFloat,step_on)

    if diff_method=="BD"
        deriv = derivBD
    elseif diff_method=="CD"
        deriv = derivCD
    end

    check = RheoTimeDataType(data)
    @assert (Int(check) == 1) || (Int(check) == 2) "Both stress and strain are already defined"

    if (modtouse == :Nothing)
        if (Int(check) == 1)
            modtouse = :G;
            controlled = data.ϵ[convert(Integer,round(length(data.ϵ)/2))]
        elseif (Int(check) == 2)
            modtouse = :J;
            controlled = data.σ[convert(Integer,round(length(data.σ)\2))]
        end
    elseif modtouse == :J
        @assert (Int(check) == 2)|| (Int(check) == 3) "Stress required"
        controlled = data.σ[convert(Integer,round(length(data.σ)\2))]
    elseif modtouse == :G
        @assert (Int(check) == 1)|| (Int(check) == 3) "Strain required"
        controlled = data.ϵ[convert(Integer,round(length(data.ϵ)/2))]
    end

    # get modulus
    modulus = getfield(model, modtouse)

    # check singularity presence at time closest to step
    stepon_el = closestindex(data.t, step_on)

    sing = singularitytest(modulus, model.parameters; t1 = convert(RheoFloat,(data.t[stepon_el] - step_on)))

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

    RheoTimeData(convert(Vector{RheoFloat},σ), convert(Vector{RheoFloat},ϵ), data.t, log)

end





function obj_dynamic(params::Vector{RheoFloat},
                     grad::Vector{RheoFloat},
                     ω::Vector{RheoFloat},
                     dataGp::Vector{RheoFloat},
                     dataGpp::Vector{RheoFloat},
                     modelGp::Function,
                     modelGpp::Function;
                     _insight::Bool = false)

    if _insight
        println("Current Parameters: ", params)
    end

    costGp = sum(0.5*(dataGp - modelGp(ω, params)).^2)
    costGpp = sum(0.5*(dataGpp - modelGpp(ω, params)).^2)

    cost = costGp + costGpp

end

function obj_dynamic_linear(params::Vector{RheoFloat},
                            grad::Vector{RheoFloat},
                            ω::Vector{RheoFloat},
                            dataGp::Vector{RheoFloat},
                            dataGpp::Vector{RheoFloat},
                            modelGp::Function,
                            modelGpp::Function,
                            meanGp::RheoFloat,
                            meanGpp::RheoFloat;
                            _insight::Bool = false)

    if _insight
        println("Current Parameters: ", params)
    end

    costGp = sum(0.5*(dataGp/meanGp - modelGp(ω, params)/meanGp).^2)
    costGpp = sum(0.5*(dataGpp/meanGpp - modelGpp(ω, params)/meanGpp).^2)

    cost = costGp + costGpp

end

function obj_dynamic_log(params::Vector{RheoFloat},
                     grad::Vector{RheoFloat},
                     ω::Vector{RheoFloat},
                     dataGp::Vector{RheoFloat},
                     dataGpp::Vector{RheoFloat},
                     modelGp::Function,
                     modelGpp::Function;
                     _insight::Bool = false)

    if _insight
        println("Current Parameters: ", params)
    end

    costGp = sum(0.5*(log.(dataGp) - log.(modelGp(ω, params))).^2)
    costGpp = sum(0.5*(log.(dataGpp) - log.(modelGpp(ω, params))).^2)

    cost = costGp + costGpp

end

function obj_dynamic_global(params::Vector{RheoFloat},
                     grad::Vector{RheoFloat},
                     ω::Vector{RheoFloat},
                     dataGp::Vector{RheoFloat},
                     dataGpp::Vector{RheoFloat},
                     modelGp::Function,
                     modelGpp::Function;
                     _insight::Bool = false)

    if _insight
        println("Current Parameters: ", params)
    end

    costGp = sum(0.5*(((dataGp - modelGp(ω, params))./dataGp).^2))
    costGpp = sum(0.5*(((dataGpp - modelGpp(ω, params))./dataGpp).^2))

    cost = costGp + costGpp

end

function obj_dynamic_manual(params::Vector{RheoFloat},
                            grad::Vector{RheoFloat},
                            ω::Vector{RheoFloat},
                            dataGp::Vector{RheoFloat},
                            dataGpp::Vector{RheoFloat},
                            modelGp::Function,
                            modelGpp::Function,
                            weights::Vector{RheoFloat};
                            _insight::Bool = false)

    if _insight
        println("Current Parameters: ", params)
    end

    costGp = sum(0.5*(dataGp - modelGp(ω, params)).^2)
    costGpp = sum(0.5*(dataGpp - modelGpp(ω, params)).^2)

    cost = weights[1]*costGp + weights[2]*costGpp
end

"""
    dynamicmodelfit(data::RheologyDynamic, model::RheologyModel; p0::Vector{T} = [-1.0], lo::Vector{T} = [-1.0], hi::Vector{T} = [-1.0], verbose::Bool = false, rel_tol = 1e-4) where T<:Real

Fits model to the frequency/loss+storage moduli data.

All arguments are as described below. The 'weights' argument some more information.
As this fitting procedure is fitting two functions simultaneously (the storage
and loss moduli), if left untransformed the fit would tend to favour the
modulus which is larger in magnitude and not fit the other modulus well. To avoid this,
RHEOS offers a number of data transforms which can be used.

# Arguments

- `data`: RheologyDynamic struct containing all data
- `model`: RheologyModel containing moduli and default (initial) parameters
- `p0`: Initial parameters to use in fit (uses 'model' parameters if none given)
- `lo`: Lower bounds for parameters
- `hi`: Upper bounds for parameters
- `verbose`: If true, prints parameters on each optimisation iteration
- `rel_tol`: Relative tolerance of optimization, see NLOpt docs for more details
- `weights`: Weighting mode for storage and loss modulus (see above)
"""
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
        meanGp = sum(data.Gp)/length(data.Gp)
        meanGpp = sum(data.Gp)/length(data.Gp)
        min_objective!(opt, (params, grad) -> obj_dynamic_linear(params, grad, data.ω, data.Gp, data.Gpp, model.Gp, model.Gpp, meanGp, meanGpp; _insight = verbose))

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

"""
    dynamicmodelpredict(data::RheologyDynamic, model::RheologyModel)

Given dynamic rheology data and model, return new dataset based on model parameters.
Returns another RheologyDynamic instance with the predicted Gp and Gpp based on the
frequencies and model parameters.
"""
function dynamicmodelpredict(data::RheologyDynamic, model::RheologyModel)

    # get results
    predGp = model.Gp(data.ω, model.parameters)
    predGpp = model.Gpp(data.ω, model.parameters)

    # store operation
    log = vcat(data.log, "Predicted data from model:", model.log)

    RheologyDynamic(predGp, predGpp, data.ω, log)

end
