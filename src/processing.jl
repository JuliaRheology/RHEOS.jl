#!/usr/bin/env julia

#############################
#~ Preprocessing Functions ~#
#############################

"""
    variableresample(self::RheologyData, refvar::Symbol, pcntdownsample::Real; mapback::Bool = false)

Convert a fixed sample rate array to a variable sample rate, with sampling points
added according to a relative change in chosen variable `refvar`, 1st derivative
of `refvar` and 2nd derivative of `refvar` (WRT time). Usually chosen as the
measured variable, so :σ for a stress relaxation test and :ϵ for a creep test.

Currently only variable downsampling supported. pcntdown sample is approximate,
works well in some cases and very poorly in others. If required, compare resampled
length vs original length after processing has finished. If data is noisy, may
benefit from sending smoothed signal to this algorithm and either using mapback
function or interpolating onto unsmoothed data.

Algorithm works as follows. 25 initial samples are generated evenly spread. After this,
the array is repeatedly sweeped, anywhere Δy, Δdy/dx, Δd2y/dx2 is greater than a threshold, α,
a new sample is created at the midpoint of the two tested points. This is allowed to
happen a maximum of 400 times, after which α is decreased and the process starts again.
This macro process continues until the desired pcntdownsample ratio has been reached.

# Arguments

- `self`: RheologyData instance
- `refvar`: The data whose derivatives will determine sample densities
- `pcntdownsample`: Approximate ratio of new samples to old samples
- `_mapback = false`: (Optional) Determines whether resampled points should 'snap' to closest original points
"""
function variableresample(self::RheologyData, refvar::Symbol, pcntdownsample::Real; _mapback::Bool = true)

    @eval import Interpolations: interpolate, Gridded, Linear

    # enforce minimum period of original period/10
    _minperiod = (self.t[2] - self.t[1])/10.0

    # get time resampled with respect to refvar
    (tᵦ, dummy)  = var_resample(self.t, getfield(self, refvar), pcntdownsample, _minperiod)

    if _mapback
        # get mapped indices wrt original
        mapback_indices = mapback(tᵦ, self.t)

        σ = self.σ[mapback_indices]
        ϵ = self.ϵ[mapback_indices]
        t = self.t[mapback_indices]

    elseif !_mapback
        # interpolate with respect to t
        σ_interp = Base.invokelatest(interpolate, (self.t,), self.σ, Base.invokelatest(Gridded, Base.invokelatest(Linear)))
        ϵ_interp = Base.invokelatest(interpolate, (self.t,), self.σ, Base.invokelatest(Gridded, Base.invokelatest(Linear)))

        # resample all using new timepoints tᵦ
        σ = σ_interp[tᵦ]
        ϵ = ϵ_interp[tᵦ]
        t = tᵦ
    end

    # change sampling type to variable
    sampling = "variable"

    # add record of operation applied
    log = vcat(self.log, "var_resample - refvar: $refvar, pcntdownsample: $pcntdownsample, mapback: $_mapback")

    self_new = RheoTimeData(σ=σ, ϵ=ϵ, t=t, log)

end

"""
    fixedresample(self::RheoTimeData, elperiods::Union{Vector{K},K}; time_boundaries::Vector{T}= [-1])

Resample data with new sample rate(s).

Fixedresample can downsample or upsample data. If the number of elperiods is negative it is going to reduce the number of samples,
viceversa if it is positive. If time boundaries are not specified, resampling is applied to the whole set of data.
"""
function fixedresample(self::RheoTimeData, elperiods::Union{Vector{K},K}; time_boundaries::Vector{T}= [-1]) where {K<:Integer,T<:Real}

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

"""
    smooth(self::RheologyData, τ::Real; pad::String="replicate")

Smooth data using a Gaussian Kernel to time scale τ (approximately half power).

Smooths both σ and ϵ. Essentially a low pass filter with frequencies of 1/τ being cut to approximately
half power. For other pad types available see ImageFiltering documentation.
"""
function smooth(self::RheologyData, τ::Real; pad::String="reflect")

    @eval import ImageFiltering: imfilter, Kernel

    # get sample-rate and Gaussian kernel (std. dev)
    samplerate = 1.0/getsampleperiod(self.t)
    Σ = getsigma(τ, samplerate)

    # smooth signal and return
    σ = Base.invokelatest(imfilter, self.σ, Base.invokelatest(Kernel.reflect, Base.invokelatest(Kernel.gaussian, (Σ,))), pad)
    ϵ = Base.invokelatest(imfilter, self.ϵ, Base.invokelatest(Kernel.reflect, Base.invokelatest(Kernel.gaussian, (Σ,))), pad)

    # add record of operation applied
    log = vcat(self.log, "smooth - τ: $τ")

    self_new = RheoTimeData(σ=σ, ϵ=ϵ, t=self.t, log)

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

    self_new = RheoTimeData(σ=σ, ϵ=ϵ,t = t, log)

end

"""
    function zerotime(self::RheologyData)

Convenience function to normalize time such that the starting time is 0.0
"""
function zerotime(self::RheologyData)

    return RheologyData(self.σ, self.ϵ, self.t .- minimum(self.t), self.sampling, vcat(self.log, ["Normalized time to start at 0.0"]))

end


# """
#     downsample(self::RheologyData, time_boundaries::Vector{T} where T<:Real, elperiods::Vector{S} where S<:Integer)
#
# Boundaries are floating point times which are then converted to the closest elements. Works by just
# reducing on indices. For example, time_boundaries could be [0.0, 10.0, 100.0] and elperiods could be
# [2, 4]. So new data would take every 2nd element from 0.0 seconds to 10.0 seconds, then every 4th element
# from 10.0 seconds to 100.0 seconds.
# """
# function downsample(self::RheologyData, time_boundaries::Vector{T} where T<:Real, elperiods::Vector{S} where S<:Integer)
#
#     # convert boundaries from times to element indicies
#     boundaries = closestindices(self.t, time_boundaries)
#
#     # get downsampled indices
#     indices = downsample(boundaries, elperiods)
#
#     # downsample data
#     σ = self.σ[indices]
#     ϵ = self.ϵ[indices]
#     t = self.t[indices]
#
#     # change to variable sampling rate if more than one section, if not then as original
#     local sampling::String
#     if length(elperiods) > 1
#         sampling = "variable"
#     else
#         sampling = self.sampling
#     end
#
#     # add record of operation applied
#     log = vcat(self.log, "downsample - boundaries: $boundaries, elperiods: $elperiods")
#
#     self_new = RheoTimeData(σ=σ, ϵ=ϵ, t=t, log)
#
# end

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
                  model::RheologyModel;
                  p0::Array{T1,1} = [-1.0],
                  lo::Array{T2,1} = [-1.0],
                  hi::Array{T3,1} = [-1.0],
                  modtouse::Symbol = :Nothing,
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
    if modtouse == :Nothing
        dσ = deriv(data.σ, data.t)
        dϵ = deriv(data.ϵ, data.t)
        max_dσ = maximum(dσ[2:end])
        max_dϵ = maximum(dϵ[2:end])
        if (max_dσ<=max_dϵ)
            modtouse = :J;
            dcontrolled = dσ;
            measured = data.ϵ
        elseif (max_dσ>max_dϵ)
            modtouse = :G;
            dcontrolled = dϵ;
            measured = data.σ
        end
        print(modtouse)
    elseif modtouse == :J
        dcontrolled = deriv(data.σ, data.t)
        measured = data.ϵ
    elseif modtouse == :G
        dcontrolled = deriv(data.ϵ, data.t)
        measured = data.σ
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
function modelpredict(data::RheoTimeData,model::RheologyModel; modtouse::Symbol=:Nothing, diff_method="BD")

    # use correct method for derivative
    if diff_method=="BD"
        deriv = derivBD
    elseif diff_method=="CD"
        deriv = derivCD
    end

    if modtouse == :Nothing
        check = RheoTimeDataType(data)
        if (Int(check) == 1)
            modtouse = :G;
            dcontrolled = deriv(data.ϵ, data.t)
        elseif (Int(check) == 2)
            modtouse = :J;
            dcontrolled = deriv(data.σ, data.t)
        elseif (Int(check) == 3)
            dσ = deriv(data.σ, data.t)
            dϵ = deriv(data.ϵ, data.t)
            max_dσ = maximum(dσ[2:end])
            max_dϵ = maximum(dϵ[2:end])
            if (max_dσ<=max_dϵ)
                modtouse = :J;
                dcontrolled = dσ;
            elseif (max_dσ>max_dϵ)
                modtouse = :G;
                dcontrolled = dϵ;
            end
        end
    elseif modtouse == :J
        dcontrolled = deriv(data.σ, data.t)
    elseif modtouse == :G
        dcontrolled = deriv(data.ϵ, data.t)
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

    # ASK LOUIS!
    # if sing
    #     if modtouse == :J
    #         sigma = data.σ
    #         epsilon = convolved
    #     elseif modtouse == :G
    #         sigma = convolved
    #         epsilon = data.ϵ
    #     end
    #     time = data.t
    #
    # elseif !sing
    #     if modtouse == :J
    #         sigma = data.σ
    #         epsilon = convolved
    #     elseif modtouse == :G
    #         sigma = convolved
    #         epsilon = data.ϵ
    #     end
    #     time = data.t
    #
    # end

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
function modelstepfit(data::RheologyData,
                  model::RheologyModel,
                  modtouse::Symbol;
                  p0::Array{T1,1} = [-1.0],
                  lo::Array{T2,1} = [-1.0],
                  hi::Array{T3,1} = [-1.0],
                  verbose::Bool = false,
                  rel_tol = 1e-4)::RheologyModel  where {T1<:Real, T2<:Real, T3<:Real}

    p0 = convert(Vector{RheoFloat},p0)
    lo = convert(Vector{RheoFloat},lo)
    hi = convert(Vector{RheoFloat},hi)
    rel_tol = convert(Vector{RheoFloat},rel_tol)
    # get modulus function
    modulus = getfield(model, modtouse)

    # use default p0 if not provided
    if quasinull(p0)
        p0 = convert(Array{RheoFloat,1},model.parameters)
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

"""
    modelsteppredict(data::RheologyData, model::RheologyModel, modtouse::Symbol; step_on::Real = 0.0)

Same as modelpredict but assumes a step loading with step starting at 'step_on'. Singularities are bypassed
by adding 1 to the index of the singular element.
"""
function modelsteppredict(data::RheologyData, model::RheologyModel, modtouse::Symbol; step_on::Real = 0.0)

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

    RheoTimeData(σ=σ, ϵ=ϵ, t=t, log)

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
