#!/usr/bin/env julia

#############################
#~ Preprocessing Functions ~#
#############################

export cutting, extract

"""
    fixedresample(self::RheoTimeData, elperiods::Union{Vector{K},K}; time_boundaries::Vector{T}= [-1])

Resample data with new sample rate(s).

Fixedresample can downsample or upsample data. If the number of elperiods is negative it is going to reduce the number of samples,
viceversa if it is positive. If time boundaries are not specified, resampling is applied to the whole set of data.
If number of elements per period (elperiods) is 1 or -1 it returns the original RheoTimeData, whilst 0 is not accepted as valid
number for elperiods.
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
    if (check == strain_only)
        (time, epsilon) = fixed_resample(self.t, self.ϵ, boundaries, elperiods)
        sigma = [];
    elseif (check == stress_only)
        (time, sigma) = fixed_resample(self.t, self.σ, boundaries, elperiods)
        epsilon = [];
    elseif (check == strain_and_stress)
        (time, sigma) = fixed_resample(self.t, self.σ, boundaries, elperiods)
        (time, epsilon) = fixed_resample(self.t, self.ϵ, boundaries, elperiods)
    end

    # add record of operation applied
    log = vcat(self.log, "fixed_resample - boundaries: $boundaries, elperiods: $elperiods")

    self_new = RheoTimeData(sigma, epsilon, time, log)

end

"""
    function cutting(self::RheoTimeData, time_on::T1,time_off::T2) where {T1<:Number, T2<:Number}

Remove the data outside a specified time interval.

By specifing a time interval (time_on, time_off), a new RheoTimeData is returned without the data
laying outside time interval.
"""
function cutting(self::RheoTimeData, time_on::T1,time_off::T2) where {T1<:Number, T2<:Number}

    @assert isreal(time_on) && isreal(time_off) "Boundaries cannot be complex numbers"

    boundary_on = closestindex(self.t, time_on);
    boundary_off = closestindex(self.t, time_off);
    time = self.t[boundary_on:boundary_off]

    check = RheoTimeDataType(self)
    if (check == strain_only)
        epsilon = self.ϵ[boundary_on:boundary_off]
        sigma = [];
    elseif (check == stress_only)
        sigma = self.σ[boundary_on:boundary_off]
        epsilon = [];
    elseif (check == strain_and_stress)
        epsilon = self.ϵ[boundary_on:boundary_off]
        sigma = self.σ[boundary_on:boundary_off]
    end
    log = vcat(self.log, "Data from $time_on to $time_off extracted.")

    return RheoTimeData(sigma,epsilon,time,log)

end

"""
function smooth(self::RheoTimeData, τ::Real; pad::String="reflect")

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
    if (check == strain_only)
        epsilon = Base.invokelatest(imfilter, self.ϵ, Base.invokelatest(Kernel.reflect, Base.invokelatest(Kernel.gaussian, (Σ,))), pad)
        sigma = [];
    elseif (check == stress_only)
        sigma = Base.invokelatest(imfilter, self.σ, Base.invokelatest(Kernel.reflect, Base.invokelatest(Kernel.gaussian, (Σ,))), pad)
        epsilon = [];
    elseif (check == strain_and_stress)
        sigma = Base.invokelatest(imfilter, self.σ, Base.invokelatest(Kernel.reflect, Base.invokelatest(Kernel.gaussian, (Σ,))), pad)
        epsilon = Base.invokelatest(imfilter, self.ϵ, Base.invokelatest(Kernel.reflect, Base.invokelatest(Kernel.gaussian, (Σ,))), pad)
    end

    # add record of operation applied
    log = vcat(self.log, "smooth - τ: $τ")

    self_new = RheoTimeData(sigma, epsilon, self.t, log)

end

"""
    function extract(self::Union{RheoTimeData,RheoFreqData}, type::Union{TimeDataType,FreqDataType,Integer})

Extract specific fields form RheoTimeData or RheoFreqData.

Extract can copy one or more fields from a given RheoXData variable into a new RheoXData one. The fields that
are copied are identified by the specified type of data.
If self is a RheoTimeData, the type that can be extracted is time_only (or 0), stress_only (or 1), strain_only (or 2).
Note that strain_and_stress (or 3) is not allowed.
If self is a RheoFreqData, the type that can be extracted is frec_only (or 0).
"""
function extract(self::Union{RheoTimeData,RheoFreqData}, type::Union{TimeDataType,FreqDataType,Integer})

    if (typeof(self)==RheoTimeData)

        type = typeof(type)==Int ? TimeDataType(type) : type
        @assert (typeof(type)==TimeDataType) || (typeof(type)==Int) "Cannot extract frequency data from RheoTimeData"
        @assert (type!= strain_and_stress) && (type!= invalid_time_data) "Cannot extract both stress and strain"
        @assert (type!= invalid_time_data) "Cannot extract information from invalid time data"
        check = RheoTimeDataType(self)

        if type == time_only
            @assert check!= invalid_time_data "Time not available"
            log = vcat(self.log,"Time extracted")
            return RheoTimeData([], [],self.t,log)
        elseif type == strain_only
            @assert (check == strain_and_stress) || (check == strain_only) "Strain not available"
            log = vcat(self.log,"Time and strain extracted")
            return RheoTimeData([], self.ϵ,self.t,log)
        elseif type == stress_only
            @assert (check == strain_and_stress) || (check == stress_only) "Stress not available"
            log = vcat(self.log,"Time and stress extracted")
            return RheoTimeData(self.σ, [], self.t,log)
        end

    elseif (typeof(self)==RheoFreqData)

        type = typeof(type)==Int ? FreqDataType(type) : type
        @assert (typeof(type)==FreqDataType) || (typeof(type)==Int) "Cannot extract time data from RheoFreqData"
        @assert (type!= with_modulus) "Cannot extract frequency with moduli"
        @assert (type!= invalid_freq_data) "Cannot extract information from invalid frequency data"
        check = RheoFreqDataType(self)
        @assert (check == with_modulus) "Frequency and modulii required"
        log = vcat(self.log,"Frequency extracted")
        return RheoFreqData([], [],self.ω,log)
    end

end




##########################
#~ Processing Functions ~#
##########################

"""
    modelfit(data::RheoTimeData, model::RheoModelClass, modloading::Symbol; p0::Union{NamedTuple,Tuple} = (), lo::Union{NamedTuple,Tuple} = (), hi::Union{NamedTuple,Tuple} = (), verbose::Bool = false, rel_tol = 1e-4, diff_method="BD")

Fit RheologyData struct to model and return a fitted model as a RheologyModel object.

# Arguments

- `data`: RheoTimeData struct containing all data
- `model`: RheoModelClass containing moduli functions and named tuple parameters
- `modloading`: strain_imposed or 1, stress_imposed or 2
- `p0`: Initial parameters to use in fit (uses 0.5 for all parameters if not defined)
- `lo`: Lower bounds for parameters
- `hi`: Upper bounds for parameters
- `verbose`: If true, prints parameters on each optimisation iteration
- `rel_tol`: Relative tolerance of optimization, see NLOpt docs for more details
- `diff_method`: Set finite difference formula to use for derivative, currently "BD" or "CD"
"""
function modelfit(data::RheoTimeData,
                  model::RheoModelClass,
                  modloading::Union{LoadingType,Integer};
                  p0::Union{NamedTuple,Tuple} = (),
                  lo::Union{NamedTuple,Tuple} = (),
                  hi::Union{NamedTuple,Tuple} = (),
                  verbose::Bool = false,
                  rel_tol = 1e-4,
                  diff_method="BD")

    if isempty(p0)
        p0a = convert(Array{RheoFloat,1}, fill(0.5,length(model.params)))
        if length(findall(x->(x==:β)||(x==:a),model.params))==2
            index = findall(x->x==:a,model.params)
            p0a[index[1]] = 0.8;
        end
        @warn "Initial values for model parameters is set to $p0a by default"
    else
        p0a = model_parameters(p0,model.params,"initial guess")
    end

    #check_ineq = model.ineq(p0)
    #@assert (all(check_ineq)==true)  "Initial guess not feasible"
    @assert model.constraint(p0a)  "Initial guess not feasible"

    if isempty(lo)
        loa = convert(Array{RheoFloat,1}, [-1])
    else
        loa = model_parameters(lo,model.params,"low bounds")
    end

    if isempty(hi)
        hia = convert(Array{RheoFloat,1}, [-1])
    else
        hia = model_parameters(hi,model.params,"high bounds")
    end

    rel_tol = convert(RheoFloat,rel_tol)

    check = RheoTimeDataType(data)
    @assert (check == strain_and_stress) "Both stress and strain are required"


    # use correct method for derivative
    if diff_method=="BD"
        deriv = derivBD
    elseif diff_method=="CD"
        deriv = derivCD
    end

    #modloading = typeof(modloading)==Int ? LoadingType(modloading) : modloading
    # get modulus function and derivative
    if modloading == stress_imposed
        dcontrolled = deriv(data.σ, data.t)
        measured = data.ϵ
        modulus = model.Ja
        modsing = (t->model.J(t,p0a))
    elseif modloading == strain_imposed
        dcontrolled = deriv(data.ϵ, data.t)
        measured = data.σ
        modulus = model.Ga
        modsing = (t->model.G(t,p0a))
    end

    @assert ~isnan(modsing(5.0)) "Modulus to use not defined"
    # get singularity presence
    sing = singularitytest(modsing)
    # get time step (only needed for convolution, which requires constant dt so t[2]-t[1] is sufficient)
    dt = data.t[2] - data.t[1]

    # TEMP - CHECK WITH ALE AND ALEXANDRE BUT IS DEFINITELY NECESSARY
    # time must start at 0 for convolution to work properly!
    t_zeroed = data.t .- minimum(data.t)

    # fit
    sampling_check = constantcheck(data.t)

    (minf, minx, ret), timetaken, bytes, gctime, memalloc = @timed leastsquares_init(p0a,
                                                                                    loa,
                                                                                    hia,
                                                                                    modulus,
                                                                                    t_zeroed,
                                                                                    dt,
                                                                                    dcontrolled,
                                                                                    measured;
                                                                                    insight = verbose,
                                                                                    sampling = sampling_check,
                                                                                    singularity = sing,
                                                                                    _rel_tol = rel_tol)

    #modulusname = string(modulus)
    #log = vcat(data.log, "Fitted $modulusname, Modulus used: $modtouse, Time: $timetaken s, Why: $ret, Parameters: $minx, Error: $minf")
    log = vcat(data.log, "Time: $timetaken s, Why: $ret, Parameters: $minx, Error: $minf")
    print("Time: $timetaken s, Why: $ret, Parameters: $minx, Error: $minf")
    nt = NamedTuple{Tuple(model.params)}(minx)

    return RheoModel(model,nt; log_add = log);

end

"""
    modelpredict(data::RheoTimeData,model::RheoModel; diff_method="BD")

Given an incomplete data set (only either stress or strain missing) and model with values substituted into
parameters (RheoModel),return a new dataset based on the model.
If data is type of stress_only, then creep modulus (:J) is used; if data type is strain_only relaxation modulus (:G).
A complete RheoTimeDatadata of type "strain_and_stress" is returned.
'diff_method' sets finite difference for calculating the derivative used in the hereditary integral and
can be either backwards difference ("BD") or central difference ("CD").
"""
function modelpredict(data::RheoTimeData,model::RheoModel; diff_method="BD")

    # use correct method for derivative
    if diff_method=="BD"
        deriv = derivBD
    elseif diff_method=="CD"
        deriv = derivCD
    end

    check = RheoTimeDataType(data)
    @assert (check == strain_only)||(check == stress_only) "Need either strain only or stress only data. Data provide: " * string(check)

    if (check == strain_only)
        modulus = model.Ga
        modsing = model.G
        dcontrolled = deriv(data.ϵ, data.t)
    elseif (check == stress_only)
        modulus = model.Ja
        modsing = model.J
        dcontrolled = deriv(data.σ, data.t)
    end

    # get singularity presence
    sing = singularitytest(modsing)

    # get time step (only needed for convolution, which requires constant so t[2]-t[1] is sufficient)
    dt = data.t[2] - data.t[1]

    # TEMP - CHECK WITH ALE AND ALEXANDRE BUT IS DEFINITELY NECESSARY
    # time must start at 0 for convolution to work properly!
    t_zeroed = data.t .- minimum(data.t)

    # get convolution
    if !sing && constantcheck(data.t)
        convolved = boltzconvolve_nonsing(modulus, t_zeroed, dt, dcontrolled)

    elseif sing && constantcheck(data.t)
        # convolved = boltzconvolve_sing(modulus, t_zeroed, dt, dcontrolled)
        t_zeroed[1] = 0.0 + (t_zeroed[2] - t_zeroed[1])/10.0
        # convolved = boltzconvolve_sing(modulus, t_zeroed, dt, dcontrolled)
        convolved = boltzconvolve(modulus, t_zeroed, dt, dcontrolled)

    elseif !sing && !constantcheck(data.t)
        convolved = boltzintegral_nonsing(modulus, t_zeroed, dcontrolled)

    elseif sing && !constantcheck(data.t)
        convolved = boltzintegral_sing(modulus, t_zeroed, dcontrolled)

    end

    if (check == stress_only)
        sigma = data.σ
        epsilon = convolved
        pred_mod = "creep function J"
    elseif (check == strain_only)
        sigma = convolved
        epsilon = data.ϵ
        pred_mod = "relaxation function G"
    end
    time = data.t

    modparam = model.params
    # store operation
    log = vcat( data.log, "Predicted data using: $pred_mod, Parameters: $modparam")


    return RheoTimeData(sigma,epsilon,time, log)

end

"""
    modelstepfit(data::RheoTimeData, model::RheoModelClass, modloading::Union{LoadingType,Integer}; step=nothing, p0::Union{NamedTuple,Tuple} = (), lo::Union{NamedTuple,Tuple} = (), hi::Union{NamedTuple,Tuple} = (), verbose::Bool = false, rel_tol = 1e-4) where T<:Real

Same as 'modelfit' except assumes a step loading. If this assumption is appropriate for the data
then fitting can be sped up greatly by use of this function. If modloading is strain_imposed, relaxation modulus is used,
then the element in the middle of the strain is assumed to be the amplitude of the step. If modloading is stress_imposed,
the creep modulus is used, then the middle element of the stress is assumed to be the amplitude of the step.
Alternatively, it is possible to define the value of the step by defining the optional "step" parameter.

# Arguments

- `data`: RheoTimeData struct containing all data
- `model`: RheoModelClass containing moduli and parameters tuples
- `modloading`: strain_imposed for relaxation modulus, stress_imposed for creep modulus
- `p0`: Named tuple of initial parameters to use in fit (uses 0.5 for all parameters if none given)
- `lo`: Named tuple of lower bounds for parameters
- `hi`: Named tuple of upper bounds for parameters
- `verbose`: If true, prints parameters on each optimisation iteration
- `rel_tol`: Relative tolerance of optimization, see NLOpt docs for more details
- `diff_method`: Set finite difference formula to use for derivative, currently "BD" or "CD"
"""
function modelstepfit(data::RheoTimeData,
                  model::RheoModelClass,
                  modloading::Union{LoadingType,Integer};
                  step = nothing,
                  p0::Union{NamedTuple,Tuple} = (),
                  lo::Union{NamedTuple,Tuple} = (),
                  hi::Union{NamedTuple,Tuple} = (),
                  verbose::Bool = false,
                  rel_tol = 1e-4,
                  diff_method="BD") where {T1<:Real, T2<:Real, T3<:Real}

    if isempty(p0)
     p0 = convert(Array{RheoFloat,1}, fill(0.5,length(model.params)))
          @warn "Initial values for model parameters is set to 0.5 by default"
    else
     p0 = model_parameters(p0,model.params,"initial guess")
    end

    if isempty(lo)
     lo = convert(Array{RheoFloat,1}, [-1])
    else
     lo = model_parameters(lo,model.params,"low bounds")
    end

    if isempty(hi)
     hi = convert(Array{RheoFloat,1}, [-1])
    else
     hi = model_parameters(hi,model.params,"high bounds")
    end

    rel_tol = convert(RheoFloat,rel_tol)

    if diff_method=="BD"
        deriv = derivBD
    elseif diff_method=="CD"
        deriv = derivCD
    end

    if  modloading == stress_imposed
        dcontrolled = deriv(data.σ, data.t)
        measured = data.ϵ
        modulus = model.Ja
    elseif modloading == strain_imposed
        dcontrolled = deriv(data.ϵ, data.t)
        measured = data.σ
        modulus = model.Ga
    end
    # get modulus function and derivative
    if (step == nothing)
        check = RheoTimeDataType(data)
        @assert (check == strain_and_stress) "Both stress and strain are required"
        if (modloading == stress_imposed)
            controlled = data.σ[convert(Integer,round(length(data.σ)/2))]
            measured = data.ϵ
            modulus = model.Ja
        elseif (modloading == strain_imposed)
            controlled = data.ϵ[convert(Integer,round(length(data.ϵ)/2))]
            measured = data.σ
            modulus = model.Ga
        end
    elseif (step != nothing)
        check = RheoTimeDataType(data)
        if (modloading == stress_imposed)
            @assert (check == strain_only)||(check == strain_and_stress) "Strain required"
            modulus = model.Ja
            controlled = convert(RheoFloat,step);
            measured = data.ϵ
        elseif (modloading == strain_imposed)
            @assert (check == stress_only)||(check == strain_and_stress) "Stress required"
            measured = data.σ
            modulus = model.Ga
            controlled =convert(RheoFloat, step);
        end
    end

    mod(t) = modulus(t,p0)

    @assert ~isnan(mod(5.0)) "Modulus to use not defined"
    # get singularity presence
    sing = singularitytest(mod)

    # TEMP - CHECK WITH ALE AND ALEXANDRE BUT IS DEFINITELY NECESSARY
    # time must start at 0 for convolution to work properly!
    t_zeroed = data.t .- minimum(data.t)


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
    log = vcat(data.log, "Fitted $modulusname, Modulus used: $modtouse, Time: $timetaken s, Why: $ret, Parameters: $minx, Error: $minf")
    nt = NamedTuple{Tuple(model.params)}(minx)

    return RheoModel(model,nt; log_add = log);

end

"""
    modelsteppredict(data::RheoTimeData, model::RheoModel; step_on::Real = 0.0, diff_method = "BD")

Same as modelpredict but assumes a step loading with step starting at 'step_on'. Singularities are bypassed
by adding 1 to the index of the singular element.
"""
function modelsteppredict(data::RheoTimeData, model::RheoModel; step_on::Real = 0.0, diff_method = "BD")

    step_on = convert(RheoFloat,step_on)

    if diff_method=="BD"
        deriv = derivBD
    elseif diff_method=="CD"
        deriv = derivCD
    end

    check = RheoTimeDataType(data)
    @assert (check == strain_only)||(check == stress_only) "Need either strain only or stress only data. Data provide: " * string(check)

    if (check == strain_only)
        modulus = model.Ga
        controlled = data.ϵ[convert(Integer,round(length(data.ϵ)/2))]
    elseif (check == stress_only)
        modulus = model.Ja
        print(round(length(data.σ)))
        controlled = data.σ[convert(Integer,round(length(data.σ)/2))]
    end

    # check singularity presence at time closest to step
    stepon_el = closestindex(data.t, step_on)

    sing = singularitytest(modulus; t1 = convert(RheoFloat,(data.t[stepon_el] - step_on)))

    # get predicted
    if !sing
        predicted = zeros(length(data.t))
        predicted[stepon_el:end] = controlled*modulus(data.t[stepon_el:end] .- step_on)

    elseif sing
        predicted = zeros(length(data.t))
        predicted[(stepon_el + 1):end] = controlled*modulus(data.t[(stepon_el + 1):end] .- step_on)

    end

    if check == stress_only
        σ = data.σ
        ϵ = predicted
    elseif check == strain_only
        σ = predicted
        ϵ = data.ϵ
    end

    # store operation
    log = vcat(data.log, "Predicted data from model:", model.log)

    RheoTimeData(convert(Vector{RheoFloat},σ), convert(Vector{RheoFloat},ϵ), data.t, log)

end





function obj_dynamic(params, grad, ω, dataGp, dataGpp, modelGp, modelGpp; _insight::Bool = false)

    if _insight
        println("Current Parameters: ", params)
    end

    modGp(ω) = modelGp(ω,params)
    modGpp(ω) = modelGpp(ω,params)
    costGp = sum(0.5*(dataGp - modGp(ω)).^2)
    costGpp = sum(0.5*(dataGpp - modGpp(ω)).^2)

    cost = costGp + costGpp

end

function obj_dynamic_linear(params, grad, ω, dataGp, dataGpp, modelGp, modelGpp, meanGp, meanGpp; _insight::Bool = false)

    if _insight
        println("Current Parameters: ", params)
    end

    modGp(ω) = modelGp(ω,params)
    modGpp(ω) = modelGpp(ω,params)
    costGp = sum(0.5*(dataGp/meanGp - modGp(ω)/meanGp).^2)
    costGpp = sum(0.5*(dataGpp/meanGpp - modGpp(ω)/meanGpp).^2)

    cost = costGp + costGpp

end

function obj_dynamic_log(params, grad, ω, dataGp, dataGpp, modelGp, modelGpp; _insight::Bool = false)

    if _insight
        println("Current Parameters: ", params)
    end
    modGp(ωa) = modelGp(ωa,params)
    modGpp(ωa) = modelGpp(ωa,params)
    costGp = sum(0.5*(log.(dataGp) - log.(modGp(ω))).^2)
    costGpp = sum(0.5*(log.(dataGpp) - log.(modGpp(ω))).^2)

    cost = costGp + costGpp

end

function obj_dynamic_global(params, grad, ω, dataGp, dataGpp, modelGp, modelGpp; _insight::Bool = false)

    if _insight
        println("Current Parameters: ", params)
    end
    modGp(ωa) = modelGp(ωa,params)
    modGpp(ωa) = modelGpp(ωa,params)
    costGp = sum(0.5*(((dataGp - modGp(ω))./dataGp).^2))
    costGpp = sum(0.5*(((dataGpp - modGpp(ω))./dataGpp).^2))

    cost = costGp + costGpp

end

function obj_dynamic_manual(params, grad, ω, dataGp, dataGpp, modelGp, modelGpp, weights; _insight::Bool = false)

    if _insight
        println("Current Parameters: ", params)
    end
    modGp(ωa) = modelGp(ωa,params)
    modGpp(ωa) = modelGpp(ωa,params)
    costGp = sum(0.5*(dataGp - modGp(ω)).^2)
    costGpp = sum(0.5*(dataGpp - modGpp(ω)).^2)

    cost = weights[1]*costGp + weights[2]*costGpp
end

"""
    dynamicmodelfit(data::RheoFreqData, model::RheoModelClass; p0::Union{NamedTuple,Tuple} = (), lo::Union{NamedTuple,Tuple} = (), hi::Union{NamedTuple,Tuple} = (), verbose::Bool = false, rel_tol = 1e-4) where T<:Real

Fits model to the frequency/loss+storage moduli data.

All arguments are as described below.
As this fitting procedure is fitting two functions simultaneously (the storage
and loss moduli), if left untransformed the fit would tend to favour the
modulus which is larger in magnitude and not fit the other modulus well. To avoid this,
RHEOS offers a number of data transforms which can be used by changing "weights" argument.

# Arguments

- `data`: RheoFreqData struct containing all data
- `model`: RheoModelClass containing moduli and symbols of parameters
- `p0`: Initial parameters to use in fit (uses 0.5 for all parameters if none given)
- `lo`: Lower bounds for parameters
- `hi`: Upper bounds for parameters
- `verbose`: If true, prints parameters on each optimisation iteration
- `rel_tol`: Relative tolerance of optimization, see NLOpt docs for more details
- `weights`: Weighting mode for storage and loss modulus (linear, log, global)
"""
function dynamicmodelfit(data::RheoFreqData,
                model::RheoModelClass;
                p0::Union{NamedTuple,Tuple} = (),
                lo::Union{NamedTuple,Tuple} = (),
                hi::Union{NamedTuple,Tuple} = (),
                verbose::Bool = false,
                rel_tol::T = 1e-4,
                weights::Union{String, Vector{T}}="log") where T<:Real


    if isempty(p0)
       p0 = convert(Array{RheoFloat,1}, fill(0.5,length(model.params)))
       @warn "Initial values for model parameters is set to 0.5 by default"
    else
       p0 = model_parameters(p0,model.params,"initial guess")
    end

    if isempty(lo)
       lo = convert(Array{RheoFloat,1}, [-1])
    else
       lo = model_parameters(lo,model.params,"low bounds")
    end

    if isempty(hi)
       hi = convert(Array{RheoFloat,1}, [-1])
    else
       hi = model_parameters(hi,model.params,"high bounds")
    end

    rel_tol = convert(RheoFloat,rel_tol)

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


    modelname = string(model)
    log = vcat(data.log, "Fitted Gp, Gpp of $modelname, Time: $timetaken s, Why: $ret, Parameters: $minx, Error: $minf")
    nt = NamedTuple{Tuple(model.params)}(minx)

    return RheoModel(model,nt; log_add = log);

end

"""
    dynamicmodelpredict(data::RheoFreqData, model::RheoModel)

Given dynamic rheology data with only frequency and model where parameters have been substituted.
Returns another RheoFreqData instance with the predicted Gp and Gpp based on the
frequencies and model given as arguments.
"""
function dynamicmodelpredict(data::RheoFreqData, model::RheoModel)

    # get results
    predGp = model.Gp(data.ω)
    predGpp = model.Gpp(data.ω)

    # store operation
    log = vcat(data.log, "Predicted data from model:", model.log)
    RheoFreqData(predGp, predGpp,data.ω, log)

end
