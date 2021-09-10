#!/usr/bin/env julia

#=
-----------------------
Preprocessing functions
-----------------------
=#

"""
    resample(d::RheoTimeData [, t, dt, scale )

Resample the data using 1D spline extrapolation (using the Dierckx.jl package).

# Arguments

- `d`: data as a `RheoTimeData` struct containing time with stress and/or strain. Without additional parameters, the data is resampled to provide a uniform sampling at constant number of timepoints.
- `t`: an array or range that determines the timepoints where the data will be provided.
- `dt`: if the array `t` is not provided, the parameter `dt` will set the timestep for a uniform resampling of the data.
- `scale`: instead of specifying particular time points or timestep, an overall multiplicator on the sampling rate can be provided. This could down-sample (`scale`<1) or upsample (`scale`>1). If timesteps are non uniform, it would interpolate values accordingly.

# Examples

Assuming `d` is a `RheoTimeData` data set:
- `resample(d)` keeps the number of sampling points the same but interpolates to set a uniform time step.
- `resample(d, t=-1:0.1:10)` resamples by interpolation to generate a new dataset with time points given by the range `t`.
- `resample(d, dt=0.1)` resamples by interpolation to generate a new dataset with uniform time step `dt`.
- `resample(d, scale=2)` resamples by multiplying the sampling rate by 2.
"""
function resample(d::RheoTimeData; t::Union{Vector{T},R}=RheoFloat[], scale::T1=0, dt::T2=0) where {T<:Real, R <: AbstractRange, T1<:Real, T2<:Real}

    @assert hastime(d) "Data without time information cannot be resampled."
    if length(t)==0
        if scale > 0
            δt = (typeof(scale)<:Union{Rational,Int}) ? 1//scale : 1. /scale
            idxt=Array(1:length(d.t))
            newidxt=Array(1:δt:length(d.t))
            t=Spline1D(idxt,d.t)(newidxt)
            keywords=(scale = scale,)
        else
            mn, mx = extrema(d.t)
            dt==0 ? δt=(mx-mn)/(length(d.t)-1) : δt=dt
            t=Array{RheoFloat}(mn:δt:mx)
            keywords=(dt = dt,)
        end
    else
        keywords=(t = t,)
        if typeof(t) <: AbstractRange
            t = Vector{RheoFloat}(t)
        end
    end
    σr = hasstress(d) ? Spline1D(d.t,d.σ)(t) : d.σ
    ϵr = hasstrain(d) ? Spline1D(d.t,d.ϵ)(t) : d.ϵ

    log = logadd_process(d, :resample, keywords=keywords, comment="Resample the data" ) 

    return RheoTimeData(σr, ϵr, t, log)

end


"""
    indexweight(self::RheoTimeData; elperiods::Vector{K}, time_boundaries::Union{Nothing, Vector{T}} = nothing, includelast=true) where {K<:Integer,T<:Real}

This function returns array indices (i.e. an array of integers) which can be sent to the `modelfit`, `modelstepfit` or `modeldiffeqfit` functions
to provide a weighted fitting whilst maintaining constant sample-rate.

Note that `time_boundaries` must have one more entry than `elperiods` so that all sections of weighting are fully defined by beginning and end points.

`indexweight` can underweight indices or overweight them. If `elperiods` in a given boundary is negative, every `abs(n)` index will be used, where `n` is the `elperiod`
corresponding to a given boundary. If `n` is positive, then indicies will be duplicated `n` times such that they are given a higher weighting during the fitting procedure.
If number of elements per period (`elperiods`) is `1` or `-1` it returns the original indicies for that boundary, whilst `0` is not accepted as a valid argument for `elperiods`.

The last element may or may not be included. By default the last element is forced to be included
but this can be negated by providing the keyword argument `includelast=false`.
"""
function indexweight(self::RheoTimeData; elperiods::Vector{K}, time_boundaries::Union{Nothing, Vector{T}} = nothing, includelast=true) where {K<:Integer,T<:Real}
    # get element-wise boundaries
    if isnothing(time_boundaries)
        boundaries = [1,length(self.t)]
    else
        boundaries = closestindices(self.t, time_boundaries)
    end

    # assert correct function signature
    @assert length(elperiods)==length(boundaries)-1 "Number of different sample periods must be 1 less than boundaries provided"
    @assert (count(iszero, elperiods)==0) "Number of elements cannot be zero"

    # initialise indices array
    indices = Integer[]

    # loop through boundaries
    for i in 1:(length(boundaries)-1)
        # upsampling, starts at each element then intepolates up N times
        if !signbit(elperiods[i])
            for j in boundaries[i]:(boundaries[i+1]-1)
                indices = vcat(indices, j*ones(Integer, elperiods[i]))
            end

        # downsampling, simply takes every N element as in downsample function
        elseif signbit(elperiods[i])
            append!(indices, collect(boundaries[i]:abs(elperiods[i]):(boundaries[i+1]-1)))
        end
    end

    # last element may be missed, so check and add if needed
    lastel = boundaries[end]
    numlastel = count(i->i==lastel, indices)
    countlastel = elperiods[end]
    if includelast && !signbit(countlastel) && numlastel<countlastel
        # up-sample case
        indices = vcat(indices, lastel*ones(Integer, countlastel-numlastel))

    elseif includelast && signbit(countlastel) && indices[end]<lastel
        append!(indices, lastel)
    end

    return indices
end

"""
    cutting(self::RheoTimeData, time_on::Real, time_off::Real)

Remove the data outside a specified time interval.

By specifing a time interval (`time_on`, `time_off`), a new `RheoTimeData` is returned without the data
lying outside time interval.
"""
function cutting(self::RheoTimeData, time_on::Real, time_off::Real)

    @assert isreal(time_on) && isreal(time_off) "Boundaries cannot be complex numbers"

    boundary_on = closestindex(self.t, time_on);
    boundary_off = closestindex(self.t, time_off);
    time = self.t[boundary_on:boundary_off]

    check = rheotimedatatype(self)
    if (check == strain_only)
        ϵ = self.ϵ[boundary_on:boundary_off]
        σ = [];
    elseif (check == stress_only)
        σ = self.σ[boundary_on:boundary_off]
        ϵ = [];
    elseif (check == strain_and_stress)
        ϵ = self.ϵ[boundary_on:boundary_off]
        σ = self.σ[boundary_on:boundary_off]
    end

    log = self.log === nothing ? nothing : [self.log; RheoLogItem( (type=:process, funct=:cutting, params=(time_on = time_on, time_off = time_off), keywords=() ),
                                        (comment="Cut section of the data between $time_on and $time_off",) ) ]

    return RheoTimeData(σ, ϵ, time, log)

end

"""
    smooth(self::RheoTimeData, τ::Real; pad::String="reflect")

Smooth data using a Gaussian Kernel to time scale `τ` (approximately half power).

Smooths both `σ` and `ϵ`. Sampling frequency must be constant as it is based on FFT convolution. Essentially a
low pass filter with frequencies of 1/τ being cut to approximately half power. For other pad types available
see ImageFiltering documentation. As of doc writing, pad options are: `"replicate"` (repeat edge values to
infinity), `"circular"` (image edges "wrap around"), `"symmetric"` (the image reflects relative to a position
between pixels), `"reflect"` (the image reflects relative to the edge itself).
"""
function smooth(self::RheoTimeData, τ::Real; pad::String="reflect")

    @eval import ImageFiltering: imfilter, Kernel

    # get sample-rate and Gaussian kernel (std. dev)
    samplerate = 1.0/getsampleperiod(self.t)
    Σ = getsigma(τ, samplerate)

    # smooth signal and return
    check = rheotimedatatype(self)
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

    log = logadd_process(self, :smooth, params=(τ = τ,), keywords=(pad = pad,), comment="Smooth data with timescale $τ" ) 

    RheoTimeData(sigma, epsilon, self.t, log)

end


function _smooth(x::Vector, window_len::Int=7, window::Symbol=:lanczos)
    w = getfield(Windows, window)(window_len)
    return filtfilt(w ./ sum(w), [1.0], x)
end

function smooth2(data, τ::Real; window::Symbol=:lanczos)

    # get sample-rate and Gaussian kernel (std. dev)
    # available from base.jl
    samplerate = 1.0/getsampleperiod(data.t)
    Σ = Int(floor(getsigma(τ, samplerate)))

    ϵ = hasstrain(data) ? _smooth(data.ϵ, Σ, window) : RheoFloat[] 
    σ = hasstress(data) ? _smooth(data.σ, Σ, window) : RheoFloat[] 

    log = data.log === nothing ? nothing : [data.log; RheoLogItem( (type=:process, funct=:smooth2, params=(τ = τ,), keywords=(window = window,) ),
                                        (comment="Smooth data with timescale $τ",) ) ]

    RheoTimeData(σ, ϵ, data.t, log)

end






"""
    extract(self::Union{RheoTimeData,RheoFreqData}, type::Union{TimeDataType,FreqDataType,Integer})

Extract specific fields form `RheoTimeData` or `RheoFreqData`.

Deprecated - prefer to use the functions `onlytime`, `onlystrain`, `onlystress`, and `onlyfreq`.

Extract can copy one or more fields from a given `RheoXData` variable into a new `RheoXData` one. The fields that
are copied are identified by the specified type of data.
If self is a `RheoTimeData`, the type that can be extracted is `time_only` (or `0`), `stress_only` (or `1`), `strain_only` (or `2`).
Note that `strain_and_stress` (or `3`) is not allowed.
If self is a `RheoFreqData`, the type that can be extracted is `freq_only` (or `0`).
"""
function extract(self::RheoTimeData, type::Union{TimeDataType,Integer})

        type = typeof(type)==Int ? TimeDataType(type) : type
        @assert (typeof(type)==TimeDataType) || (typeof(type)==Int) "Cannot extract frequency data from RheoTimeData"
        @assert (type!= strain_and_stress) && (type!= invalid) "Cannot extract both stress and strain"
        @assert (type!= invalid) "Cannot extract information from invalid time data"
        check = rheotimedatatype(self)


        log = logadd_process(self, :extract, params=(type=type,), comment="Data field extraction, $type from $check") 


        if type == time_only
            @assert check!= invalid "Time not available"
            return RheoTimeData([], [], self.t,log)
        elseif type == strain_only
            @assert (check == strain_and_stress) || (check == strain_only) "Strain not available"
            return RheoTimeData([], self.ϵ, self.t,log)
        elseif type == stress_only
            @assert (check == strain_and_stress) || (check == stress_only) "Stress not available"
            return RheoTimeData(self.σ, [], self.t,log)
        end

end

function extract(self::RheoFreqData, type::Union{FreqDataType,Integer})

        type = typeof(type)==Int ? FreqDataType(type) : type
        @assert (type!= with_modulus) "Cannot extract frequency with moduli"
        @assert (type!= invalid_freq_data) "Cannot extract information from invalid frequency data"
        check = rheofreqdatatype(self)
        @assert (check == with_modulus) "Frequency and modulii required"

        log = logadd_process(self, :extract, params=(type=type,), comment="Frequency field extraction")  
        return RheoFreqData([], [], self.ω,log)
end



"""
    onlytime(d::RheoTimeData)

Return a `RheoTimeData` object with only the timeline of the parameter `d`.

"""
function onlytime(d::RheoTimeData)
  extract(d,time_only)       
end



"""
    onlystress(d::RheoTimeData)

Return a `RheoTimeData` object with only the time and stress data of the parameter `d`.

"""
function onlystress(d::RheoTimeData)
  extract(d,stress_only)       
end



"""
    onlystrain(d::RheoTimeData)

Return a `RheoTimeData` object with only the time and strain data of the parameter `d`.

"""
function onlystrain(d::RheoTimeData)
  extract(d,strain_only)       
end



"""
    onlyfreq(d::RheoFreqData)

Return a `RheoFreqData` object with only the frequency data of the parameter `d`.

"""
function onlyfreq(d::RheoFreqData)
  extract(d,freq_only)       
end





#=
--------------------------------
Fitting and predicting functions
--------------------------------
=#
function fill_init_params(model::RheoModelClass, p0::Union{NamedTuple,Nothing})

    springpotnumber = length(findall(x->(x==:a)||(x==:β), model.freeparams))

    if isnothing(p0)
        # fill all initial guesses as 0.5, if two springpots then
        # make sure dissipation coefficients conform to constraints
        p0a = ones(RheoFloat, length(model.freeparams))/2
        if springpotnumber==2
            index = findfirst(model.freeparams.==:a)
            p0a[index] = RheoFloat(0.8)
        end
        @warn "Initial values for model parameters are set to $p0a by default"

    elseif length(p0)<length(model.freeparams)
        # check all provided symbols are part of the model, fill all parameters not provided
        for i in keys(p0)
            @assert (i in model.freeparams) "Incorrect parameter name provided for model initial conditions:" * i
        end
        # propagate list
        p0a = Vector{RheoFloat}(undef, length(model.freeparams))
        for (i,p) in enumerate(model.freeparams)
            p0a[i] = (p in keys(p0)) ? p0[p] : RheoFloat(1/2)
        end
        # check for second springpot specification
        indexa = findfirst(model.freeparams.==:a)
        indexβ = findfirst(model.freeparams.==:β)
        if !(:a in keys(p0)) && !(:β in keys(p0))
            p0a[indexa] = RheoFloat(0.8)
        elseif !(:β in keys(p0))
            p0a[indexβ] = p0[:a]/2
        elseif !(:a in keys(p0))
            p0a[indexa] = minimum((p0[:β]*2, RheoFloat(0.999)))
        end

        @warn "Unspecified initial guesses filled by default values: $p0a"

    else
        p0a = check_and_reorder_parameters(p0, model.freeparams,  err_string="initial guess")
    end

    @assert model._constraint(p0a)  "Initial guess not feasible"

    return p0a
end

function fill_lower_bounds(model::RheoModelClass, lo::Union{NamedTuple,Nothing})
    if isnothing(lo)
        loa = nothing

    elseif length(lo)<length(model.freeparams)
        loa = Vector{RheoFloat}(undef, length(model.freeparams))
        for (i,p) in enumerate(model.freeparams)
            loa[i] = (p in keys(lo)) ? lo[p] : zero(RheoFloat)
        end

    else
        loa = check_and_reorder_parameters(lo, model.freeparams,  err_string="low bounds")
    end

    return loa
end

function fill_upper_bounds(model::RheoModelClass, hi::Union{NamedTuple,Nothing})
    if isnothing(hi)
        hia = nothing

    elseif length(hi)<length(model.freeparams)
        hia = Vector{RheoFloat}(undef, length(model.freeparams))
        for (i,p) in enumerate(model.freeparams)
            hia[i] = (p in keys(hi)) ? hi[p] : RheoFloat(Inf)
        end

    else
        hia = check_and_reorder_parameters(hi, model.freeparams, err_string="high bounds")
    end

    return hia
end

"""
    modelfit(data::RheoTimeData, model::RheoModelClass, modloading::Symbol; p0::Union{NamedTuple,Nothing} = nothing, lo::Union{NamedTuple,Nothing} = nothing, hi::Union{NamedTuple,Nothing} = nothing, verbose::Bool = false, rel_tol = 1e-4, diff_method="BD", weights::Union{Vector{Integer},Nothing}=nothing)

Fit `RheologyData` struct to model and return a fitted model as a `RheologyModel` object.
For the fitting process RHEOS relies on the optimistion package [NLopt](https://nlopt.readthedocs.io/en/latest/).
By default, RHEOS makes use of a local derivative free algorithm, specifically the Tom Rowan's "Subplex"
It is possible to specify which algorithm NLopt should use with the keyword parameter `optmethod`, by providing the relevant symbols as defined by NLopt.jl.
Suitable options include `:LN_SBPLX` (default),  `:LN_COBYLA`, `:LN_BOBYQA` for local derivative free optimisation. Global optimisation methods are also available, but require all parameters to have lower and upper bounds set. More information about these algorithms is available on the [NLopt website](https://nlopt.readthedocs.io/en/latest/NLopt_Algorithms/).


# Arguments

- `data`: `RheoTimeData` struct containing all data
- `model`: `RheoModelClass` containing moduli functions and named tuple parameters
- `modloading`: `strain_imposed` or `1`, `stress_imposed` or `2`
- `p0`: Initial parameters to use in fit (uses 0.5 for all parameters if not defined), provided as a named tupple
- `lo`: Lower bounds for parameters, provided as a named tupple
- `hi`: Upper bounds for parameters, provided as a named tupple
- `verbose`: If true, prints parameters on each optimisation iteration
- `rel_tol`: Relative tolerance of optimization, see NLOpt docs for more details
- `diff_method`: Set finite difference formula to use for derivative, currently `"BD"` or `"CD"`
- `weights`: Vector of indices weighted according to importance, can be generated by `indexweight` function
- `optmethod`: optimisation algorithm used by NLOpt (passed as symbol or string) 
"""
function modelfit(data::RheoTimeData,
                  model::RheoModelClass,
                  modloading::LoadingType;
                  p0::Union{NamedTuple,Nothing,Dict} = nothing,
                  lo::Union{NamedTuple,Nothing,Dict} = nothing,
                  hi::Union{NamedTuple,Nothing,Dict} = nothing,
                  verbose::Bool = false,
                  rel_tol = 1e-4,
                  diff_method="BD",
                  weights::Union{Vector{Integer},Nothing,Vector{Int64}} = nothing,
                  optmethod::Union{Symbol,String}= :LN_SBPLX)

    p0a = fill_init_params(model, symbol_to_unicode(p0))
    loa = fill_lower_bounds(model, symbol_to_unicode(lo))
    hia = fill_upper_bounds(model, symbol_to_unicode(hi))

    rel_tol = convert(RheoFloat,rel_tol)

    check = rheotimedatatype(data)
    @assert (check == strain_and_stress) "Both stress and strain are required"

    # check provided weights are all valid
    if !isnothing(weights)
        @assert isempty(weights[weights.<1]) "Invalid weighting indices provided"
    end

    # use correct method for derivative
    if diff_method=="BD"
        deriv = derivBD
    elseif diff_method=="CD"
        deriv = derivCD
    end

    # get modulus function and derivative
    if modloading == stress_imposed
        dcontrolled = deriv(data.σ, data.t)
        measured = data.ϵ
        modulus = model._Ja
        modsing = (t->model._J(t,p0a))
        modused = "J"
    elseif modloading == strain_imposed
        dcontrolled = deriv(data.ϵ, data.t)
        measured = data.σ
        modulus = model._Ga
        modsing = (t->model._G(t,p0a))
        modused ="G"
    end

    # check necessary modulus is defined
    @assert modulusexists(modsing) "Modulus to use not defined"
    # get singularity presence
    sing = singularitytest(modsing)
    # get time step (only needed for convolution, which requires constant dt so t[2]-t[1] is sufficient)
    dt = data.t[2] - data.t[1]

    # time must start at 0 for convolution to work properly!
    t_zeroed = data.t .- minimum(data.t)

    # fit
    is_constant = constantcheck(data.t)

    # indices weighting only used for constant sample-rate data
    if !is_constant && !isnothing(weights)
        @warn "Indices weighting not used as variable sample-rate data has been provided"
    end

    (minf, minx, ret), timetaken, bytes, gctime, memalloc =
                    @timed leastsquares_init(   p0a,
                                                loa,
                                                hia,
                                                modulus,
                                                t_zeroed,
                                                dt,
                                                dcontrolled,
                                                measured;
                                                insight = verbose,
                                                constant_sampling = is_constant,
                                                singularity = sing,
                                                _rel_tol = rel_tol,
                                                indweights = weights,
                                                optmethod = Symbol(optmethod))

    println("Time: $timetaken s, Why: $ret, Parameters: $minx, Error: $minf")

    nt = NamedTuple{Tuple(model.freeparams)}(minx)

    if data.log !== nothing
        # Preparation of data for log item
        info=(comment="Fiting rheological model to data", model_name=model.name, model_params=nt, time_taken=timetaken, stop_reason=ret, error=minf)
        params=(model=model, modloading=modloading)
        keywords=(p0=p0, lo=lo, hi=hi, rel_tol=rel_tol, diff_method=diff_method)
        # Add data to the log
        push!(data.log, RheoLogItem( (type=:analysis, funct=:modelfit, params=params, keywords=keywords), info))
    end

    return RheoModel(model, nt);

end

"""
    modelpredict(data::RheoTimeData, model::RheoModel; diff_method="BD")

Given an incomplete data set (only either stress or strain missing) and model with values substituted into
parameters (`RheoModel`), return a new dataset based on the model.
If data is type of `stress_only`, then creep modulus (`:J`) is used; if data type is `strain_only` relaxation modulus (`:G`).
A complete `RheoTimeData` of type `strain_and_stress` is returned.
`diff_method` sets finite difference for calculating the derivative used in the hereditary integral and
can be either backwards difference (`"BD"`) or central difference (`"CD"`).
"""
function modelpredict(data::RheoTimeData, model::RheoModel; diff_method="BD")

    # use correct method for derivative
    if diff_method=="BD"
        deriv = derivBD
    elseif diff_method=="CD"
        deriv = derivCD
    end

    check = rheotimedatatype(data)
    @assert (check == strain_only)||(check == stress_only) "Need either strain only or stress only data. Data provide: " * string(check)

    if (check == strain_only)
        modulus = model._Ga
        modsing = model._G
        dcontrolled = deriv(data.ϵ, data.t)
    elseif (check == stress_only)
        modulus = model._Ja
        modsing = model._J
        dcontrolled = deriv(data.σ, data.t)
    end

    # get singularity presence
    sing = singularitytest(modsing)

    # get time step (only needed for convolution, which requires constant so t[2]-t[1] is sufficient)
    dt = data.t[2] - data.t[1]

    # time must start at 0 for convolution to work properly
    t_zeroed = data.t .- minimum(data.t)

    # get convolution
    if !sing && constantcheck(data.t)
        convolved = boltzconvolve(modulus, t_zeroed, dt, dcontrolled)

    elseif sing && constantcheck(data.t)
        t_zeroed[1] = 0.0 + (t_zeroed[2] - t_zeroed[1])/singularity_offset
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

    log = data.log === nothing ? nothing : [ data.log;
            RheoLogItem( (type=:process, funct=:modelpredict, params=(model::RheoModel,), keywords=(diff_method = diff_method,)),
                         (comment="Predicted data - modulus: $pred_mod, parameters:$(model.fixedparams)",) ) ]

    return RheoTimeData(sigma, epsilon, data.t, log)

end

"""
    modelstepfit(data::RheoTimeData, model::RheoModelClass, modloading::Union{LoadingType,Integer}; step=nothing, p0::Union{NamedTuple,Nothing} = nothing, lo::Union{NamedTuple,Nothing} = nothing, hi::Union{NamedTuple,Nothing} = nothing, verbose::Bool = false, rel_tol = 1e-4, weights::Union{Vector{Integer},Nothing} = nothing)

Same as `modelfit` except assumes a step loading. If this assumption is appropriate for the data
then fitting can be sped up greatly by use of this function. If modloading is `strain_imposed`, relaxation modulus is used,
then the element in the middle of the strain is assumed to be the amplitude of the step. If modloading is `stress_imposed`,
the creep modulus is used, then the middle element of the stress is assumed to be the amplitude of the step.
Alternatively, it is possible to define the value of the step by defining the optional `step` parameter.

# Arguments

- `data`: `RheoTimeData` struct containing all data
- `model`: `RheoModelClass` containing moduli and parameters tuples
- `modloading`: `strain_imposed` for relaxation modulus, `stress_imposed` for creep modulus
- `step`: Optional amplitude for step
- `p0`: Named tuple of initial parameters to use in fit (uses 0.5 for all parameters if none given)
- `lo`: Named tuple of lower bounds for parameters
- `hi`: Named tuple of upper bounds for parameters
- `verbose`: If true, prints parameters on each optimisation iteration
- `rel_tol`: Relative tolerance of optimization, see NLOpt docs for more details
- `weights`: Vector of indices weighted according to importance, can be generated by `indexweight` function
"""
function modelstepfit(data::RheoTimeData,
                  model::RheoModelClass,
                  modloading::Union{LoadingType,Integer};
                  step = nothing,
                  p0::Union{NamedTuple,Nothing} = nothing,
                  lo::Union{NamedTuple,Nothing} = nothing,
                  hi::Union{NamedTuple,Nothing} = nothing,
                  verbose::Bool = false,
                  rel_tol = 1e-4,
                  weights::Union{Vector{Integer},Nothing} = nothing)

    p0a = fill_init_params(model, symbol_to_unicode(p0))
    loa = fill_lower_bounds(model, symbol_to_unicode(lo))
    hia = fill_upper_bounds(model, symbol_to_unicode(hi))

    rel_tol = convert(RheoFloat, rel_tol)

    # get data type provided
    check = rheotimedatatype(data)
    # if step amplitude not provided, use loading data specified
    if isnothing(step)
        @assert (check == strain_and_stress) "Both stress and strain are required if step amplitude not provided."
        if (modloading == stress_imposed)
            # step amplitude is set to the middle value of the 'loading' data array
            controlled = data.σ[round(Integer, length(data.σ)/2)]
            measured = data.ϵ
            modulus = model._Ja
            modsing = (t->model._J(t, p0a))
            modused = "J"
        elseif (modloading == strain_imposed)
            # step amplitude is set to the middle value of the 'loading' data array
            controlled = data.ϵ[round(Integer, length(data.σ)/2)]
            measured = data.σ
            modulus = model._Ga
            modsing = (t->model._G(t, p0a))
            modused = "G"
        end
    # if step provided then can use that value
    elseif !isnothing(step)
        if (modloading == stress_imposed)
            @assert (check == strain_only)||(check == strain_and_stress) "Strain required"
            modulus = model._Ja
            modsing = (t->model._J(t, p0a))
            controlled = convert(RheoFloat, step);
            measured = data.ϵ
            modused = "J"
        elseif (modloading == strain_imposed)
            @assert (check == stress_only)||(check == strain_and_stress) "Stress required"
            measured = data.σ
            modulus = model._Ga
            modsing = (t->model._G(t, p0a))
            controlled =convert(RheoFloat, step);
            modused = "G"
        end
    end

    # step is assumed to be at 0, warn if data does not start at 0
    if data.t[1] != 0
        @warn "Step fitting assumes that step occurs at time 0.0s but data does not start at 0.0s. If this mismatch is unintended, this fitting may return erroneous results."
    end

    # check necessary modulus is defined
    @assert modulusexists(modsing) "Modulus to use not defined"
    # get singularity presence
    sing = singularitytest(modsing)

    # time must start at 0 for modulus to be defined properly
    t_zeroed = data.t .- minimum(data.t)

    # start fit
    (minf, minx, ret), timetaken, bytes, gctime, memalloc =
                        @timed leastsquares_stepinit(p0a,
                                                    loa,
                                                    hia,
                                                    modulus,
                                                    t_zeroed,
                                                    controlled,
                                                    measured;
                                                    insight = verbose,
                                                    singularity = sing,
                                                    _rel_tol = rel_tol,
                                                    indweights = weights)

    nt = NamedTuple{Tuple(model.freeparams)}(minx)

    if data.log !== nothing
        # Preparation of data for log item
        info=(comment="Fiting rheological model to data (step input assumed)", model_name=model.name, model_params=nt, time_taken=timetaken, stop_reason=ret, error=minf)
        params = (model = model, modloading = modloading)
        keywords = (step = step, p0 = p0, lo = lo, hi = hi, rel_tol = rel_tol)
        # Add data to the log
        push!(data.log, RheoLogItem( (type=:analysis, funct=:modelstepfit, params=params, keywords=keywords), info))
    end

    print("Time: $timetaken s, Why: $ret, Parameters: $minx, Error: $minf")

    return RheoModel(model, nt)
end

"""
    modelsteppredict(data::RheoTimeData, model::RheoModel; step_on::Real = 0.0)

Same as `modelpredict` but assumes a step loading with step starting at `step_on` or closest actual value to that
specified. If the loading data is variable, the magnitude in the middle of the array is used.
"""
function modelsteppredict(data::RheoTimeData, model::RheoModel; step_on::Real = 0.0)

    check = rheotimedatatype(data)
    @assert (check == strain_only)||(check == stress_only) "Need either strain only or stress only data. Data provide: " * string(check)

    if (check == strain_only)
        modulus = model._Ga
        modsing = model._G
        controlled = data.ϵ[round(Integer, length(data.ϵ)/2)]
    elseif (check == stress_only)
        modulus = model._Ja
        modsing = model._J
        controlled = data.σ[round(Integer, length(data.σ)/2)]
    end

    # get closest index and time to desired step-on time point
    stepon_el = closestindex(data.t, step_on)
    stepon_closest = data.t[stepon_el]

    # check singularity presence at time closest to step
    sing = singularitytest(modsing)

    # get predicted
    predicted = zeros(RheoFloat, length(data.t))
    if !sing
        predicted[stepon_el:end] = controlled*modulus(data.t[stepon_el:end] .- stepon_closest)
    elseif sing
        predicted[stepon_el:end] = controlled*modulus(data.t[stepon_el:end] .- stepon_closest)
        @warn "Singularity (Inf value) is present in modulus used for modelsteppredict"
    end

    if check == stress_only
        σ = data.σ
        ϵ = predicted
    elseif check == strain_only
        σ = predicted
        ϵ = data.ϵ
    end

    # store operation
    if (check == stress_only)
        sigma = data.σ
        epsilon = predicted
        pred_mod = "creep function J"
    elseif (check == strain_only)
        sigma = predicted
        epsilon = data.ϵ
        pred_mod = "relaxation function G"
    end
    time = data.t

    log = data.log === nothing ? nothing : [ data.log;
            RheoLogItem( (type=:process, funct=:modelsteppredict, params=(model::RheoModel,), keywords=(step_on = step_on,)),
                         (comment="Predicted step response - modulus: $pred_mod, parameters:$(model.fixedparams)",) ) ]


    return RheoTimeData(sigma,epsilon,time, log)
end

function obj_dynamic(params, grad, ω, dataGp, dataGpp, modelGp, modelGpp; _insight::Bool = false)

    if _insight
        println("Current Parameters: ", params)
    end

    costGp = sum((dataGp - modelGp(ω, params)).^2)
    costGpp = sum((dataGpp - modelGpp(ω, params)).^2)

    cost = costGp + costGpp
end

function obj_dynamic_mean(params, grad, ω, dataGp, dataGpp, modelGp, modelGpp, meanGp, meanGpp; _insight::Bool = false)

    if _insight
        println("Current Parameters: ", params)
    end

    costGp = sum((dataGp/meanGp - modelGp(ω, params)/meanGp).^2)
    costGpp = sum((dataGpp/meanGpp - modelGpp(ω, params)/meanGpp).^2)

    cost = costGp + costGpp

end

function obj_dynamic_log(params, grad, ω, dataGp, dataGpp, modelGp, modelGpp; _insight::Bool = false)

    if _insight
        println("Current Parameters: ", params)
    end

    costGp = sum((log.(dataGp) - log.(modelGp(ω, params))).^2)
    costGpp = sum((log.(dataGpp) - log.(modelGpp(ω, params))).^2)

    cost = costGp + costGpp

end

function obj_dynamic_local(params, grad, ω, dataGp, dataGpp, modelGp, modelGpp; _insight::Bool = false)

    if _insight
        println("Current Parameters: ", params)
    end

    costGp = sum(((dataGp - modelGp(ω, params))./dataGp).^2)
    costGpp = sum(((dataGpp - modelGpp(ω, params))./dataGpp).^2)

    cost = costGp + costGpp

end

function obj_dynamic_manual(params, grad, ω, dataGp, dataGpp, modelGp, modelGpp, weights; _insight::Bool = false)

    if _insight
        println("Current Parameters: ", params)
    end

    costGp = sum((dataGp - modelGp(ω, params)).^2)
    costGpp = sum((dataGpp - modelGpp(ω, params)).^2)

    cost = weights[1]*costGp + weights[2]*costGpp
end

"""
    dynamicmodelfit(data::RheoFreqData, model::RheoModelClass; p0::Union{NamedTuple,Nothing} = nothing, lo::Union{NamedTuple,Nothing} = nothing, hi::Union{NamedTuple,Nothing} = nothing, verbose::Bool = false, rel_tol = 1e-4) where T<:Real

Fits model to the frequency/loss+storage moduli data.

All arguments are as described below.
As this fitting procedure is fitting two functions simultaneously (the storage
and loss moduli), if left untransformed the fit would tend to favour the
modulus which is larger in magnitude and not fit the other modulus well. To avoid this,
RHEOS offers a number of data transforms which can be used by changing `weights` argument.

# Arguments

- `data`: `RheoFreqData` struct containing all data
- `model`: `RheoModelClass` containing moduli and symbols of parameters
- `p0`: Initial parameters to use in fit (uses 0.5 for all parameters if none given)
- `lo`: Lower bounds for parameters
- `hi`: Upper bounds for parameters
- `verbose`: If true, prints parameters on each optimisation iteration
- `rel_tol`: Relative tolerance of optimization, see NLOpt docs for more details
- `weights`: Weighting mode for storage and loss modulus (`"none"`, `"mean"`, `"log"`, `"local"` or manually specified)
"""
function dynamicmodelfit(data::RheoFreqData,
                model::RheoModelClass;
                p0::Union{NamedTuple,Nothing} = nothing,
                lo::Union{NamedTuple,Nothing} = nothing,
                hi::Union{NamedTuple,Nothing} = nothing,
                verbose::Bool = false,
                rel_tol::T = 1e-4,
                weights::Union{String, Vector{T}}="local") where T<:Real

    p0a = fill_init_params(model, symbol_to_unicode(p0))
    loa = fill_lower_bounds(model, symbol_to_unicode(lo))
    hia = fill_upper_bounds(model, symbol_to_unicode(hi))

    # check necessary moduli are defined
    modsingGp = (ω->model._Gp(ω,p0a))
    modsingGpp = (ω->model._Gpp(ω,p0a))
    @assert modulusexists(modsingGp) "Storage modulus not defined for this model"
    @assert modulusexists(modsingGpp) "Loss modulus not defined for this model"

    # initialise NLOpt.Opt object with :LN_SBPLX Subplex algorithm
    opt = Opt(:LN_SBPLX, length(p0a))

    if !isnothing(loa)
       lower_bounds!(opt, loa)
    end

    if !isnothing(hia)
       upper_bounds!(opt, hia)
    end

    # set relative tolerance
    rel_tol = convert(RheoFloat, rel_tol)
    xtol_rel!(opt, rel_tol)

    # set objective/cost function
    if weights=="none"
        min_objective!(opt, (params, grad) -> obj_dynamic(params, grad, data.ω, data.Gp, data.Gpp, model._Gpa, model._Gppa; _insight = verbose))

    elseif weights=="mean"
        meanGp = sum(data.Gp)/length(data.Gp)
        meanGpp = sum(data.Gp)/length(data.Gp)
        min_objective!(opt, (params, grad) -> obj_dynamic_mean(params, grad, data.ω, data.Gp, data.Gpp, model._Gpa, model._Gppa, meanGp, meanGpp; _insight = verbose))

    elseif weights=="log"
        @warn "Note that a logarithmic rescaling will fail if Gp or Gpp data contain 0.0 values as it will result in -Inf cost. Trying a different rescaling scheme, or not fitting around ω≈0.0 may alleviate the issue."
        min_objective!(opt, (params, grad) -> obj_dynamic_log(params, grad, data.ω, data.Gp, data.Gpp, model._Gpa, model._Gppa; _insight = verbose))

    elseif weights=="local"
        @warn "Note that a local rescaling will fail if Gp or Gpp data contain 0.0 values as it will result in division by 0.0. Trying a different rescaling scheme, or not fitting around ω≈0.0 may alleviate the issue."
        min_objective!(opt, (params, grad) -> obj_dynamic_local(params, grad, data.ω, data.Gp, data.Gpp, model._Gpa, model._Gppa; _insight = verbose))

    elseif typeof(weights)==Vector{T} && length(weights)==2
        min_objective!(opt, (params, grad) -> obj_dynamic_manual(params, grad, data.ω, data.Gp, data.Gpp, model._Gpa, model._Gppa, weights; _insight = verbose))

    end

    # timed fitting
    (minf, minx, ret), timetaken, bytes, gctime, memalloc = @timed NLopt.optimize(opt, p0a)

    nt = NamedTuple{Tuple(model.freeparams)}(minx)

    if data.log !== nothing
        # Preparation of data for log item
        info = (comment="Fiting rheological model to frequency spectrum", model_name=model.name, model_params=nt, time_taken=timetaken, stop_reason=ret, error=minf)
        params = (model=model, )
        keywords = (p0=p0, lo=lo, hi=hi, rel_tol=rel_tol, weights=weights)
        # Add data to the log
        push!(data.log, RheoLogItem( (type=:analysis, funct=:dynamicmodelfit, params=params, keywords=keywords), info))
    end

    print("Time: $timetaken s, Why: $ret, Parameters: $minx, Error: $minf")

    return RheoModel(model, nt)
end

"""
    dynamicmodelpredict(data::RheoFreqData, model::RheoModel)

Given dynamic rheology data with only frequency and model where parameters have been substituted.
Returns another `RheoFreqData` instance with the predicted `Gp` and `Gpp` based on the
frequencies and model given as arguments.
"""
function dynamicmodelpredict(data::RheoFreqData, model::RheoModel)

    predGp = model._Gpa(data.ω)
    predGpp = model._Gppa(data.ω)

    log = data.log === nothing ? nothing : [ data.log;
            RheoLogItem( (type=:process, funct=:dynamicmodelpredict, params=(model::RheoModel,), keywords=() ),
                         (comment="Calculated frequency spectrum - parameters:$(model.fixedparams)",) ) ]

    return RheoFreqData(predGp, predGpp, data.ω, log)

end

