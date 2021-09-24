#!/usr/bin/env julia

#=
-----------------
Utility Functions
-----------------
=#
"""
    trapz(y, x)

Array based trapezoidal integration of y with respect to x.

Limits of integration defined by the start and end points of the arrays.
"""
function trapz(y::Vector{RheoFloat}, x::Vector{RheoFloat})

    n = length(x)

    @assert n==length(y) "X and Y array length must match."
    n==1 && return zero(RheoFloat)

    r = zero(RheoFloat)
    @inbounds for i in 2:n
        r += (y[i-1] + y[i])*(x[i] - x[i-1])
    end

    r/2

end

"""
    derivCD(y, x)

Given two arrays of data, x and y, calculate dy/dx using central difference
method and forward and backward difference for array boundaries.
"""
function derivCD(y::Vector{RheoFloat}, x::Vector{RheoFloat})

    # get length
    N = length(x)

    # assert y and x arrays are same length
    @assert length(y)==N "X and Y Array lengths must match."

    # initialise zero array of length y
    ydot = similar(y)

    # assume 'imaginary' previous point is 0.0, and Δx is the same as the next one ahead
    # this is a physical assumption that material is at rest before first data point.
    # Could be problematic in some cases if sudden jump as we are actually missing
    # important information about how quickly that jump happened.
    @inbounds ydot[1] = y[1]/(x[2] - x[1])

    # central difference with uneven spacing for general case of constant or variable sample rate
    @inbounds for i in 2:(N-1)
        Δx₁ = x[i] - x[i-1]
        Δx₂ = x[i+1] - x[i]
        ydot[i] = (y[i+1]*Δx₁^2 + (Δx₂^2 - Δx₁^2)*y[i] - y[i-1]*Δx₂^2)/(Δx₁*Δx₂*(Δx₁ + Δx₂))
    end

    # 1st order backwards difference for last element
    ydot[N] = (y[N] - y[N-1])/(x[N] - x[N-1])

    return ydot

end

"""
    derivBD(y, x)

Given two arrays of data, x and y, calculate dy/dx using 1st order
backward difference. Assumes y==0 at a previous point, i.e.
y is 'at rest'. Captures instantaneous loading where derivCD will smooth.
"""
function derivBD(y::Vector{RheoFloat}, x::Vector{RheoFloat})

    # get length
    N = length(x)

    # assert y and x arrays are same length
    @assert length(y)==N "X and Y Array lengths must match."

    # initialise zero array of length y
    ydot = similar(y)

    # assume 'imaginary' previous point is 0.0, and Δx is the same as the next one ahead
    # this is a physical assumption that material is at rest before first data point.
    # Could be problematic in some cases if sudden jump as we are actually missing
    # important information about how quickly that jump happened.
    @inbounds ydot[1] = y[1]/(x[2] - x[1])

    # backwards difference method for rest of points
    @inbounds for i in 2:N
        ydot[i] = (y[i] - y[i-1])/(x[i] - x[i-1])
    end

    return ydot

end

"""
    doublederivCD(y, x)

Given two arrays of data, x and y, calculate d^2(y)/dx^2 using 2nd order difference. Data
should conform to the following constraints:
- No step loading (signal should start at 0)
"""
function doublederivCD(y::Vector{RheoFloat}, x::Vector{RheoFloat})
    # get length
    N = length(x)

    # assert y and x arrays are same length
    @assert length(y)==N "X and Y Array lengths must match."

    # initialise zero array of length y
    yddot = similar(y)

    # no physical assumptions made, unlike for first order BD derivative.
    # forward difference for first element, backwards difference for last element
    @inbounds yddot[1] = (y[3] - 2*y[2] + y[1])/(x[2] - x[1])^2
    @inbounds yddot[N] = (y[N] - 2*y[N-1] + y[N-2])/(x[N] - x[N-1])^2

    # central difference method for rest of points
    @inbounds for i in 2:(N-1)
        yddot[i] = (y[i+1] - 2*y[i] + y[i-1])/(x[i] - x[i-1])^2
    end

    return yddot
end

"""
    derivFbasic(α, ydot, yddot, t, dt)

Compute the Caputo fractional derivative using a simple convolution on the integrated-by-parts
form which does not have a singularity in the kernel. Sample rate must be constant, and this
attribute is not checked in the function body for computational efficiency.
"""
function derivFbasic(α, ydot, yddot, t, dt)
    I1 = ydot[1]*t.^(1-α)
    I2 = conv(t.^(1-α), yddot)[1:length(t)]*dt
    return (I1 + I2)/gamma(2 - α)
end

function constantcheck(t::Vector{RheoFloat})
    # get array of backward differences
    diff = t[2:end] - t[1:end-1]
    # check if any element is not approximately equal to 1st element
    check = all(x -> x ≈ diff[1], diff)
end

function getsampleperiod(t::Vector{RheoFloat})
    # check sample rate is constant, otherwise sample period varies
    @assert constantcheck(t) "Sample-rate must be constant"
    # return sample period
    sampleperiod = t[2] - t[1]
end

"""
    closestindex(x, val)

Find the index of the array element closest to val.
"""
function closestindex(x::Vector{T}, val::Real) where {T<:Real}

    # intialise closest match variable, assuming best match is index 1
    ibest = 1

    # diff between value and current element
    dxbest = abs(x[ibest]-val)

    # loop through all elements, looking for smallest difference
    @inbounds for I in eachindex(x)
        dx = abs(x[I]-val)
        if dx < dxbest
            dxbest = dx
            ibest = I
        end
    end

    ibest
end

"""
    closestindices(x, vals)

Uses `closestindex` iteratively to find closest index for each value in `vals` array,
returns array of indices.
"""
closestindices(x::Vector{T1}, vals::Vector{T2}) where {T1<:Real, T2<:Real} = broadcast(closestindex, (x,), vals)

#=
--------------------------------
Preprocessing base functionality
--------------------------------
=#
"""
    getsigma(τ, samplerate)

Generate sigma/std deviation for gaussian smoothing kernel.

Acts as a low pass filter. Information of time scale τ will be half power,
faster will be increasingly cut. Called by smooth function.
"""
function getsigma(τ::Real, samplerate::Real)

    # get freq, reduce to half power and generate gaussian std. deviation
    smoothfreq = 1.0/τ

    sF_halfpower = smoothfreq/sqrt(2.0*log(2.0))

    σ = samplerate/(2.0*π*sF_halfpower)

    RheoFloat(σ)
end


#=
-----------------------------------------
Fitting and predicting base functionality
-----------------------------------------
=#
function singularitytest(modulus; t1::RheoFloat=zero(RheoFloat))

    startval = modulus(t1)

    if isnan(startval) || startval == Inf
        return true
    else
        return false
    end

end

"""
    boltzintegral_nonsing(modulus, time_series, prescribed_dot)

Calculate Boltzmann Superposition integral using direct integration method.

This is much slower and slightly less accurate (depending on sample resolution)
than the convolution method. However, it works for variable sample rate.

# Arguments

- `modulus`: Viscoelastic modulus function in which the parameters have been substituted
- `time_series`: The array of times
- `prescribed_dot`: Derivative of (usually prescribed) variable inside the integration kernel
"""
function boltzintegral_nonsing(modulus, time_series::Vector{Float64}, prescribed_dot::Vector{Float64})

    # need to add an additional 'previous' time point to capture any instantaneous loading
    time_previous = time_series[1] - (time_series[2] - time_series[1])
    time_mod = vcat([time_previous], time_series)
    # material is assumed at rest at this 'previous' time
    prescribed_dot_mod = vcat([0.0], prescribed_dot)

    I = zeros(length(time_mod))
    @inbounds for (i,v) in enumerate(time_mod)
        # generate integral for each time step
        τ = time_mod[1:i]
        Modulus_arg = v .- τ
        Modulusᵢ = modulus(Modulus_arg)
        df_dtᵢ = prescribed_dot_mod[1:i]
        intergrand = Modulusᵢ.*df_dtᵢ

        I[i] = trapz(intergrand, τ)
    end
    # fix initial point
    I[2] = (prescribed_dot[1]*modulus(time_series)*(time_series[2] - time_series[1]))[1]
    I[2:end]

end

"""
    obj_var_nonsing(params, grad, modulus, time_series, prescribed_dot, measured; _insight::Bool = false)

Generate the sum-of-squares of the difference between `measured` and the Boltzmann
convolution integral of the selected `modulus` and `prescribed_dot`. Used when
sample rate is variable and no singularity in model.

# Arguments

- `params`: Array of parameters sent to `modulus`
- `grad`: Gradient argument used by NLOpt
- `modulus`: Viscoelastic modulus function
- `time_series`: Array of time data
- `prescribed_dot`: Convolved with `modulus` over `time_series`, usually dϵ/dt for stress relaxation and dσ/dt for creep
- `measured`: Data for comparison against, usually σ for stress relaxation and ϵ for creep
- `_insight`: Declare whether insight info should be shown when this function is called, true or false
"""
function obj_var_nonsing(params, grad, modulus, time_series, prescribed_dot, measured; _insight::Bool = false)

    _insight && println("Current Parameters: ", params)

    mod = (t->modulus(t,params))
    convolved = boltzintegral_nonsing(mod, time_series, prescribed_dot)

    cost = sum((measured - convolved).^2)

end

"""
    boltzintegral_sing(modulus, time_series, prescribed_dot)

Calculate Boltzmann Superposition integral using direct integration method.

This is much slower and slightly less accurate (depending on sample resolution)
than the convolution method. However, it works for variable sample rate.

# Arguments

- `modulus`: Viscoelastic modulus function in which the parameters have been substituted
- `time_series`: The array of times
- `prescribed_dot`: Derivative of (usually prescribed) variable inside the integration kernel
"""
function boltzintegral_sing(modulus, time_series::Vector{Float64}, prescribed_dot::Vector{Float64})

    offset_t0 = (time_series[2] - time_series[1])/singularity_offset
    # need to add an additional 'previous' time point to capture any instantaneous loading
    time_previous = time_series[1] - (time_series[2] - time_series[1])
    time_mod = vcat([time_previous], time_series)
    # material is assumed at rest at this 'previous' time
    prescribed_dot_mod = vcat([0.0], prescribed_dot)

    offset = offset_t0
    I = zeros(length(time_mod))
    for (i,v) in enumerate(time_mod)
        if i>1
            offset = (time_mod[i] - time_mod[i-1])/singularity_offset
        end
        τ = time_mod[1:i]
        Modulus_arg = v .- τ
        Modulus_arg[end] = offset
        Modulusᵢ = modulus(Modulus_arg)
        df_dtᵢ = prescribed_dot_mod[1:i]

        intergrand = Modulusᵢ.*df_dtᵢ

        I[i] = trapz(intergrand, τ)
    end
    # fix initial point
    I[2] = (prescribed_dot[1]*modulus([offset_t0])*(time_series[2] - time_series[1]))[1]
    I[2:end]

end

"""
    obj_var_sing(params, grad, modulus, time_series, prescribed_dot, measured; _insight::Bool = false)

Generate the sum-of-squares of the difference between `measured` and the Boltzmann
convolution integral of the selected `modulus` and `prescribed_dot`. Used when
sample rate is variable and there IS singularity in model.

# Arguments

- `params`: Array of parameters sent to `modulus`
- `grad`: Gradient argument used by NLOpt
- `modulus`: Viscoelastic modulus function
- `time_series`: Array of time data
- `prescribed_dot`: Convolved with `modulus` over `time_series`, usually dϵ/dt for stress relaxation and dσ/dt for creep
- `measured`: Data for comparison against, usually σ for stress relaxation and ϵ for creep
- `_insight`: Declare whether insight info should be shown when this function is called, true or false
"""
function obj_var_sing(params, grad, modulus, time_series, prescribed_dot, measured; _insight::Bool = false)

    _insight && println("Current Parameters: ", params)

    mod = (t->modulus(t,params))
    convolved = boltzintegral_sing(mod, time_series, prescribed_dot)

    # skip first point for error computation, workaround for singularity approximation error
    cost = sum((measured[2:end] - convolved[2:end]).^2)

end

"""
    boltzconvolve(modulus, time_series, dt::RheoFloat,prescribed_dot)

Calculate Boltzmann Superposition integral using convolution method.

This is much faster and slightly more accurate (depending on sample resolution)
than the integral method. However, it only works for constant sample rate.

# Arguments

- `modulus`: Viscoelastic modulus function in which the parameters have been substituted
- `time_series`: The array of times
- `dt`: Constant time step (sample period)
- `prescribed_dot`: Derivative of (usually prescribed) variable inside the integration kernel
"""
function boltzconvolve(modulus, time_series, dt, prescribed_dot)

    Modulus = modulus(time_series)
    # fast convolution
    β = conv(Modulus, prescribed_dot)
    # pick out relevant elements (1st half) and multiply by dt
    β = β[1:length(time_series)]*dt

end

"""
    obj_const_nonsing(params, grad, modulus, time_series, dt, prescribed_dot, measured; _insight::Bool = false)

Generate the sum-of-squares of the difference between `measured` and the Boltzmann
convolution integral of the selected `modulus` and `prescribed_dot`. Used when
sample rate is constant and model does not feature singularity.

# Arguments

- `params`: Array of parameters sent to `modulus`
- `grad`: Gradient argument used by NLOpt
- `modulus`: Viscoelastic modulus function
- `time_series`: Array of time data
- `dt`: Constant time step (sample period)
- `prescribed_dot`: Convolved with `modulus` over `time_series`, usually dϵ/dt for stress relaxation and dσ/dt for creep
- `measured`: Data for comparison against, usually σ for stress relaxation and ϵ for creep
- `_insight`: Declare whether insight info should be shown when this function is called, true or false
"""
function obj_const_nonsing(params, grad, modulus, time_series, dt, prescribed_dot, measured; _insight::Bool = false)

    _insight && println("Current Parameters: ", params)

    # Replace parameters inside the viscoelastic modulus
    mod = (t->modulus(t,params))
    convolved = boltzconvolve(mod, time_series, dt, prescribed_dot)

    cost = sum((measured - convolved).^2)
    return cost

end

"""
    obj_const_sing(params, grad, modulus, time_series, dt, prescribed_dot, measured; _insight::Bool = false)

Generate the sum-of-squares of the difference between `measured` and the Boltzmann
convolution integral of the selected `modulus` and `prescribed_dot`. Used when
sample rate is constant and model does feature singularity.

# Arguments

- `params`: Array of parameters sent to `modulus`
- `grad`: Gradient argument used by NLOpt
- `modulus`: Viscoelastic modulus function
- `time_series`: Array of time data
- `dt`: Constant time step (sample period)
- `prescribed_dot`: Convolved with `modulus` over `time_series`, usually dϵ/dt for stress relaxation and dσ/dt for creep
- `measured`: Data for comparison against, usually σ for stress relaxation and ϵ for creep
- `_insight`: Declare whether insight info should be shown when this function is called, true or false
"""
function obj_const_sing(params, grad,modulus, time_series,dt, prescribed_dot, measured; _insight::Bool = false)

    _insight && println("Current Parameters: ", params)

    mod = (t->modulus(t,params))
    # convolved = boltzconvolve_sing(modulus, time_series, dt, params, prescribed_dot)
    convolved = boltzconvolve(mod, time_series, dt, prescribed_dot)

    # skip first point for error computation, workaround for singularity approximation error
    cost = sum((measured[2:end] - convolved[2:end]).^2)

end

"""
    obj_const_weighted(params, grad, modulus, time_series, dt, prescribed_dot, measured_weighted, weights; _insight::Bool = false)

Post-convolution weighted constant sample-rate cost function for emphasising regions without introducing variable sample-rates.
"""
function obj_const_weighted(params, grad, modulus, time_series, dt, prescribed_dot, measured_weighted, weights; _insight::Bool = false)

    _insight && println("Current Parameters: ", params)

    # Replace parameters inside the viscoelastic modulus
    mod = (t->modulus(t,params))
    convolved = boltzconvolve(mod, time_series, dt, prescribed_dot)

    cost = sum((measured_weighted - convolved[weights]).^2)
    return cost
end

"""
    leastsquares_init(params_init, low_bounds, hi_bounds, modulus, time_series, dt, prescribed_dot, measured; insight = false, sampling = "constant", singularity = false)

Initialise then begin a least squares fitting of the supplied data.

# Arguments

- `params_init`: Initial parameters to be used (starting guess)
- `low_bounds`: Lower bounds for parameters
- `hi_bounds`: Higher bounds for parameters
- `modulus`: Viscoelastic modulus function
- `time_series`: Array of time data
- `dt`: Constant time step (sample period)
- `prescribed_dot`: Convolved with `modulus` over `time_series`, usually dϵ/dt for stress relaxation and dσ/dt for creep
- `measured`: Data for comparison against, usually σ for stress relaxation and ϵ for creep
- `constant_sampling`: true if sample rate is constant, false otherwise
- `insight`: Declare whether insight info should be shown when this function is called, true or false
- `singularity`: Presence of singularity in model
- `indweight`: indices weighting, see `modelfit` for more discussion
- `optmethod`: optimisation algorithm used by NLOpt. 
- `opttimeout`: allows user to set a wall clock timeout on optimisation 
"""
function leastsquares_init(params_init::Vector{RheoFloat},
                            low_bounds::RheovecOrNone,
                            hi_bounds::RheovecOrNone, 
                            modulus,
                            time_series::Vector{RheoFloat},
                            dt::RheoFloat,
                            prescribed_dot::Vector{RheoFloat},
                            measured::Vector{RheoFloat};
                            insight::Bool = false,
                            constant_sampling::Bool=true,
                            singularity::Bool = false,
                            rel_tol_x::Union{Real,Nothing} = nothing,
                            rel_tol_f::Union{Real,Nothing} = nothing,
                            indweights = nothing,
                            optmethod::Symbol = :LN_SBPLX,
                            opttimeout::Union{Real,Nothing} = nothing,
                            optmaxeval::Union{Integer,Nothing} = nothing)
                           

    # initialise NLOpt.Opt object with :LN_SBPLX Subplex algorithm
    opt = Opt(optmethod, length(params_init))
    # opt = Opt(:LN_BOBYQA, length(params_init))    # Passing tests
    # opt = Opt(:LN_COBYLA, length(params_init))    # Failing test - not precise enough?

    # set optimiser stopping criteria

    # wall clock timeout
    if !isnothing(opttimeout)
        opttimeout = convert(RheoFloat, opttimeout)
        maxtime!(opt, opttimeout)
    end

    # evaluation cycle ceiling
    if !isnothing(optmaxeval)
        maxeval!(opt, optmaxeval)
    end

    # input parameter change tolerance
    if !isnothing(rel_tol_x)
        rel_tol = convert(RheoFloat, rel_tol_x)
        xtol_rel!(opt, rel_tol_x)
    end

    # objective function change tolerance 
    if !isnothing(rel_tol_f)
        rel_tol_f = convert(RheoFloat, rel_tol_f)
        ftol_rel!(opt, rel_tol_f)
    end

    # set lower bounds and upper bounds unless they take null value
    if !isnothing(low_bounds)
        low_bounds = convert(Vector{Float64},low_bounds)
        lower_bounds!(opt, low_bounds)
    end

    if !isnothing(hi_bounds)
        hi_bounds = convert(Vector{Float64}, hi_bounds)
        upper_bounds!(opt, hi_bounds)
    end

    # Convert to float64 to avoid conversion by NLOpt
    params_init = convert(Vector{Float64},params_init)
    time_series = convert(Vector{Float64},time_series)
    prescribed_dot = convert(Vector{Float64},prescribed_dot)
    measured = convert(Vector{Float64},measured)
    dt = convert(Float64, dt)

    # set Opt object as a minimisation objective. Use a closure for additional
    # arguments sent to object objectivefunc
    if constant_sampling && !singularity && isnothing(indweights)
        min_objective!(opt, (params, grad) -> obj_const_nonsing(params, grad, modulus,
                                                            time_series, dt,
                                                            prescribed_dot, measured;
                                                            _insight = insight))

    elseif constant_sampling && singularity && isnothing(indweights)
        # remove singularity, just go close to it, 1/10th over first sample period
        time_series[1] = 0.0 + (time_series[2] - time_series[1])/singularity_offset
        min_objective!(opt, (params, grad) -> obj_const_sing(params, grad, modulus,
                                                        time_series, dt,
                                                        prescribed_dot, measured;
                                                        _insight = insight))

    elseif constant_sampling && !singularity && !isnothing(indweights)
        measured_weighted = measured[indweights]
        min_objective!(opt, (params, grad) -> obj_const_weighted(params, grad, modulus,
                                                            time_series, dt,
                                                            prescribed_dot, measured_weighted, indweights;
                                                            _insight = insight))

    elseif constant_sampling && singularity && !isnothing(indweights)
        indsingular = indweights[indweights.>1]
        measured_weighted = measured[indsingular]
        # remove singularity, just go close to it, 1/10th over first sample period
        time_series[1] = 0.0 + (time_series[2] - time_series[1])/singularity_offset
        min_objective!(opt, (params, grad) -> obj_const_weighted(params, grad, modulus,
                                                        time_series, dt,
                                                        prescribed_dot, measured_weighted, indsingular;
                                                        _insight = insight))

    elseif !constant_sampling && !singularity
        min_objective!(opt, (params, grad) -> obj_var_nonsing(params, grad, modulus,
                                                        time_series, prescribed_dot,
                                                        measured; _insight = insight))

    elseif !constant_sampling && singularity
        @warn "Note that large changes in sample rate, particularly from low sample rate to
        much higher sample rate can introduce significant innacuracies due to the singularity
        approximation method being used. To avoid this, please resample your data to a constant
        sample rate before fitting."
        min_objective!(opt, (params, grad) -> obj_var_sing(params, grad, modulus,
                                                        time_series, prescribed_dot,
                                                        measured; _insight = insight))

    end

    # minimise objective func, minx are the parameters resulting in minimum
    (minf, minx, ret) = NLopt.optimize(opt, params_init)

    # return all
    return (convert(RheoFloat,minf), convert(Vector{RheoFloat},minx), ret)

end

#=
----------------------------------------------
Step fitting and predicting base functionality
----------------------------------------------
=#
function obj_step_nonsing(params, grad, modulus, t, prescribed::Float64, measured::Vector{Float64}; _insight=false)
    if _insight
        println("Current Parameters: ", params)
    end

    mod = (t->modulus(t,params))
    estimated = prescribed*mod(t)

    cost = sum((measured - estimated).^2)
end

function obj_step_weighted(params, grad, modulus, t, prescribed::Float64, measured::Vector{Float64}, weights; _insight=false)
    if _insight
        println("Current Parameters: ", params)
    end

    mod = (t->modulus(t,params))
    estimated = prescribed*mod(t)

    cost = sum((measured - estimated[weights]).^2)
end

function leastsquares_stepinit(params_init::Vector{RheoFloat},
                                low_bounds::RheovecOrNone,
                                hi_bounds::RheovecOrNone,
                                modulus,
                                time_series::Vector{RheoFloat},
                                prescribed::RheoFloat,
                                measured::Vector{RheoFloat};
                                insight::Bool = false,
                                singularity::Bool = false,
                                rel_tol_x::Union{Real,Nothing} = nothing,
                                rel_tol_f::Union{Real,Nothing} = nothing,
                                indweights=nothing,
                                optmethod::Symbol = :LN_SBPLX,
                                opttimeout::Union{Real,Nothing} = nothing,
                                optmaxeval::Union{Integer,Nothing} = nothing)

    # initialise NLOpt.Opt object with :LN_SBPLX Subplex algorithm
    opt = Opt(optmethod, length(params_init))

    # set optimiser stopping criteria

    # wall clock timeout
    if !isnothing(opttimeout)
        opttimeout = convert(RheoFloat, opttimeout)
        maxtime!(opt, opttimeout)
    end

    # evaluation cycle ceiling
    if !isnothing(optmaxeval)
        maxeval!(opt, optmaxeval)
    end

    # input parameter change tolerance
    if !isnothing(rel_tol_x)
        rel_tol = convert(RheoFloat, rel_tol_x)
        xtol_rel!(opt, rel_tol_x)
    end

    # objective function change tolerance 
    if !isnothing(rel_tol_f)
        rel_tol = convert(RheoFloat, rel_tol_f)
        ftol_rel!(opt, rel_tol_f)
    end

    # set lower bounds and upper bounds unless they take null value
    if !isnothing(low_bounds)
        low_bounds = convert(Vector{Float64}, low_bounds)
        lower_bounds!(opt, low_bounds)
    end

    if !isnothing(hi_bounds)
        hi_bounds = convert(Vector{Float64}, hi_bounds)
        upper_bounds!(opt, hi_bounds)
    end

    # set Opt object as a minimisation objective. Use a closure for additional
    # arguments sent to object objectivefunc

    params_init = convert(Vector{Float64},params_init)
    time_series = convert(Vector{Float64},time_series)
    prescribed = convert(Float64,prescribed)
    measured = convert(Vector{Float64},measured)

    if isnothing(indweights) && !singularity
        min_objective!(opt, (params, grad) -> obj_step_nonsing(params, grad, modulus,
                                                            time_series, prescribed, measured;
                                                            _insight = insight))

    elseif isnothing(indweights) && singularity
        indices = collect(Integer, 2:length(measured))
        measured_weighted = measured[indices]
        min_objective!(opt, (params, grad) -> obj_step_weighted(params, grad, modulus,
                                                        time_series, prescribed, measured_weighted, indices;
                                                        _insight = insight))

    elseif !isnothing(indweights) && !singularity
        measured_weighted = measured[indweights]
        min_objective!(opt, (params, grad) -> obj_step_weighted(params, grad, modulus,
                                                            time_series, prescribed, measured_weighted, indweights;
                                                            _insight = insight))

    elseif !isnothing(indweights) && singularity
        indsingular = indweights[indweights.>1]
        measured_weighted = measured[indsingular]
        min_objective!(opt, (params, grad) -> obj_step_weighted(params, grad, modulus,
                                                        time_series, prescribed, measured_weighted, indsingular;
                                                        _insight = insight))

    end

    # minimise objective func, minx are the parameters resulting in minimum
    (minf, minx, ret) = NLopt.optimize(opt, params_init)

    # return all
    return (convert(RheoFloat,minf), convert(Vector{RheoFloat},minx), ret)

end
