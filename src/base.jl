#!/usr/bin/env julia

#######################
#~ Utility Functions ~#
#######################

"""
    trapz(y, x; init=0.0)

Array based trapezoidal integration of y with respect to x.

Limits of integration defined by the start and end points of the arrays. 'init'
keyword argument is used for setting an initial condition.
"""
function trapz(y, x; init=0.0)

    n = length(x)

    @assert n==length(y) "X and Y array length must match."

    # init*2 to simplify final division by 2
    r = 2.0*init

    if n==1; return r; end

    # trapezoidal rule
    @inbounds for i in 2:length(x)
        r += (y[i-1] + y[i])*(x[i] - x[i-1])
    end

    # return summation
    r/2.0

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
    @inbounds ydot[1] = (y[1] - 0.0)/(x[2] - x[1])

    # central difference with uneven spacing for general case of constant or variable sample rate
    @inbounds for i in 2:(N-1)
        Δx₁ = x[i] - x[i-1]
        Δx₂ = x[i+1] - x[i]
        ydot[i] = (y[i+1]*Δx₁^2 + (Δx₂^2 - Δx₁^2)*y[i] - y[i-1]*Δx₂^2)/(Δx₁*Δx₂*(Δx₁ + Δx₂))
    end

    # 1st order backwards difference for last element
    ydot[N] = (y[N] - y[N-1])/(x[N] - x[N-1])

    return convert(Vector{RheoFloat},ydot)

end

"""
    derivBD(y, x)

Given two arrays of data, x and y, calculate dy/dx using 1st order
backward difference. Assumes y==0 at a previous point, i.e.
y is 'at rest'. Captures instantaneous loading where derivCD will not.
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
    @inbounds ydot[1] = (y[1] - 0.0)/(x[2] - x[1])

    # backwards difference method for rest of points
    @inbounds for i in 2:N
        ydot[i] = (y[i] - y[i-1])/(x[i] - x[i-1])
    end

    return convert(Vector{RheoFloat},ydot)
end

function quasinull(x::Array{RheoFloat,1})

    if x == convert(Array{RheoFloat,1},[-1.0])
        return true
    else
        return false
    end

end

function constantcheck(t::Vector{RheoFloat} )

    diff = round.(t[2:end]-t[1:end-1]; digits=4)
    # check if any element is not equal to 1st element
    check = !any(x -> x != diff[1], diff)

end

function getsampleperiod(t::Vector{RheoFloat})

    @assert constantcheck(t) "Sample-rate must be constant"

    rate = t[2]-t[1]

end

#############################
#~ Preprocessing Functions ~#
#############################

"""
    getsigma(τ, samplerate)

Generate sigma/std deviation for gaussian smoothing kernel.

Acts as a low pass filter. Information of time scale τ will be half power,
faster will be increasingly cut. Called by smoothgauss function.
"""
function getsigma(τ::Real, samplerate::Real)

    # get freq, reduce to half power and generate gaussian std. deviation
    smoothfreq = 1.0/τ

    sF_halfpower = smoothfreq/sqrt(2.0*log(2.0))

    σ = samplerate/(2.0*π*sF_halfpower)

end

"""
    closestindex(x, val)

Find the index of the array element closest to val.
"""
function closestindex(x::Vector{T1}, val::T2) where {T1<:Real,T2<:Real}

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

"""
    fixed_resample(x::Array{RheoFloat,1}, y::Array{RheoFloat,1}, boundaries::Array{Int64,1}, elperiods::Array{Int64,1}, direction::Array{String,1})

Resample two arrays with new sample rate(s).

Whereas downsample can only reduce the sample rate by not taking every array element,
fixed_resample can also upsample. This is why it does not work on and return indices
but rather the new resampled arrays themselves.

# Example

```julia-repl
julia> x = collect(1.0:1.0:7.0);

julia> (x1,y1,z1) = fixed_resample(x, x, x, [1,3,5,length(x)], [2,1,4], ["up","down","up"]);

julia> println(x)
[1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0]

julia> println(x1)
[1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 5.25, 5.5, 5.75, 6.0, 6.25, 6.5, 6.75, 7.0]
```
"""
function fixed_resample(x::Vector{T}, y::Vector{T},
                        boundaries::Vector{U}, elperiods::Union{Vector{U},U}) where T<:RheoFloat where U<:Integer

    #@eval using Interpolations

    # assert correct function signature
    @assert length(x)==length(y) "X and Y arrays must have same length."
    @assert length(elperiods)==length(boundaries)-1 "Number of different sample periods must be 1 less than boundaries provided"

    # y as callable interpolations, used for upsampled regions
    yInterp = Base.invokelatest(interpolate, (x,), y, Base.invokelatest(Gridded, Base.invokelatest(Linear)))

    # initialise resampled arrays as empty
    xᵦ = zeros(RheoFloat,0)
    yᵦ = zeros(RheoFloat,0)

    # loop through boundaries
    for i in 1:length(boundaries)-1

        # upsampling, starts at each element then intepolates up N times
        if !signbit(elperiods[i])
            for k in boundaries[i]:(boundaries[i+1]-1)
                # starting element
                append!(xᵦ, x[k])
                append!(yᵦ, yInterp[x[k]])

                # dx increment up
                Δ = (x[k+1]-x[k])/elperiods[i]
                for N in 1:(elperiods[i]-1)
                    # add N increments to starting element
                    append!(xᵦ, x[k]+N*Δ)
                    append!(yᵦ, yInterp[x[k]+N*Δ])
                end
            end

        # downsampling, simply takes every N element as in downsample function
        elseif signbit(elperiods[i]) #under/sampling
            append!(xᵦ,x[boundaries[i]:abs(elperiods[i]):(boundaries[i+1]-1)])
            append!(yᵦ,y[boundaries[i]:abs(elperiods[i]):(boundaries[i+1]-1)])
        end
    end

    # last element is missed out due to algorithm spec, needs to be added in after
    lastel = boundaries[end]
    # safety check then append
    if xᵦ[end]<x[lastel]
        append!(xᵦ,x[lastel])
        append!(yᵦ,y[lastel])
    end

    return xᵦ, yᵦ
end

##########################
#~ Processing Functions ~#
##########################

function singularitytest(modulus::Function, params::Array{RheoFloat, 1}; t1::RheoFloat=convert(RheoFloat,0.0))


    startval = modulus([t1], params)[1]

    if isnan(startval) || startval == Inf
        return true
    else
        return false
    end

end

"""
    boltzintegral_nonsing(modulus::Function, time_series::Array{RheoFloat,1}, params::Array{RheoFloat,1}, prescribed_dot::Array{RheoFloat,1})

Calculate Boltzmann Superposition integral using direct integration method.

This is much slower and slightly less accurate (depending on sample resolution)
than the convolution method. However, it works for variable sample rate.

# Arguments

- `modulus`: Viscoelastic modulus function
- `time_series`: The array of times
- `params`: Parameters passed to viscoelastic modulus
- `prescribed_dot`: Derivative of (usually prescribed) variable inside the integration kernel
"""
function boltzintegral_nonsing(modulus, time_series, params,prescribed_dot)

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
        Modulusᵢ = modulus(Modulus_arg, params)
        df_dtᵢ = prescribed_dot_mod[1:i]
        intergrand = Modulusᵢ.*df_dtᵢ

        I[i] = trapz(intergrand, τ)
    end
    # fix initial point
    # I[2] = (prescribed_dot[1]*modulus([0.0], params)*(time_series[2] - time_series[1]))[1]
    # to catch weird bug in InverseLaplace
    I[2] = (prescribed_dot[1]*modulus(time_series, params)*(time_series[2] - time_series[1]))[1]

    I[2:end]


end

# """
#     boltzintegral_sing(modulus::Function, time_series::Array{RheoFloat,1}, params::Array{RheoFloat,1}, prescribed_dot::Array{RheoFloat,1})

# Calculate Boltzmann Superposition integral using direct integration method.

# This is much slower and slightly less accurate (depending on sample resolution)
# than the convolution method. However, it works for variable sample rate.

# Should be used when viscoelastic model contains a singularity and should be compared
# with [2:end] of reference array when fitting.

# # Arguments

# - `modulus`: Viscoelastic modulus function
# - `time_series`: The array of times
# - `params`: Parameters passed to viscoelastic modulus
# - `prescribed_dot`: Derivative of (usually prescribed) variable inside the integration kernel
# """
# function boltzintegral_sing(modulus::Function, time_series::Array{RheoFloat,1}, params::Array{RheoFloat,1},
#                     prescribed_dot::Array{RheoFloat,1})::Array{RheoFloat,1}

#     I = zeros(length(time_series)-1)
#     # only traverse after time t=0, first element
#     @inbounds for (i, v) in enumerate(time_series[2:end])
#         # generate integral for each time step
#         τ = time_series[1:i]
#         Modulus_arg = v - τ
#         Modulusᵢ = modulus(Modulus_arg, params)
#         df_dtᵢ = prescribed_dot[1:i]
#         intergrand = Modulusᵢ.*df_dtᵢ
#         I[i] = trapz(intergrand, τ)
#     end

#     I

# end

function boltzintegral_sing(modulus, time_series, params,prescribed_dot)


    # init time diff, used to cope with singularity
    init_offset = (time_series[2] - time_series[1])/10.0

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
        Modulus_arg[end] = init_offset
        Modulusᵢ = modulus(Modulus_arg, params)
        df_dtᵢ = prescribed_dot_mod[1:i]

        intergrand = Modulusᵢ.*df_dtᵢ

        I[i] = trapz(intergrand, τ)

    end

    # fix initial point
    I[2] = (prescribed_dot[1]*modulus([init_offset], params)*(time_series[2] - time_series[1]))[1]

    I[2:end]

end



# """
#     boltzconvolve_sing(modulus::Function, time_series::Array{RheoFloat,1}, dt::RheoFloat, params::Array{RheoFloat,1}, prescribed_dot::Array{RheoFloat,1})

# Calculate Boltzmann Superposition integral using convolution method.

# This is much faster and slightly more accurate (depending on sample resolution)
# than the integral method. However, it works for constant sample rate.

# Should be used when singularity exists in viscoelastic model and should be compared
# with [2:end] of reference array when fitting.

# # Arguments

# - `modulus`: Viscoelastic modulus function
# - `time_series`: The array of times
# - `dt`: Constant time step (sample period)
# - `params`: Parameters passed to viscoelastic modulus
# - `prescribed_dot`: Derivative of (usually prescribed) variable inside the integration kernel
# """
# function boltzconvolve_sing(modulus::Function, time_series::Array{RheoFloat,1}, dt::RheoFloat,
#                         params::Array{RheoFloat,1}, prescribed_dot::Array{RheoFloat,1})::Array{RheoFloat,1}

#     # convolved length will be original_length-1
#     len = length(time_series)-1
#     Modulus = modulus(time_series, params)
#     # fast convolution, ignoring initial singularity
#     β = convn(Modulus[2:end], prescribed_dot[1:end])
#     # pick out relevant elements (1st half) and multiply by dt
#     β = β[1:len]*dt

# end

function boltzconvolve(modulus, time_series, dt,params, prescribed_dot)


    Modulus = modulus(time_series, params)
    # fast convolution
    β = convn(Modulus, prescribed_dot)
    # pick out relevant elements (1st half) and multiply by dt
    β = β[1:length(time_series)]*dt
    return convert(Vector{RheoFloat},β)

end

"""
    boltzconvolve_nonsing(modulus::Function, time_series::Array{RheoFloat,1}, dt::RheoFloat, params::Array{RheoFloat,1}, prescribed_dot::Array{RheoFloat,1})

Calculate Boltzmann Superposition integral using convolution method.

This is much faster and slightly more accurate (depending on sample resolution)
than the integral method. However, it works for constant sample rate.

# Arguments

- `modulus`: Viscoelastic modulus function
- `time_series`: The array of times
- `dt`: Constant time step (sample period)
- `params`: Parameters passed to viscoelastic modulus
- `prescribed_dot`: Derivative of (usually prescribed) variable inside the integration kernel
"""
function boltzconvolve_nonsing(modulus, time_series, dt,params, prescribed_dot)

    Modulus = modulus(time_series, params)

    β = convn(Modulus, prescribed_dot)
    # pick out relevant elements (1st half) and multiply by dt
    β = β[1:length(time_series)]*dt

end

"""
    obj_const_nonsing(params::Array{RheoFloat,1}, grad::Array{RheoFloat,1}, modulus::Function, time_series::Array{RheoFloat,1}, dt::RheoFloat, prescribed_dot::Array{RheoFloat,1}, measured::Array{RheoFloat,1}; _insight::Bool = false)

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
function obj_const_nonsing(params, grad,modulus, time_series,dt, prescribed_dot,measured; _insight::Bool = false)

    if _insight
        println("Current Parameters: ", params)
    end
    convolved = boltzconvolve_nonsing(modulus, time_series, dt, params, prescribed_dot)

    cost = sum(0.5*(measured - convolved).^2)
    return convert(RheoFloat,cost)

end

"""
    obj_const_sing(params::Array{RheoFloat,1}, grad::Array{RheoFloat,1}, modulus::Function, time_series::Array{RheoFloat,1}, dt::RheoFloat, prescribed_dot::Array{RheoFloat,1}, measured::Array{RheoFloat,1}; _insight::Bool = false)

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
function obj_const_sing(params, grad,modulus, time_series,dt, prescribed_dot,measured; _insight::Bool = false)

    if _insight
        println("Current Parameters: ", params)
    end

    # convolved = boltzconvolve_sing(modulus, time_series, dt, params, prescribed_dot)
    convolved = boltzconvolve(modulus, time_series, dt, params, prescribed_dot)

    # do as this has been taken care of in convolution!
    cost = sum(0.5*(measured - convolved).^2)

end

"""
    obj_var_nonsing(params::Array{RheoFloat,1}, grad::Array{RheoFloat,1}, modulus::Function, time_series::Array{RheoFloat,1}, prescribed_dot::Array{RheoFloat,1}, measured::Array{RheoFloat,1}; _insight::Bool = false)

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
function obj_var_nonsing(params, grad,modulus, time_series, prescribed_dot, measured;_insight::Bool = false)

    if _insight
        println("Current Parameters: ", params)
    end
    #print("entrato1")
    convolved = boltzintegral_nonsing(modulus, time_series, params, prescribed_dot)
    #print("entrato2")
    #print(typeof(convolved))

    cost = sum(0.5*(measured - convolved).^2)

end

"""
    obj_var_sing(params::Array{RheoFloat,1}, grad::Array{RheoFloat,1}, modulus::Function, time_series::Array{RheoFloat,1}, prescribed_dot::Array{RheoFloat,1}, measured::Array{RheoFloat,1}; _insight::Bool = false)

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
function obj_var_sing(params, grad,modulus, time_series,prescribed_dot, measured;_insight::Bool = false)

    if _insight
        println("Current Parameters: ", params)
    end

    convolved = boltzintegral_sing(modulus, time_series, params, prescribed_dot)

    # don't use first element as singularity exists in model - Edit: INCORRECT COMMENT!
    # cost = sum(0.5*(measured[2:end] - convolved).^2)

    # do as this has been taken care of in convolution!
    cost = sum(0.5*(measured - convolved).^2)

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
- `sampling`: Declare whether sample rate is `constant` or `variable` so that convolution or integration is used respectively
- `insight`: Declare whether insight info should be shown when this function is called, true or false
- `singularity`: Presence of singularity in model
"""
function leastsquares_init(params_init::Vector{RheoFloat}, low_bounds::Vector{RheoFloat},
                           hi_bounds::Vector{RheoFloat}, modulus::Function,
                           time_series::Vector{RheoFloat}, dt::RheoFloat,
                           prescribed_dot::Vector{RheoFloat}, measured::Vector{RheoFloat};
                           insight::Bool = false, sampling::Bool=true,
                           singularity::Bool = false, _rel_tol = 1e-4)

    # initialise NLOpt.Opt object with :LN_SBPLX Subplex algorithm
    opt = Opt(:LN_SBPLX, length(params_init))

    # set lower bounds and upper bounds unless they take null value of [-1.0]
    if !quasinull(low_bounds)
        lower_bounds!(opt, low_bounds)
    end

    if !quasinull(hi_bounds)
        upper_bounds!(opt, hi_bounds)
    end

    # set relative tolerance
    xtol_rel!(opt, _rel_tol)

    # set Opt object as a minimisation objective. Use a closure for additional
    # arguments sent to object objectivefunc

    params_init = convert(Vector{Float64},params_init)
    low_bounds = convert(Vector{Float64},low_bounds)
    hi_bounds = convert(Vector{Float64}, hi_bounds)
    time_series = convert(Vector{Float64},time_series)
    prescribed_dot = convert(Vector{Float64},prescribed_dot)
    measured = convert(Vector{Float64},measured)


    if !singularity && sampling
        min_objective!(opt, (params, grad) -> obj_const_nonsing(params, grad, modulus,
                                                            time_series, dt,
                                                            prescribed_dot, measured;
                                                            _insight = insight))

    elseif singularity && sampling
        # remove singularity, just go close to it, 1/10th over first sample period
        time_series[1] = 0.0 + (time_series[2] - time_series[1])/10.0

        min_objective!(opt, (params, grad) -> obj_const_sing(params, grad, modulus,
                                                        time_series, dt,
                                                        prescribed_dot, measured;
                                                        _insight = insight))

    elseif !singularity && !sampling
        min_objective!(opt, (params, grad) -> obj_var_nonsing(params, grad, modulus,
                                                        time_series, prescribed_dot,
                                                        measured; _insight = insight))

    elseif singularity && !sampling
        min_objective!(opt, (params, grad) -> obj_var_sing(params, grad, modulus,
                                                        time_series, prescribed_dot,
                                                        measured; _insight = insight))

    end

    # minimise objective func, minx are the parameters resulting in minimum
    (minf, minx, ret) = NLopt.optimize(opt, params_init)

    # return all
    return (convert(RheoFloat,minf), convert(Array{RheoFloat,1},minx), ret)


end

function obj_step_nonsing(params, grad, modulus, t, prescribed, measured; _insight=false)
    if _insight
        println("Current Parameters: ", params)
    end

    estimated = prescribed*modulus(t, params)

    cost = sum(0.5*(measured - estimated).^2)
end

function obj_step_sing(params, grad, modulus, t, prescribed, measured; _insight=false)
    if _insight
        println("Current Parameters: ", params)
    end

    estimated = prescribed*modulus(t, params)[2:end]

    cost = sum(0.5*(measured[2:end] - estimated).^2)
end

function leastsquares_stepinit(params_init::Array{RheoFloat,1}, low_bounds::Array{RheoFloat,1},
                           hi_bounds::Array{RheoFloat,1}, modulus::Function,
                           time_series::Array{RheoFloat,1}, prescribed::RheoFloat,
                           measured::Array{RheoFloat,1}; insight::Bool = false,
                           singularity::Bool = false, _rel_tol = 1e-4)

    # initialise NLOpt.Opt object with :LN_SBPLX Subplex algorithm
    opt = Opt(:LN_SBPLX, length(params_init))

    # set lower bounds and upper bounds unless they take null value of [-1.0]
    if !quasinull(low_bounds)
        lower_bounds!(opt, low_bounds)
    end

    if !quasinull(hi_bounds)
        upper_bounds!(opt, hi_bounds)
    end

    # set relative tolerance
    xtol_rel!(opt, _rel_tol)

    # set Opt object as a minimisation objective. Use a closure for additional
    # arguments sent to object objectivefunc

    params_init = convert(Vector{Float64},params_init)
    low_bounds = convert(Vector{Float64},low_bounds)
    hi_bounds = convert(Vector{Float64}, hi_bounds)
    time_series = convert(Vector{Float64},time_series)
    prescribed = convert(Float64,prescribed)
    measured = convert(Vector{Float64},measured)

    if !singularity
        min_objective!(opt, (params, grad) -> obj_step_nonsing(params, grad, modulus,
                                                            time_series, prescribed, measured;
                                                            _insight = insight))

    elseif singularity
        min_objective!(opt, (params, grad) -> obj_step_sing(params, grad, modulus,
                                                        time_series, prescribed, measured;
                                                        _insight = insight))

    end

    # minimise objective func, minx are the parameters resulting in minimum
    (minf, minx, ret) = NLopt.optimize(opt, params_init)

    # return all
    return (convert(RheoFloat,minf), convert(Array{RheoFloat,1},minx), ret)

end
