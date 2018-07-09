#!/usr/bin/env julia

#######################
#~ Utility Functions ~#
#######################

"""
    trapz(y::Array{Float64,1}, x::Array{Float64,1}[; init::Float64 = 0.0])

Array based trapezoidal integration of y with respect to x.

Limits of integration defined by the start and end points of the arrays. 'init'
keyword argument is used for setting an initial condition.
"""
function trapz(y::Array{Float64,1}, x::Array{Float64,1}; init::Float64=0.0)::Float64

    @assert length(x)==length(y) "X and Y array length must match."

    # init*2 to simplify final division by 2
    r = 2.0*init

    # trapezoidal rule
    for i in 2:length(x)
        @inbounds r += (y[i-1] + y[i])*(x[i] - x[i-1])
    end

    # return summation
    r/2.0
end


"""
    deriv(y::Array{Float64,1}[, x::Array{Float64,1}])

Given two arrays of data, x and y, calculate dy/dx using central difference
method and forward and backward difference for array boundaries.

If only one array passed in, y, then calculates derivative y with respect to
array element indices.

# Example

```jldoctest
x = collect(0.0:0.1:5.0)

y = x.^2

dy_dx = deriv(y, x) # ~= 2*x
```
"""
function deriv(y::Array{Float64,1}, x::Array{Float64,1})::Array{Float64,1}

    # assert y and x arrays are same length
    @assert length(y)==length(x) "X and Y Array lengths must match."

    # initialise zero array of length y
    ydot = zeros(Float64, length(y))
    for i in 1:length(y)
        if i==1
            # right handed difference for lower boundary
            @inbounds ydot[i] = (y[i+1] - y[i])/(x[i+1] - x[i])
        elseif i==length(y)
            # left handed difference for upper boundary
            @inbounds ydot[i] = (y[i] - y[i-1])/(x[i] - x[i-1])
        else
            # central difference with uneven spacing
            @inbounds Δx₁ = x[i] - x[i-1]
            @inbounds Δx₂ = x[i+1] - x[i]
            @inbounds ydot[i] = (y[i+1]*Δx₁^2 + (Δx₂^2 - Δx₁^2)*y[i] - y[i-1]*Δx₂^2)/(Δx₁*Δx₂*(Δx₁ + Δx₂))
        end
    end
    ydot
end

function deriv(y::Array{Float64,1})::Array{Float64,1}

    # initialise zero array of length y
    ydot = zeros(Float64, length(y))
    for i in 1:length(y)
        if i==1
            # right handed difference for lower boundary
            @inbounds ydot[i] = (y[i+1] - y[i])
        elseif i==length(y)
            # left handed difference for upper boundary
            @inbounds ydot[i] = (y[i] - y[i-1])
        else
            # central difference with even unit 1 spacing
            @inbounds ydot[i] = (y[i+1] - y[i-1])/2.0
        end
    end
    ydot
end

"""
    mittleff(α::Float64[, β::Float64], x::Array{Float64,1})

Threaded convience wrapper around MittagLeffler.mittleff for array arguments. X
array is limited to float type in accordance with the requirements of RHEOS.
"""
function mittleff(α::Float64, x_array::Array{Float64,1})::Array{Float64,1}

    # initialise array
    y = Array{Float64}(length(x_array))
    # do static scheduled threaded loop with no bounds checking
    @threads for i = 1:length(x_array)
        @inbounds y[i] = MittagLeffler.mittleff(α, x_array[i])
    end
    # return
    y
end

function mittleff(α::Float64, β::Float64, x_array::Array{Float64,1})::Array{Float64,1}

    # initialise array
    y = Array{Float64}(length(x_array))
    # do static scheduled threaded loop with no bounds checking
    @threads for i = 1:length(x_array)
        @inbounds y[i] = MittagLeffler.mittleff(α, β, x_array[i])
    end
    # return
    y
end
# NON THREADED VERSIONS
# function mittleff(α::Float64, xList::Array{Float64,1})::Array{Float64,1}
#     # call mittagleffler within an array comprehension
#     y = [MittagLeffler.mittleff(α, x) for x in xList]
# end
#
# function mittleff(α::Float64, β::Float64, xList::Array{Float64,1})::Array{Float64,1}
#     # call mittagleffler within an array comprehension
#     y = [MittagLeffler.mittleff(α, β, x) for x in xList]
# end

function quasinull(x::Array{Float64,1})

    if x == [-1.0]
        return true
    else
        return false
    end

end

#############################
#~ Preprocessing Functions ~#
#############################


"""
    getsigma(τ::Float64, samplerate::Float64)::Float64

Generate sigma/std deviation for gaussian smoothing kernel.

Acts as a low pass filter. Information of time scale τ will be half power,
faster will be increasingly cut. Called by smoothgauss function.
"""
function getsigma(τ::Float64, samplerate::Float64)::Float64

    # smoothing frequency to reduce to half power
    smoothfreq = 1.0/τ

    # reduce to half power
    sF_halfpower = smoothfreq/sqrt(2.0*log(2.0))

    # calculate std deviation and return
    σ = samplerate/(2.0*π*sF_halfpower)
end

"""
    smoothgauss(yArray::Array{Float64,1}, τ::Float64, samplerate::Float64[; pad::String="replicate"])

Smooth a signal using a Gaussian kernel.

Essentially a low pass filter with frequencies of 1/τ being cut to approximately
half power. For other pad types available see ImageFiltering documentation.
"""
function smoothgauss(yArray::Array{Float64,1}, τ::Float64, samplerate::Float64; pad::String="replicate")::Array{Float64,1}

    # get standard deviation for Gaussian kernel
    σ = getsigma(τ, samplerate)

    # smooth signal and return
    smooth = ImageFiltering.imfilter(yArray, ImageFiltering.Kernel.reflect(ImageFiltering.Kernel.gaussian((σ,))), pad)
end

"""
    var_resample(tᵢ::Array{Float64,1}, yᵢ::Array{Float64,1}, pcntdownsample::Float64, minperiod::Float64[; minsamplenum::Int64 = 25])

Convert a fixed sample rate array to a variable sample rate, with sampling points
added according to a relative change in y, 1st derivative of y and 2nd derivative
of y (WRT x).

Currently only variable downsampling supported. pcntdown sample is approximate,
works well in some cases and very poorly in others. If required, compare resampled
length vs original length after processing has finished. If data is noisy, may
benefit from sending smoothed signal to this algorithm and either using mapback
function or interpolating onto unsmoothed data.

Algorithm works as follows. Initial samples are generated evenly spread, the number
of which is determined by the minsamplenum argument. After this, the array is
repeatedly sweeped, anywhere Δy, Δdy/dx, Δd2y/dx2 is greater than a threshold, α,
a new sample is created at the midpoint of the two tested points. This is allowed to
happen a maximum of 400 times, after which α is decreased and the process starts again.
This macro process continues until the desired pcntdownsample ratio has been reached.
See source code for more implementation details.

# Arguments

- `tᵢ`: x array (usually time)
- `yᵢ`: y array, this is the data whose derivatives will determine sample densities. (Usually prescribed variable e.g. stress or strain)
- `pcntdownsample`: The ratio of new samples to old samples. Algorithm will halt ASAP once ratio greater than this number.
- `minperiod`: Minimum allowed distance between x array points, recommended to set > 0.0 to something like dx/10.0 to avoid algorithm over-focusing on a particular region.
- `minsamplenum = 25`: (Optional) number of initial, equally spaced seed samples required for algorithm to initialise.
"""
function var_resample(tᵢ::Array{Float64,1}, yᵢ::Array{Float64,1}, pcntdownsample::Float64, minperiod::Float64; minsamplenum::Int64 = 25)::Tuple{Array{Float64,1},Array{Float64,1}}

    @assert length(tᵢ)==length(yᵢ) "X and Y arrays must have same length."

    # length of original arrays
    originalSize = length(tᵢ)

    # interpolated, callable versions of y and z arrays. No need to interpolate
    # t as it is linear (fixed dt, monotonically increasing)
    yInterp = Interpolations.interpolate((tᵢ,), yᵢ, Interpolations.Gridded(Interpolations.Linear()))

    # initialise arrays for resampled data
    xInit = zeros(Float64, minsamplenum)
    yInit = zeros(Float64, minsamplenum)

    # variables for generating correct intervals during initial sweep
    minsamplenum -= 1 # correction due integer rounding
    modNum = div(originalSize, minsamplenum)

    # initial sweep
    initCounter = 1
    @inbounds for i in 1:originalSize
        if mod(i-1, modNum) == 0
            xInit[initCounter] = tᵢ[i]
            yInit[initCounter] = yᵢ[i]
            initCounter += 1
        end
    end

    # send initial resampled arrays to a dict for more convenient mutability
    xyDict = Dict{Float64, Float64}(xInit[i] => yInit[i] for i in 1:length(xInit))

    ## set up while loop variables
    maxIter = 400 # max number of iterations/sweeps for a particular alpha
    α = 0.4 # starting value of α, resampling sensitivity, decreases from here.

    # outer while loop checks to see if resampled array is large enough, if not
    # reduces alpha and starts inner while loop sweeps again. Min α required
    # to avoid strange bug
    @inbounds while (length(collect(keys(xyDict))) <= pcntdownsample*originalSize) && α>=0.002

        iₒ = 0 # counter for inner loop to ensure doesn't excedd maxIter,
        # exceeding maxIter causes bugs, need to investigate, may be overflow.

        testFailed = true # assume test failed to start inner while loop,
        # true if two succesive elements fail test for a particular alpha

        # inner while loop for each particular α value.
        while testFailed && iₒ<maxIter

            # assume test passes
            testFailed = false

            # get x and y values
            x = sort(collect(keys(xyDict)))
            y = [xyDict[i] for i in x]

            # get y, dy/dx, d2y/dx2 normalised for sweeps
            yNorm = y/maximum(abs.(y))

            dy = deriv(y, x)
            dyNorm = dy/maximum(abs.(dy))

            ddy = deriv(dy, x)
            ddyNorm = ddy/maximum(abs.(ddy))

            # sweep through array, testing every consecutive pair of elements
            for i in 1:(length(x)-1)
                logic0 = abs(yNorm[i+1] - yNorm[i])>α # 0th derivative test
                logic1 = abs(dyNorm[i+1] - dyNorm[i])>α # 1st derivative test
                logic2 = abs(ddyNorm[i+1] - ddyNorm[i])>α # 2nd derivative test

                # test failing criterion
                if  (x[i+1]-x[i]>=0.95*minperiod) && (logic0 || logic1 || logic2) && (length(collect(keys(xyDict))) <= pcntdownsample*originalSize)

                    # create new x point and y point between two points which failed test
                    xNew = (x[i+1]+x[i])/2
                    yNew = yInterp[xNew]

                    xyDict[xNew]=yNew

                    testFailed = true

                # break criterion if resampled array has enough points
                elseif (length(collect(keys(xyDict))) > pcntdownsample*originalSize)
                    break

                end # if
            end # for
            iₒ += 1 # increment inner while loop counter
        end # inner while

        α -= 0.001 # decrement α for outer while loop

    end #outer while

    finalX = sort(collect(keys(xyDict)))
    finalY = [xyDict[i] for i in finalX]

    return finalX, finalY
end

"""
    closestindex(x::Array{Float64,1}, val::Float64)::Int32

Find the index of the array element closest to val.
"""
function closestindex(x::Array{Float64,1}, val::Float64)::Int32

    # intialise closest match variable, assuming best match is index 1
    ibest = start(eachindex(x))

    # diff between value and current element
    dxbest = abs(x[ibest]-val)

    # loop through all elements, looking for smallest difference
    for I in eachindex(x)
        @inbounds dx = abs(x[I]-val)
        if dx < dxbest
            dxbest = dx
            ibest = I
        end
    end

    ibest
end

"""
    closestindices(x::Array{Float64,1}, vals::Array{Float64,1})

Uses `closestindex` iteratively to find closest index for all values in `vals` array.
"""
function closestindices(x::Array{Float64,1}, vals::Array{Float64,1})::Array{Int64,1}

    indicesbest = Array{Int64,1}(length(vals))
    # call closest index iteratively
    for (i, v) in enumerate(vals)
        indicesbest[i] = closestindex(x, v)
    end

    indicesbest

end

"""
    mapback(xᵦ::Array{Float64,1}, x::Array{Float64,1})

Match a variable downsampled array xᵦ to its closest possible elements in original
array x. Returns the array of unique indices corresponding to these matched elements.

Can be used after variable resampling to maintain resampling priorities whilst
not interfering with data of the fidelity by use of interpolations. Can also be
useful to ensure no region of the original data has been oversampled as duplicate
indices are deleted.
"""
function mapback(xᵦ::Array{Float64,1}, x::Array{Float64,1})

    # array to store indices
    indices = zeros(Int32, length(xᵦ))

    # loop through all values in resampled array, finding closest index in original array
    for i in 1:length(xᵦ)
        @inbounds iClosest = closestindex(x, xᵦ[i])
        @inbounds indices[i] = iClosest
    end

    # remove duplicates
    indices = unique(indices)

    # return indices, can be then used like xᵩ = x[indices]
    indices
end

"""
    downsample(boundaries::Array{Int64,1}, elperiods::Array{Int64,1})

Reduce indices.

Can be applied to the whole array using boundaries
[1,length(array_to_downsample)] or using different spacers for different parts
of the data.

# Example

```julia-repl
julia> x = collect(1.0:1.0:15.0);

julia> newIndices = downsample([1,5,11,length(x)], [2,3,1]);

julia> println(x)
[1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0]

julia> println(x[newIndices])
[1.0, 3.0, 5.0, 8.0, 11.0, 12.0, 13.0, 14.0, 15.0]
```
"""
function downsample(boundaries::Array{Int64,1}, elperiods::Array{Int64,1})

    # assert correct function signature
    @assert length(elperiods)==length(boundaries)-1 "Number of different sample periods must be 1 less than boundaries provided"

    # initialise indices array
    indices = zeros(Int64,0)

    # loop through, skipping elements as required
    for i in 1:length(boundaries)-1

        # get indices for this 'section'
        indicesCur = boundaries[i]:elperiods[i]:boundaries[i+1]

        # append
        append!(indices, collect(indicesCur))
    end

    # return unique indices
    indices = unique(indices)
end

"""
    fixed_resample(x::Array{Float64,1}, y::Array{Float64,1}, boundaries::Array{Int64,1}, elperiods::Array{Int64,1}, direction::Array{String,1})

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
function fixed_resample(x::Array{Float64,1}, y::Array{Float64,1},
                        boundaries::Array{Int64,1}, elperiods::Array{Int64,1},
                        direction::Array{String,1})

    # assert correct function signature
    @assert length(x)==length(y) "X and Y arrays must have same length."
    @assert length(elperiods)==length(boundaries)-1 "Number of different sample periods must be 1 less than boundaries provided"

    # y as callable interpolations, used for upsampled regions
    yInterp = Interpolations.interpolate((x,), y, Interpolations.Gridded(Interpolations.Linear()))

    # initialise resampled arrays as empty
    xᵦ = zeros(Float64,0)
    yᵦ = zeros(Float64,0)

    # loop through boundaries
    for i in 1:length(boundaries)-1

        # upsampling, starts at each element then intepolates up N times
        if direction[i]=="up"
            for k in boundaries[i]:(boundaries[i+1]-1)
                # starting element
                append!(xᵦ,x[k])
                append!(yᵦ,yInterp[x[k]])

                # dx increment up
                Δ = (x[k+1]-x[k])/elperiods[i]
                for N in 1:(elperiods[i]-1)
                    # add N increments to starting element
                    append!(xᵦ,x[k]+N*Δ)
                    append!(yᵦ,yInterp[x[k]+N*Δ])
                end
            end

        # downsampling, simply takes every N element as in downsample function
        elseif direction[i]=="down" #under/sampling
            append!(xᵦ,x[boundaries[i]:elperiods[i]:(boundaries[i+1]-1)])
            append!(yᵦ,y[boundaries[i]:elperiods[i]:(boundaries[i+1]-1)])
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

"""
    boltzintegral_nonsing(modulus::Function, time_series::Array{Float64,1}, params::Array{Float64,1}, prescribed_dot::Array{Float64,1})

Calculate Boltzmann Superposition integral using direct integration method.

This is much slower and slightly less accurate (depending on sample resolution)
than the convolution method. However, it works for variable sample rate.

# Arguments

- `modulus`: Viscoelastic modulus function
- `time_series`: The array of times
- `params`: Parameters passed to viscoelastic modulus
- `prescribed_dot`: Derivative of (usually prescribed) variable inside the integration kernel
"""
function boltzintegral_nonsing(modulus::Function, time_series::Array{Float64,1}, params::Array{Float64,1},
                    prescribed_dot::Array{Float64,1})::Array{Float64,1}

    I = zeros(length(time_series))
    for (i,v) in enumerate(time_series[1:end])
        # generate integral for each time step
        @inbounds τ = time_series[1:i]
        Modulus_arg = v - τ
        Modulusᵢ = modulus(Modulus_arg, params)
        @inbounds df_dtᵢ = prescribed_dot[1:i]
        intergrand = Modulusᵢ.*df_dtᵢ
        @inbounds I[i] = trapz(intergrand, τ)
    end

    I

end

"""
    boltzintegral_sing(modulus::Function, time_series::Array{Float64,1}, params::Array{Float64,1}, prescribed_dot::Array{Float64,1})

Calculate Boltzmann Superposition integral using direct integration method.

This is much slower and slightly less accurate (depending on sample resolution)
than the convolution method. However, it works for variable sample rate.

Should be used when viscoelastic model contains a singularity and should be compared
with [2:end] of reference array when fitting.

# Arguments

- `modulus`: Viscoelastic modulus function
- `time_series`: The array of times
- `params`: Parameters passed to viscoelastic modulus
- `prescribed_dot`: Derivative of (usually prescribed) variable inside the integration kernel
"""
function boltzintegral_sing(modulus::Function, time_series::Array{Float64,1}, params::Array{Float64,1},
                    prescribed_dot::Array{Float64,1})::Array{Float64,1}

    I = zeros(length(time_series)-1)
    # only traverse after time t=0, first element
    for (i,v) in enumerate(time_series[2:end])
        # generate integral for each time step
        @inbounds τ = time_series[1:i]
        Modulus_arg = v - τ
        Modulusᵢ = modulus(Modulus_arg, params)
        @inbounds df_dtᵢ = prescribed_dot[1:i]
        intergrand = Modulusᵢ.*df_dtᵢ
        @inbounds I[i] = trapz(intergrand, τ)
    end

    I

end

"""
    boltzconvolve_nonsing(modulus::Function, time_series::Array{Float64,1}, dt_series::Array{Float64,1}, params::Array{Float64,1}, prescribed_dot::Array{Float64,1})

Calculate Boltzmann Superposition integral using convolution method.

This is much faster and slightly more accurate (depending on sample resolution)
than the integral method. However, it works for constant sample rate.

# Arguments

- `modulus`: Viscoelastic modulus function
- `time_series`: The array of times
- `dt_series`: Array of gradient of time series, can be found using deriv(time_series)
- `params`: Parameters passed to viscoelastic modulus
- `prescribed_dot`: Derivative of (usually prescribed) variable inside the integration kernel
"""
function boltzconvolve_nonsing(modulus::Function, time_series::Array{Float64,1}, dt_series::Array{Float64,1},
                        params::Array{Float64,1}, prescribed_dot::Array{Float64,1})::Array{Float64,1}

    Modulus = modulus(time_series, params)
    β = FastConv.convn(Modulus, prescribed_dot)
    # pick out relevant elements (1st half) and multiply by dt
    β = β[1:length(dt_series)].*dt_series

end

"""
    boltzconvolve_sing(modulus::Function, time_series::Array{Float64,1}, dt_series::Array{Float64,1}, params::Array{Float64,1}, prescribed_dot::Array{Float64,1})

Calculate Boltzmann Superposition integral using convolution method.

This is much faster and slightly more accurate (depending on sample resolution)
than the integral method. However, it works for constant sample rate.

Should be used when singularity exists in viscoelastic model and should be compared
with [2:end] of reference array when fitting.

# Arguments

- `modulus`: Viscoelastic modulus function
- `time_series`: The array of times
- `dt_series`: Array of gradient of time series, can be found using deriv(time_series)
- `params`: Parameters passed to viscoelastic modulus
- `prescribed_dot`: Derivative of (usually prescribed) variable inside the integration kernel
"""
function boltzconvolve_sing(modulus::Function, time_series::Array{Float64,1}, dt_series::Array{Float64,1},
                        params::Array{Float64,1}, prescribed_dot::Array{Float64,1})::Array{Float64,1}

    # convolved length will be original_length-1
    len = length(dt_series)-1
    Modulus = modulus(time_series, params)
    # fast convolution, ignoring initial singularity
    β = FastConv.convn(Modulus[2:end], prescribed_dot[1:end])
    # pick out relevant elements (1st half) and multiply by dt
    β = β[1:len].*dt_series[1:len]

end

"""
    obj_const_nonsing(params::Array{Float64,1}, grad::Array{Float64,1}, modulus::Function, time_series::Array{Float64,1}, dt_series::Array{Float64,1}, prescribed_dot::Array{Float64,1}, measured::Array{Float64,1}; _insight::Bool = false)

Generate the sum-of-squares of the difference between `measured` and the Boltzmann
convolution integral of the selected `modulus` and `prescribed_dot`. Used when
sample rate is constant and model does not feature singularity.

# Arguments

- `params`: Array of parameters sent to `modulus`
- `grad`: Gradient argument used by NLOpt
- `modulus`: Viscoelastic modulus function
- `time_series`: Array of time data
- `dt_series`: Array of time data differences,  equal to deriv(time_series)
- `prescribed_dot`: Convolved with `modulus` over `time_series`, usually dϵ/dt for stress relaxation and dσ/dt for creep
- `measured`: Data for comparison against, usually σ for stress relaxation and ϵ for creep
- `_insight`: Declare whether insight info should be shown when this function is called, true or false
"""
function obj_const_nonsing(params::Array{Float64,1}, grad::Array{Float64,1},
                            modulus::Function, time_series::Array{Float64,1},
                            dt_series::Array{Float64,1}, prescribed_dot::Array{Float64,1},
                            measured::Array{Float64,1}; _insight::Bool = false)::Float64

    if _insight
        println("Current Parameters: ", params)
    end

    convolved = boltzconvolve_nonsing(modulus, time_series, dt_series, params, prescribed_dot)

    cost = sum(0.5*(measured - convolved).^2)

end

"""
    obj_const_sing(params::Array{Float64,1}, grad::Array{Float64,1}, modulus::Function, time_series::Array{Float64,1}, dt_series::Array{Float64,1}, prescribed_dot::Array{Float64,1}, measured::Array{Float64,1}; _insight::Bool = false)

Generate the sum-of-squares of the difference between `measured` and the Boltzmann
convolution integral of the selected `modulus` and `prescribed_dot`. Used when
sample rate is constant and model does feature singularity.

# Arguments

- `params`: Array of parameters sent to `modulus`
- `grad`: Gradient argument used by NLOpt
- `modulus`: Viscoelastic modulus function
- `time_series`: Array of time data
- `dt_series`: Array of time data differences,  equal to deriv(time_series)
- `prescribed_dot`: Convolved with `modulus` over `time_series`, usually dϵ/dt for stress relaxation and dσ/dt for creep
- `measured`: Data for comparison against, usually σ for stress relaxation and ϵ for creep
- `_insight`: Declare whether insight info should be shown when this function is called, true or false
"""
function obj_const_sing(params::Array{Float64,1}, grad::Array{Float64,1},
                            modulus::Function, time_series::Array{Float64,1},
                            dt_series::Array{Float64,1}, prescribed_dot::Array{Float64,1},
                            measured::Array{Float64,1}; _insight::Bool = false)::Float64

    if _insight
        println("Current Parameters: ", params)
    end

    convolved = boltzconvolve_sing(modulus, time_series, dt_series, params, prescribed_dot)

    # don't use first element as singularity exists in model
    cost = sum(0.5*(measured[2:end] - convolved).^2)

end

"""
    obj_var_nonsing(params::Array{Float64,1}, grad::Array{Float64,1}, modulus::Function, time_series::Array{Float64,1}, prescribed_dot::Array{Float64,1}, measured::Array{Float64,1}; _insight::Bool = false)

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
function obj_var_nonsing(params::Array{Float64,1}, grad::Array{Float64,1},
                            modulus::Function, time_series::Array{Float64,1},
                            prescribed_dot::Array{Float64,1}, measured::Array{Float64,1};
                            _insight::Bool = false)::Float64

    if _insight
        println("Current Parameters: ", params)
    end

    convolved = boltzintegral_nonsing(modulus, time_series, params, prescribed_dot)

    cost = sum(0.5*(measured - convolved).^2)

end

"""
    obj_var_sing(params::Array{Float64,1}, grad::Array{Float64,1}, modulus::Function, time_series::Array{Float64,1}, prescribed_dot::Array{Float64,1}, measured::Array{Float64,1}; _insight::Bool = false)

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
function obj_var_sing(params::Array{Float64,1}, grad::Array{Float64,1},
                            modulus::Function, time_series::Array{Float64,1},
                            prescribed_dot::Array{Float64,1}, measured::Array{Float64,1};
                            _insight::Bool = false)::Float64

    if _insight
        println("Current Parameters: ", params)
    end

    convolved = boltzintegral_sing(modulus, time_series, params, prescribed_dot)

    # don't use first element as singularity exists in model
    cost = sum(0.5*(measured[2:end] - convolved).^2)

end

"""
    leastsquares_init(params_init, low_bounds, hi_bounds, modulus, time_series, dt_series, prescribed_dot, measured; insight = false, sampling = "constant", singularity = false)

Initialise then begin a least squares fitting of the supplied data.

# Arguments

- `params_init`: Initial parameters to be used (starting guess)
- `low_bounds`: Lower bounds for parameters
- `hi_bounds`: Higher bounds for parameters
- `modulus`: Viscoelastic modulus function
- `time_series`: Array of time data
- `dt_series`: Array of time data differences,  equal to deriv(time_series)
- `prescribed_dot`: Convolved with `modulus` over `time_series`, usually dϵ/dt for stress relaxation and dσ/dt for creep
- `measured`: Data for comparison against, usually σ for stress relaxation and ϵ for creep
- `sampling`: Declare whether sample rate is `constant` or `variable` so that convolution or integration is used respectively
- `insight`: Declare whether insight info should be shown when this function is called, true or false
- `singularity`: Presence of singularity in model
"""
function leastsquares_init(params_init::Array{Float64,1}, low_bounds::Array{Float64,1},
                           hi_bounds::Array{Float64,1}, modulus::Function,
                           time_series::Array{Float64,1}, dt_series::Array{Float64,1},
                           prescribed_dot::Array{Float64,1}, measured::Array{Float64,1};
                           insight::Bool = false, sampling::String = "constant",
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
    if !singularity && sampling == "constant"
        min_objective!(opt, (params, grad) -> obj_const_nonsing(params, grad, modulus,
                                                            time_series, dt_series,
                                                            prescribed_dot, measured;
                                                            _insight = insight))

    elseif singularity && sampling == "constant"
        min_objective!(opt, (params, grad) -> obj_const_sing(params, grad, modulus,
                                                        time_series, dt_series,
                                                        prescribed_dot, measured;
                                                        _insight = insight))

    elseif !singularity && sampling == "variable"
        min_objective!(opt, (params, grad) -> obj_var_nonsing(params, grad, modulus,
                                                        time_series, prescribed_dot,
                                                        measured; _insight = insight))

    elseif singularity && sampling == "variable"
        min_objective!(opt, (params, grad) -> obj_var_sing(params, grad, modulus,
                                                        time_series, prescribed_dot,
                                                        measured; _insight = insight))

    end

    # minimise objective func, minx are the parameters resulting in minimum
    (minf, minx, ret) = NLopt.optimize(opt, params_init)

    # return all
    return (minf, minx, ret)

end
