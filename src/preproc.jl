#!/usr/bin/env julia

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
    var_resample(tᵢ::Array{Float64,1}, yᵢ::Array{Float64,1}[, zᵢ::Array{Float64,1}], pcntdownsample::Float64, minperiod::Float64[; minsamplenum::Int64 = 25])

Convert a fixed sample rate array to a variable sample rate, with sampling points
added according to a relative change in y, 1st derivative of y and 2nd derivative
of y (WRT x).

Currently only variable downsampling supported. pcntdown sample is approximate,
works well in some cases and very poorly in others. If required, compare resampled
length vs original length after processing has finished. If data is noisy, may
benefit from sending smoothed signal to this algorithm and either using mapback
function or interpolating onto unsmoothed data.

If a third array, zᵢ, is provided then resample according to the points found by
the above process (which may be useful for resampling stress, strain and time with
one function call.)

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
- `zᵢ`: (Optional) z array. Array of data to interpolate and resample according to points found from y array information. (Usually non-prescribed variable e.g. strain or stress respectively)
- `pcntdownsample`: The ratio of new samples to old samples. Algorithm will halt ASAP once ratio greater than this number.
- `minperiod`: Minimum allowed distance between x array points, recommended to set > 0.0 to something like dx/10.0 to avoid algorithm over-focusing on a particular region.
- `minsamplenum = 25`: (Optional) number of initial, equally spaced seed samples required for algorithm to initialise.
"""
function var_resample(tᵢ::Array{Float64,1}, yᵢ::Array{Float64,1}, zᵢ::Array{Float64,1}, pcntdownsample::Float64, minperiod::Float64; minsamplenum::Int64 = 25)::Tuple{Array{Float64,1},Array{Float64,1}, Array{Float64,1}}

    # check array lengths are equal
    @assert length(tᵢ)==length(yᵢ) "X and Y arrays must have same length."
    @assert length(tᵢ)==length(zᵢ) "X and Z arrays must have same length."

    # length of original arrays
    originalSize = length(tᵢ)

    # interpolated, callable versions of y and z arrays. No need to interpolate
    # t as it is linear (fixed dt, monotonically increasing)
    yInterp = Interpolations.interpolate((tᵢ,), yᵢ, Interpolations.Gridded(Interpolations.Linear()))
    zInterp = Interpolations.interpolate((tᵢ,), zᵢ, Interpolations.Gridded(Interpolations.Linear()))

    # initialise arrays for resampled data
    xInit = zeros(Float64, minsamplenum)
    yInit = zeros(Float64, minsamplenum)

    # variables for generating correct intervals during initial sweep
    minsamplenum -= 1 # correction due integer rounding
    modNum = div(originalSize, minsamplenum)

    # initial sweep
    initCounter = 1
    @inbounds for i in 1:originalSize
        if mod(i-1, modNum)==0
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
    # reduces alpha and starts inner while loop sweeps again.
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

                end # end if
            end # end for loop
            iₒ += 1 # increment inner while loop counter
        end # inner while

        α -= 0.001 # decrement α for outer while loop

    end # outer while

    finalX = sort(collect(keys(xyDict)))
    finalY = [xyDict[i] for i in finalX]
    finalZ = [zInterp[i] for i in finalX]

    return finalX, finalY, finalZ
end

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
    mapback(xᵦ::Array{Float64,1},x::Array{Float64,1})

Match a variable downsampled array to its closest possible elements in original
array. Returns the array of unique indices corresponding to these matched elements.

Can be used after variable resampling to maintain resampling priorities whilst
not interfering with data of the fidelity by use of interpolations. Can also be
useful to ensure no region of the original data has been oversampled as duplicate
indices are deleted.
"""
function mapback(xᵦ::Array{Float64,1},x::Array{Float64,1})

    @assert length(xᵦ)<length(x) "xᵦ must be downsampled from x, i.e. have less elements."

    # array to store indices, do not know length a-priori so must mutate in place.
    indices = zeros(Int32,0)

    # loop through all values in resampled array, finding closest index in original array
    for i in 1:length(xᵦ)
        @inbounds iClosest = closestindex(x, xᵦ[i])
        append!(indices, iClosest)
    end

    # remove duplicates
    indices = unique(indices)

    # return indices, can be then used like xᵩ = x[indices]
    indices
end

"""
    downsample(boundaries::Array{Int64,1}, elPeriods::Array{Int64,1})

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
function downsample(boundaries::Array{Int64,1}, elPeriods::Array{Int64,1})

    # assert correct function signature
    @assert length(elPeriods)==length(boundaries)-1 "Number of different sample periods must be 1 less than boundaries provided"

    # initialise indices array
    indices = zeros(Int64,0)

    # loop through, skipping elements as required
    for i in 1:length(boundaries)-1

        # get indices for this 'section'
        indicesCur = boundaries[i]:elPeriods[i]:boundaries[i+1]

        # append
        append!(indices, collect(indicesCur))
    end

    # return unique indices
    indices = unique(indices)
end

"""
    fixed_resample(x::Array{Float64,1}, y::Array{Float64,1}[, z::Array{Float64,1}], boundaries::Array{Int64,1}, elPeriods::Array{Int64,1}, direction::Array{String,1})

Resample three (or two) arrays with new sample rate(s).

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
function fixed_resample(x::Array{Float64,1}, y::Array{Float64,1}, z::Array{Float64,1},
                        boundaries::Array{Int64,1}, elPeriods::Array{Int64,1},
                        direction::Array{String,1})

    # assert correct function signature
    @assert length(x)==length(y) "X and Y arrays must have same length."
    @assert length(x)==length(z) "X and Z arrays must have same length."
    @assert length(elPeriods)==length(boundaries)-1 "Number of different sample periods must be 1 less than boundaries provided"
    @assert length(elPeriods)==length(direction) "Number of sections must equal number of up/down regions"

    # y and z as callable interpolations, used for upsampled regions
    yInterp = Interpolations.interpolate((x,), y, Interpolations.Gridded(Interpolations.Linear()))
    zInterp = Interpolations.interpolate((x,), z, Interpolations.Gridded(Interpolations.Linear()))

    # initialise resampled arrays as empty
    xᵦ = zeros(Float64,0)
    yᵦ = zeros(Float64,0)
    zᵦ = zeros(Float64,0)

    # loop through boundaries
    for i in 1:length(boundaries)-1

        # upsampling, starts at each element then intepolates up N times
        if direction[i]=="up"
            for k in boundaries[i]:(boundaries[i+1]-1)
                # starting element
                append!(xᵦ,x[k])
                append!(yᵦ,yInterp[x[k]])
                append!(zᵦ,zInterp[x[k]])

                # dx increment up
                Δ = (x[k+1]-x[k])/elPeriods[i]
                for N in 1:(elPeriods[i]-1)
                    # add N increments to starting element
                    append!(xᵦ,x[k]+N*Δ)
                    append!(yᵦ,yInterp[x[k]+N*Δ])
                    append!(zᵦ,zInterp[x[k]+N*Δ])
                end
            end

        # downsampling, simply takes every N element as in downsample function
        elseif direction[i]=="down" #under/sampling
            append!(xᵦ,x[boundaries[i]:elPeriods[i]:(boundaries[i+1]-1)])
            append!(yᵦ,y[boundaries[i]:elPeriods[i]:(boundaries[i+1]-1)])
            append!(zᵦ,z[boundaries[i]:elPeriods[i]:(boundaries[i+1]-1)])
        end
    end

    # last element is missed out due to algorithm spec, needs to be added in after
    lastel = boundaries[end]
    # safety check then append
    if xᵦ[end]<x[lastel]
        append!(xᵦ,x[lastel])
        append!(yᵦ,y[lastel])
        append!(zᵦ,z[lastel])
    end

    return xᵦ, yᵦ, zᵦ
end

function fixed_resample(x::Array{Float64,1}, y::Array{Float64,1},
                        boundaries::Array{Int64,1}, elPeriods::Array{Int64,1},
                        direction::Array{String,1})

    # assert correct function signature
    @assert length(x)==length(y) "X and Y arrays must have same length."
    @assert length(elPeriods)==length(boundaries)-1 "Number of different sample periods must be 1 less than boundaries provided"

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
                Δ = (x[k+1]-x[k])/elPeriods[i]
                for N in 1:(elPeriods[i]-1)
                    # add N increments to starting element
                    append!(xᵦ,x[k]+N*Δ)
                    append!(yᵦ,yInterp[x[k]+N*Δ])
                end
            end

        # downsampling, simply takes every N element as in downsample function
        elseif direction[i]=="down" #under/sampling
            append!(xᵦ,x[boundaries[i]:elPeriods[i]:(boundaries[i+1]-1)])
            append!(yᵦ,y[boundaries[i]:elPeriods[i]:(boundaries[i+1]-1)])
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
