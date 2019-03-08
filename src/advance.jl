# Processing

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
function variableresample(self::RheoTimeData, refvar::Symbol, pcntdownsample::Real; _mapback::Bool = true)

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

###############################################################################
#           BASE
###############################################################################


function sampleratecompare(t1::Vector{RheoFloat}, t2::Vector{RheoFloat})

    @assert constantcheck(t1) "Sample-rate of both arguments must be constant, first argument is non-constant"
    @assert constantcheck(t2) "Sample-rate of both arguments must be constant, second argument is non-constant"

    diff1 = getsampleperiod(t1)
    diff2 = getsampleperiod(t2)

    diff1 == diff2

end

"""
    var_resample(tᵢ, yᵢ, pcntdownsample, minperiod[; minsamplenum = 25])

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
function var_resample(tᵢ::Vector{T}, yᵢ::Vector{T}, pcntdownsample::T, minperiod::T; minsamplenum::Integer = 25) where T<:Real

    #@eval using Interpolations

    @assert length(tᵢ)==length(yᵢ) "X and Y arrays must have same length."

    # length of original arrays
    originalSize = length(tᵢ)

    # interpolated, callable versions of y and z arrays. No need to interpolate
    # t as it is linear (fixed dt, monotonically increasing)
    yInterp = Base.invokelatest(interpolate, (tᵢ,), yᵢ, Base.invokelatest(Gridded, Base.invokelatest(Linear)))

    # initialise arrays for resampled data
    xInit = zeros(T, minsamplenum + 1)
    yInit = zeros(T, minsamplenum + 1)

    # variables for generating correct intervals during initial sweep
    modNum = div(originalSize, minsamplenum)

    # initial sweep
    initCounter = 1
    # @inbounds for i in 1:originalSize
    for i in 1:originalSize
        if mod(i-1, modNum) == 0
            xInit[initCounter] = tᵢ[i]
            yInit[initCounter] = yᵢ[i]
            initCounter += 1
        end
    end

    # send initial resampled arrays to a dict for more convenient mutability
    xyDict = Dict{T, T}(xInit[i] => yInit[i] for i in 1:length(xInit))

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

            dy = derivBD(y, x)
            dyNorm = dy/maximum(abs.(dy))

            ddy = derivBD(dy, x)
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
    mapback(x₁, x₀)

Match a resampled array x₁ to its closest possible elements in original
array x₀. Returns the array of unique indices corresponding to these matched elements.

Can be used after variable resampling to maintain resampling priorities whilst
not interfering with data of the fidelity by use of interpolations. Can also be
useful to ensure no region of the original data has been oversampled as duplicate
indices are deleted.
"""
function mapback(x₁::Vector{T}, x₀::Vector{T}; return_indices=true) where T<:Real

    # get indices
    indices = closestindices(x₀, x₁)

    # remove duplicates
    indices = unique(indices)

    if return_indices
        return indices
    else
        return x₀[indices]
    end
end

function downsample(boundaries::Vector{T}, elperiods::Vector{T}) where T<:Integer

    # assert correct function signature
    @assert length(elperiods)==length(boundaries)-1 "Number of different sample periods must be 1 less than boundaries provided"

    # initialise indices array
    indices = zeros(Integer, 0)

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
