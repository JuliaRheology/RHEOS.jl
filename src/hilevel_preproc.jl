#!/usr/bin/env julia

"""
    var_resample(self::RheologyData, refvar::Symbol, pcntdownsample::Float64; _mapback = true)

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
function var_resample(self::RheologyData, refvar::Symbol, pcntdownsample::Float64; _mapback = false)

    # enforce minimum period of original period/10
    _minperiod = (self.t[2] - self.t[1])/10.0

    # interpolate all with respect to t
    ϵ_interp = Interpolations.interpolate((self.t,), self.ϵ, Interpolations.Gridded(Interpolations.Linear()))
    σ_interp = Interpolations.interpolate((self.t,), self.σ, Interpolations.Gridded(Interpolations.Linear()))
    dϵ_interp = Interpolations.interpolate((self.t,), self.dϵ, Interpolations.Gridded(Interpolations.Linear()))
    dσ_interp = Interpolations.interpolate((self.t,), self.dσ, Interpolations.Gridded(Interpolations.Linear()))

    # get time resampled with respect to refvar
    (tᵦ, dummy)  = var_resample(self.t, getfield(self, refvar), pcntdownsample, _minperiod)

    # get resampled
    ϵᵦ = ϵ_interp[tᵦ]
    σᵦ = σ_interp[tᵦ]
    dϵᵦ = dϵ_interp[tᵦ]
    dσᵦ = dσ_interp[tᵦ]

    # if mapback required then mapback
    if _mapback
        # get mapped indices wrt original
        mapback_indices = mapback(tᵦ, self.t)
        #map
        tᵦ = self.t[mapback_indices]
        ϵᵦ = self.ϵ[mapback_indices]
        σᵦ = self.σ[mapback_indices]
        dϵᵦ = self.dϵ[mapback_indices]
        dσᵦ = self.dσ[mapback_indices]
    end

    # plot if insight required
    if self.insight
        # strain subplot
        subplot(121)
        plot(self.t, self.ϵ, "o")
        plot(tᵦ, ϵᵦ, "x")
        title("Strain")
        # stress subplot
        subplot(122)
        plot(self.t, self.σ, "o")
        plot(tᵦ, σᵦ, "x")
        title("Stress")
        # show
        show()
    end

    # construct new RheologyData struct with resampled members
    RheologyData(self.insight, "variable", self.test_type, σᵦ, ϵᵦ, tᵦ, dσᵦ, dϵᵦ)

end

"""
    downsample(self::RheologyData, boundaries::Array{Int64,1}, elPeriods::Array{Int64,1})

Use indices from `downsample` to downsample all data in the RheologyData struct.

See help docstring for `downsample` for more info on use of boundaries and elPeriods
arguments.
"""
function downsample(self::RheologyData, boundaries::Array{Int64,1}, elPeriods::Array{Int64,1})

    # get downsampled indices
    indices = downsample(boundaries::Array{Int64,1}, elPeriods::Array{Int64,1})

    # downsample data
    tᵦ = self.t[indices]
    σᵦ = self.σ[indices]
    ϵᵦ = self.ϵ[indices]
    dϵᵦ = self.dϵ[indices]
    dσᵦ = self.dσ[indices]

    # plot if insight required
    if self.insight
        # strain subplot
        subplot(121)
        plot(self.t, self.ϵ)
        plot(tᵦ, ϵᵦ, "x")
        title("Strain")
        # stress subplot
        subplot(122)
        plot(self.t, self.σ)
        plot(tᵦ, σᵦ, "x")
        title("Stress")
        # show
        show()
    end

    # change to variable sampling rate if more than one section
    if length(elPeriods) > 1
        _sampling = "variable"
    # otherwise, remain as before
    else
        _sampling = self.sampling
    end

    # construct new RheologyData struct with resampled members
    RheologyData(self.insight, _sampling, self.test_type, σᵦ, ϵᵦ, tᵦ, dσᵦ, dϵᵦ)

end

"""
    fixed_resample(self::RheologyData, boundaries::Array{Int64,1}, elPeriods::Array{Int64,1}, direction::Array{String,1})

Resample three (or two) arrays with new sample rate(s).

Whereas `downsample` can only reduce the sample rate by not taking every array element,
fixed_resample can also upsample. This is why it cannot return indices as upsampled
sections will be interpolated versions of the original data.

See docstring for `fixed_resample` for more info and example on calling signature.
"""
function fixed_resample(self::RheologyData, boundaries::Array{Int64,1}, elPeriods::Array{Int64,1}, direction::Array{String,1})

    # resample main data
    (tᵦ, σᵦ, ϵᵦ) = fixed_resample(self.t, self.σ, self.ϵ,
                                                 boundaries, elPeriods, direction)

    # resample derivatives
    (tᵦ, dϵᵦ) = fixed_resample(self.t, self.dϵ, boundaries,
                                          elPeriods, direction)

    (tᵦ, dσᵦ) = fixed_resample(self.t, self.dσ, boundaries,
                                          elPeriods, direction)

    # plot if insight required
    if self.insight
        # strain subplot
        subplot(121)
        plot(self.t, self.ϵ)
        plot(tᵦ, ϵᵦ, "x")
        title("Strain")
        # stress subplot
        subplot(122)
        plot(self.t, self.σ)
        plot(tᵦ, σᵦ, "x")
        title("Stress")
        # show
        show()
    end

    # change to variable sampling rate if more than one section
    if length(elPeriods) > 1
        _sampling = "variable"
    # otherwise, remain as before
    else
        _sampling = self.sampling
    end

    # construct new RheologyData struct with resampled members
    RheologyData(self.insight, _sampling, self.test_type, σᵦ, ϵᵦ, tᵦ, dσᵦ, dϵᵦ)

end
