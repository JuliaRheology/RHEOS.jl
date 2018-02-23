#!/usr/bin/env julia

"""
    var_resample!(self::RheologyData, refvar::Symbol, pcntdownsample::Float64; _mapback = true)

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
function var_resample!(self::RheologyData, refvar::Symbol, pcntdownsample::Float64; _mapback = false)

    # set sampling rate as variable
    self.sampling = "variable"

    # enforce minimum period of original period/10
    _minperiod = (self.t[2] - self.t[1])/10.0

    # interpolate all with respect to t
    ϵ_interp = Interpolations.interpolate((self.tᵦ,), self.ϵᵦ, Interpolations.Gridded(Interpolations.Linear()))
    σ_interp = Interpolations.interpolate((self.tᵦ,), self.σᵦ, Interpolations.Gridded(Interpolations.Linear()))
    dϵ_interp = Interpolations.interpolate((self.tᵦ,), self.dϵᵦ, Interpolations.Gridded(Interpolations.Linear()))
    dσ_interp = Interpolations.interpolate((self.tᵦ,), self.dσᵦ, Interpolations.Gridded(Interpolations.Linear()))

    # get time resampled with respect to refvar
    (self.tᵦ, dummy)  = var_resample(self.tᵦ, getfield(self, refvar), pcntdownsample, _minperiod)

    # get resampled
    self.ϵᵦ = ϵ_interp[self.tᵦ]
    self.σᵦ = σ_interp[self.tᵦ]
    self.dϵᵦ = dϵ_interp[self.tᵦ]
    self.dσᵦ = dσ_interp[self.tᵦ]

    # if mapback required then mapback
    if _mapback
        # get mapped indices wrt original
        mapback_indices = mapback(self.tᵦ, self.t)
        #map
        self.tᵦ = self.t[mapback_indices]
        self.ϵᵦ = self.ϵ[mapback_indices]
        self.σᵦ = self.σ[mapback_indices]
        self.dϵᵦ = self.dϵ[mapback_indices]
        self.dσᵦ = self.dσ[mapback_indices]
    end

    # plot if insight required
    if self.insight
        # strain subplot
        subplot(121)
        plot(self.t, self.ϵ, "o")
        plot(self.tᵦ, self.ϵᵦ, "x")
        title("Strain")

        # stress subplot
        subplot(122)
        plot(self.t, self.σ, "o")
        plot(self.tᵦ, self.σᵦ, "x")
        title("Stress")

        show()
    end

end

"""
    downsample!(self::RheologyData, boundaries::Array{Int64,1}, elPeriods::Array{Int64,1})

Use indices from `downsample` to downsample all data in the RheologyData struct.

See help docstring for `downsample` for more info on use of boundaries and elPeriods
arguments.
"""
function downsample!(self::RheologyData, boundaries::Array{Int64,1}, elPeriods::Array{Int64,1})

    # change to variable sampling rate if more than one section
    if length(elPeriods) > 1
        self.sampling = "variable"
    end

    # get downsampled indices
    indices = downsample(boundaries::Array{Int64,1}, elPeriods::Array{Int64,1})

    # downsample data
    self.tᵦ = self.t[indices]
    self.σᵦ = self.σ[indices]
    self.ϵᵦ = self.ϵ[indices]
    self.dϵᵦ = self.dϵ[indices]
    self.dσᵦ = self.dσ[indices]

    # plot if insight required
    if self.insight
        # strain subplot
        subplot(121)
        plot(self.t, self.ϵ)
        plot(self.tᵦ, self.ϵᵦ, "x")
        title("Strain")

        # stress subplot
        subplot(122)
        plot(self.t, self.σ)
        plot(self.tᵦ, self.σᵦ, "x")
        title("Stress")

        show()
    end

end

"""
    fixed_resample!(self::RheologyData, boundaries::Array{Int64,1}, elPeriods::Array{Int64,1}, direction::Array{String,1})

Resample three (or two) arrays with new sample rate(s).

Whereas `downsample` can only reduce the sample rate by not taking every array element,
fixed_resample can also upsample. This is why it cannot return indices as upsampled
sections will be interpolated versions of the original data.

See docstring for `fixed_resample` for more info and example on calling signature.
"""
function fixed_resample!(self::RheologyData, boundaries::Array{Int64,1}, elPeriods::Array{Int64,1}, direction::Array{String,1})

    # change to variable sampling rate if more than one section
    if length(elPeriods) > 1
        self.sampling = "variable"
    end

    # temp variable, think of a better way round this
    dummyt = self.tᵦ

    # resample main data
    (self.tᵦ, self.σᵦ, self.ϵᵦ) = fixed_resample(dummyt, self.σᵦ, self.ϵᵦ,
                                                 boundaries, elPeriods, direction)

    # resample derivatives
    (self.tᵦ, self.dϵᵦ) = fixed_resample(dummyt, self.dϵᵦ, boundaries,
                                          elPeriods, direction)

    (self.tᵦ, self.dσᵦ) = fixed_resample(dummyt, self.dσᵦ, boundaries,
                                          elPeriods, direction)

    # plot if insight required
    if self.insight
        # strain subplot
        subplot(121)
        plot(self.t, self.ϵ)
        plot(self.tᵦ, self.ϵᵦ, "x")
        title("Strain")

        # stress subplot
        subplot(122)
        plot(self.t, self.σ)
        plot(self.tᵦ, self.σᵦ, "x")
        title("Stress")

        show()
    end
end
