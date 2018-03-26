#!/usr/bin/env julia

function var_resample(self::AFMData, refvar::Symbol, pcntdownsample::Float64; _mapback = false)

    # enforce minimum period of original period/10
    _minperiod = (self.t[2] - self.t[1])/10.0

    # interpolate all with respect to t
    f_interp = Interpolations.interpolate((self.t,), self.f, Interpolations.Gridded(Interpolations.Linear()))
    δ_interp = Interpolations.interpolate((self.t,), self.δ, Interpolations.Gridded(Interpolations.Linear()))
    cm_measured_interp = Interpolations.interpolate((self.t,), self.cm_measured, Interpolations.Gridded(Interpolations.Linear()))
    cm_dcontrolled_interp = Interpolations.interpolate((self.t,), self.cm_dcontrolled, Interpolations.Gridded(Interpolations.Linear()))

    # get time resampled with respect to refvar
    (tᵦ, dummy)  = var_resample(self.t, getfield(self, refvar), pcntdownsample, _minperiod)

    # get resampled
    fᵦ = f_interp[tᵦ]
    δᵦ = δ_interp[tᵦ]
    cm_measuredᵦ = cm_measured_interp[tᵦ]
    cm_dcontrolledᵦ = cm_dcontrolled_interp[tᵦ]

    # if mapback required then mapback
    if _mapback
        # get mapped indices wrt original
        mapback_indices = mapback(tᵦ, self.t)
        #map
        tᵦ = self.t[mapback_indices]
        δᵦ = self.δ[mapback_indices]
        fᵦ = self.f[mapback_indices]
        cm_measuredᵦ = self.cm_measured[mapback_indices]
        cm_dcontrolledᵦ = self.cm_dcontrolled[mapback_indices]
    end

    # plot if insight required
    if self.insight
        # strain subplot
        subplot(121)
        plot(self.t, self.δ, "o")
        plot(tᵦ, δᵦ, "x")
        title("Displacement")
        # stress subplot
        subplot(122)
        plot(self.t, self.f, "o")
        plot(tᵦ, fᵦ, "x")
        title("Force")
        # show
        show()
    end

    # add record of operation applied
    appliedops = vcat(self.appliedops, "var_resample - refvar: $refvar, pcntdownsample: $pcntdownsample, _mapback: $_mapback")

    # construct new AFMData struct with resampled members
    AFMData(self.filedir, self.insight, "variable", self.test_type, self.R, fᵦ, δᵦ, tᵦ, cm_measuredᵦ, cm_dcontrolledᵦ, Dict(), appliedops)

end

function downsample(self::AFMData, boundaries::Array{Int64,1}, elperiods::Array{Int64,1})

    # get downsampled indices
    indices = downsample(boundaries::Array{Int64,1}, elperiods::Array{Int64,1})

    # downsample data
    tᵦ = self.t[indices]
    δᵦ = self.δ[indices]
    fᵦ = self.f[indices]
    cm_measuredᵦ = self.cm_measured[indices]
    cm_dcontrolledᵦ = self.cm_dcontrolled[indices]

    # plot if insight required
    if self.insight
        # strain subplot
        subplot(121)
        plot(self.t, self.δ)
        plot(tᵦ, δᵦ, "x")
        title("Displacement")
        # stress subplot
        subplot(122)
        plot(self.t, self.f)
        plot(tᵦ, fᵦ, "x")
        title("Force")
        # show
        show()
    end

    # change to variable sampling rate if more than one section
    if length(elperiods) > 1
        _sampling = "variable"
    # otherwise, remain as before
    else
        _sampling = self.sampling
    end

    # add record of operation applied
    appliedops = vcat(self.appliedops, "downsample - boundaries: $boundaries, elperiods: $elperiods")

    # construct new AFMData struct with resampled members
    AFMData(self.filedir, self.insight, _sampling, self.test_type, self.R, fᵦ, δᵦ, tᵦ, cm_measuredᵦ, cm_dcontrolledᵦ, Dict(), appliedops)

end

function fixed_resample(self::AFMData, boundaries::Array{Int64,1}, elperiods::Array{Int64,1}, direction::Array{String,1})

    # resample main data
    (tᵦ, fᵦ, δᵦ) = fixed_resample(self.t, self.f, self.δ,
                                                 boundaries, elperiods, direction)

    # resample derivatives
    (tᵦ, cm_measuredᵦ) = fixed_resample(self.t, self.cm_measured, boundaries,
                                          elperiods, direction)

    (tᵦ, cm_dcontrolledᵦ) = fixed_resample(self.t, self.cm_dcontrolled, boundaries,
                                          elperiods, direction)

    # plot if insight required
    if self.insight
        # strain subplot
        subplot(121)
        plot(self.t, self.δ)
        plot(tᵦ, δᵦ, "x")
        title("Strain")
        # stress subplot
        subplot(122)
        plot(self.t, self.f)
        plot(tᵦ, fᵦ, "x")
        title("Stress")
        # show
        show()
    end

    # change to variable sampling rate if more than one section
    if length(elperiods) > 1
        _sampling = "variable"
    # otherwise, remain as before
    else
        _sampling = self.sampling
    end

    # add record of operation applied
    appliedops = vcat(self.appliedops, "fixed_resample - boundaries: $boundaries, elperiods: $elperiods, direction: $direction")

    # construct new AFMData struct with resampled members
    AFMData(self.filedir, self.insight, _sampling, self.test_type, self.R, fᵦ, δᵦ, tᵦ, cm_measuredᵦ, cm_dcontrolledᵦ, Dict(), appliedops)

end

function smooth(self::AFMData, τ::Float64; to_smooth = "all")

    @assert self.sampling=="constant" "Sample rate must be constant for Gaussian smoothing kernel to function properly"

    # check what to smooth
    if to_smooth == "all"
        tosmooth = [:f, :δ, :cm_measured, :cm_dcontrolled]
    else
        tosmooth = to_smooth
    end

    # copy self
    self_copy = deepcopy(self)

    # get constant sample period
    dt = self.t[2] - self.t[1]

    # sample rate
    samplerate = 1.0/dt

    # loop through fields to smooth
    for i in tosmooth
        # set iteratively
        setfield!(self_copy, i, smoothgauss(getfield(self, i), τ, samplerate))
    end

    # plot if insight required
    if self.insight
        # loop through symbols in tosmooth and plot
        for i in tosmooth
            plot(self.t, getfield(self, i))
            plot(self_copy.t, getfield(self_copy, i), "--")
            title(String(i))
            show()
        end
    end

    # add record of operation applied
    self_copy.appliedops = vcat(self.appliedops, "smooth - τ: $τ, to_smooth: $to_smooth")

    # return
    self_copy
end

function mapbackdata(self_new::AFMData, self_original::AFMData)

    # get copy
    self_out = deepcopy(self_new)

    # symbols of data
    to_map = [:f, :δ, :cm_measured, :cm_dcontrolled]

    # get mapped back indices
    indices = mapback(self_new.t, self_original.t)

    # loop through to set variables
    for i in to_map
        setfield!(self_out, i, getfield(self_original, i)[indices])
    end

    # plot if insight required
    if self_new.insight
        # loop through symbols in tosmooth and plot
        for i in to_map
            plot(self_original.t, getfield(self_original, i), "o")
            plot(self_out.t, getfield(self_out, i), "x")
            title(String(i))
            show()
        end
    end

    # add record of operation applied
    self_out.appliedops = vcat(self_new.appliedops, "mapped back")

    # return
    self_out
end
