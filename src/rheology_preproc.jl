#!/usr/bin/env julia

"""
    var_resample(self::RheologyType, refvar::Symbol, pcntdownsample::Float64; _mapback::Bool = false)

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
function var_resample(self::RheologyType, refvar::Symbol, pcntdownsample::Float64; _mapback::Bool = false)

    # enforce minimum period of original period/10
    _minperiod = (self.t[2] - self.t[1])/10.0

    self_new = deepcopy(self)

    # get time resampled with respect to refvar
    (tᵦ, dummy)  = var_resample(self.t, getfield(self, refvar), pcntdownsample, _minperiod)

    for variable in self_new.numericdata
        # interpolate with respect to t
        var_interp = Interpolations.interpolate((self.t,), getfield(self, variable), Interpolations.Gridded(Interpolations.Linear()))

        # set field of new struct
        setfield!(self_new, variable, var_interp[tᵦ])
    end

    # if mapback required then mapback
    if _mapback
        # get mapped indices wrt original
        mapback_indices = mapback(tᵦ, self.t)
        
        for variable in self_new.numericdata

            setfield!(self_new, variable, getfield(self, variable)[mapback_indices])

        end
    end

    # plot if insight required
    if self.insight
        # plot all numeric data
        for variable in self_new.numericdata[1:3]
            # original data
            plot(self.t, getfield(self, variable), "o")
            # new data
            plot(self_new.t, getfield(self_new, variable), "x")
            title(String(variable))
            show()
        end
    end

    # change sampling type to variable
    self_new.sampling = "variable"
    
    # add record of operation applied
    self_new.appliedops = vcat(self_new.appliedops, "var_resample - refvar: $refvar, pcntdownsample: $pcntdownsample, _mapback: $_mapback")

    self_new

end

"""
    downsample(self::RheologyType, boundaries::Array{Int64,1}, elperiods::Array{Int64,1})

Use indices from `downsample` to downsample all data in the RheologyType struct.

See help docstring for `downsample` for more info on use of boundaries and elperiods
arguments.
"""
function downsample(self::RheologyType, boundaries::Array{Int64,1}, elperiods::Array{Int64,1})

    # get downsampled indices
    indices = downsample(boundaries::Array{Int64,1}, elperiods::Array{Int64,1})

    self_new = deepcopy(self)

    # downsample data
    for variable in self_new.numericdata
        # set field of new struct
        setfield!(self_new, variable, getfield(self, variable)[indices])
    end

    # plot if insight required
    if self.insight
        # plot all numeric data
        for variable in self_new.numericdata[1:3]
            # original data
            plot(self.t, getfield(self, variable), "o")
            # new data
            plot(self_new.t, getfield(self_new, variable), "x")
            title(String(variable))
            show()
        end
    end

    # change to variable sampling rate if more than one section
    if length(elperiods) > 1
        self_new.sampling = "variable"
    end

    # add record of operation applied
    appliedops = vcat(self.appliedops, "downsample - boundaries: $boundaries, elperiods: $elperiods")

    self_new

end

"""
    fixed_resample(self::RheologyType, boundaries::Array{Int64,1}, elperiods::Array{Int64,1}, direction::Array{String,1})

Resample three (or two) arrays with new sample rate(s).

Whereas `downsample` can only reduce the sample rate by not taking every array element,
fixed_resample can also upsample. This is why it cannot return indices as upsampled
sections will be interpolated versions of the original data.

See docstring for `fixed_resample` for more info and example on calling signature.
"""
function fixed_resample(self::RheologyType, boundaries::Array{Int64,1}, elperiods::Array{Int64,1}, direction::Array{String,1})

    self_new = deepcopy(self)

    # resample all data
    for variable in self_new.numericdata
        # resample derivatives
        (tᵦ, tempvarᵦ) = fixed_resample(self.t, getfield(self, variable), boundaries, elperiods, direction)

        setfield!(self_new, variable, tempvarᵦ)
    end

    # plot if insight required
    if self.insight
        # plot all numeric data
        for variable in self_new.numericdata[1:3]
            # original data
            plot(self.t, getfield(self, variable), "o")
            # new data
            plot(self_new.t, getfield(self_new, variable), "x")
            title(String(variable))
            show()
        end
    end

    # change to variable sampling rate if more than one section
    if length(elperiods) > 1
        self_new.sampling = "variable"
    end

    # add record of operation applied
    appliedops = vcat(self.appliedops, "fixed_resample - boundaries: $boundaries, elperiods: $elperiods, direction: $direction")

    self_new

end

"""
    smooth(self::RheologyType, τ::Float64; to_smooth = "all")

Smooth data using a Gaussian Kernel to time scale τ (approximately half power).

to_smooth can be an array of symbols so for example if you want to smooth ϵ and dϵ
then to_smooth = [:ϵ, :dϵ]
"""
function smooth(self::RheologyType, τ::Float64; to_smooth = "all")

    @assert self.sampling=="constant" "Sample rate must be constant for Gaussian smoothing kernel to function properly"

    # check what to smooth
    if to_smooth == "all"
        temp_variables = deepcopy(self.numericdata)
        filter!(i -> i≠:t, temp_variables)
        tosmooth = temp_variables
    else
        tosmooth = to_smooth
    end

    # copy self
    self_new = deepcopy(self)

    # get constant sample period
    dt = self.t[2] - self.t[1]

    # sample rate
    samplerate = 1.0/dt

    # loop through fields to smooth
    for i in tosmooth
        # set iteratively
        setfield!(self_new, i, smoothgauss(getfield(self, i), τ, samplerate))
    end

    # plot if insight required
    if self.insight
        # loop through symbols in tosmooth and plot
        for i in tosmooth[1:3]
            plot(self.t, getfield(self, i))
            plot(self_new.t, getfield(self_new, i), "--")
            title(String(i))
            show()
        end
    end

    # add record of operation applied
    self_new.appliedops = vcat(self.appliedops, "smooth - τ: $τ, to_smooth: $to_smooth")

    self_new

end

"""
    function mapbackdata(self_new::RheologyType, self_original::RheologyType)

Map back elements (WRT closest time elements) of all data from self_new to
self_original. See mapback help docstring for more info on how algorithm works.
"""
function mapbackdata(self_new::RheologyType, self_original::RheologyType)

    # get copy
    self_out = deepcopy(self_new)

    # symbols of data
    to_map = self_out.numericdata

    # get mapped back indices
    indices = mapback(self_new.t, self_original.t)

    # loop through to set variables
    for variable in to_map
        setfield!(self_out, variable, getfield(self_original, variable)[indices])
    end

    # plot if insight required
    if self_new.insight
        # loop through symbols in tosmooth and plot
        for variable in to_map[1:3]
            plot(self_original.t, getfield(self_original, variable), "o")
            plot(self_out.t, getfield(self_out, variable), "x")
            title(String(variable))
            show()
        end
    end

    # add record of operation applied
    self_out.appliedops = vcat(self_new.appliedops, "mapped back")

    self_out

end
