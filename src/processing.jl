#!/usr/bin/env julia

#############################
#~ Preprocessing Functions ~#
#############################

"""
    var_resample(self::RheologyData, refvar::Symbol, pcntdownsample::Float64; mapback::Bool = false)

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
function var_resample(self::RheologyData, refvar::Symbol, pcntdownsample::Float64; _mapback::Bool = false)

    # enforce minimum period of original period/10
    _minperiod = (self.t[2] - self.t[1])/10.0

    # get time resampled with respect to refvar
    (tᵦ, dummy)  = var_resample(self.t, getfield(self, refvar), pcntdownsample, _minperiod)

    local σ::Array{Float64,1}
    local ϵ::Array{Float64,1}
    local t::Array{Float64,1}

    if _mapback
        # get mapped indices wrt original
        mapback_indices = mapback(tᵦ, self.t)

        σ = self.σ[mapback_indices]
        ϵ = self.ϵ[mapback_indices]
        t = self.t[mapback_indices]
        
    elseif !_mapback
        # interpolate with respect to t
        σ_interp = Interpolations.interpolate((self.t,), self.σ, Interpolations.Gridded(Interpolations.Linear()))
        ϵ_interp = Interpolations.interpolate((self.t,), self.ϵ, Interpolations.Gridded(Interpolations.Linear()))
        t_interp = Interpolations.interpolate((self.t,), self.t, Interpolations.Gridded(Interpolations.Linear()))

        # resample all using new timepoints tᵦ
        σ = σ_interp[tᵦ]
        ϵ = ϵ_interp[tᵦ]
        t = t_interp[tᵦ]
    end

    # change sampling type to variable
    sampling = "variable"
    
    # add record of operation applied
    log = vcat(self.log, "var_resample - refvar: $refvar, pcntdownsample: $pcntdownsample, mapback: $_mapback")

    self_new = RheologyData(σ, ϵ, t, sampling, log)

end

"""
    downsample(self::RheologyData, boundaries::Array{Int64,1}, elperiods::Array{Int64,1})

Use indices from `downsample` to downsample all data in the RheologyData struct.

See help docstring for `downsample` for more info on use of boundaries and elperiods
arguments.
"""
function downsample(self::RheologyData, boundaries::Array{Int64,1}, elperiods::Array{Int64,1})

    # get downsampled indices
    indices = downsample(boundaries::Array{Int64,1}, elperiods::Array{Int64,1})

    # downsample data
    σ = self.σ[indices]
    ϵ = self.ϵ[indices]
    t = self.t[indices]

    # change to variable sampling rate if more than one section, if not then as original
    local sampling::String
    if length(elperiods) > 1
        sampling = "variable"
    else
        sampling = self.sampling
    end

    # add record of operation applied
    log = vcat(self.log, "downsample - boundaries: $boundaries, elperiods: $elperiods")

    self_new = RheologyData(σ, ϵ, t, sampling, log)

end

"""
    fixed_resample(self::RheologyData, boundaries::Array{Int64,1}, elperiods::Array{Int64,1}, direction::Array{String,1})

Resample three (or two) arrays with new sample rate(s).

Whereas `downsample` can only reduce the sample rate by not taking every array element,
fixed_resample can also upsample. This is why it cannot return indices as upsampled
sections will be interpolated versions of the original data.

See docstring for `fixed_resample` for more info and example on calling signature.
"""
function fixed_resample(self::RheologyData, boundaries::Array{Int64,1}, elperiods::Array{Int64,1}, direction::Array{String,1})

    # resample all data
    (t, σ) = fixed_resample(self.t, self.σ, boundaries, elperiods, direction)
    (t, ϵ) = fixed_resample(self.t, self.ϵ, boundaries, elperiods, direction)

    # change to variable sampling rate if more than one section
    local sampling::String
    if length(elperiods) > 1
        sampling = "variable"
    else
        sampling = self.sampling
    end

    # add record of operation applied
    log = vcat(self.log, "fixed_resample - boundaries: $boundaries, elperiods: $elperiods, direction: $direction")

    self_new = RheologyData(σ, ϵ, t, sampling, log)

end

"""
    smooth(self::RheologyData, τ::Float64)

Smooth data using a Gaussian Kernel to time scale τ (approximately half power).

Smooths both σ and ϵ.
"""
function smooth(self::RheologyData, τ::Float64)

    @assert self.sampling=="constant" "Sample rate must be constant for Gaussian smoothing kernel to function properly"

    t = self.t

    # get constant sample period
    dt = t[2] - t[1]

    # sample rate
    samplerate = 1.0/dt;

    σ = smoothgauss(self.σ, τ, samplerate)
    ϵ = smoothgauss(self.ϵ, τ, samplerate)

    # add record of operation applied
    log = vcat(self.log, "smooth - τ: $τ")

    self_new = RheologyData(σ, ϵ, t, self.sampling, log)

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

    self_new = RheologyData(σ, ϵ, t, self_new.sampling, log)

end

##########################
#~ Processing Functions ~#
##########################

"""
    modelfit(data::RheologyData, modulus::Function[, p0::Array{Float64,1}, lo::Array{Float64,1}, hi::Array{Float64,1}; verbose::Bool = false])

Fit RheologyData struct to model and return a fitted model as a RheologyModel object.

# Arguments

- `data`: RheologyData struct containing all data
- `modulus`: E.g. G_SLS, J_springpot etc. See base_models.jl for full list
- `p0`: Initial parameters to use in fit
- `lo`: Lower bounds for parameters
- `hi`: Higher bounds for parameters
- `verbose`: If true, prints parameters on each optimisation iteration
"""
function modelfit(data::RheologyData,
                  modulus::Function;
                  p0::Array{Float64,1} = [-1.0],
                  lo::Array{Float64,1} = [-1.0],
                  hi::Array{Float64,1} = [-1.0],
                  verbose::Bool = false)::RheologyModel

    # get modulus info from database
    (controlledvar, p0_default) = modeldatabase(modulus)

    # get singularity presence
    sing = singularitytest(modulus)

    # use default p0 if not provided
    if quasinull(p0)
        p0 = p0_default
    end 

    # generate time series difference array (for convolution)
    dt_series = deriv(data.t)

    # get derivative of controlled variable and measured variable
    if controlledvar == "σ"
        dcontrolled = deriv(data.σ, data.t)
        measured = data.ϵ
    elseif controlledvar == "ϵ"
        dcontrolled = deriv(data.ϵ, data.t)
        measured = data.σ
    end

    # start fit
    tic()
    (minf, minx, ret) = leastsquares_init(p0,
                                          lo, 
                                          hi,
                                          modulus, 
                                          data.t, 
                                          dt_series, 
                                          dcontrolled,
                                          measured; 
                                          insight = verbose,
                                          sampling = data.sampling, 
                                          singularity = sing)
    timetaken = toq()

    #
    modulusname = string(:modulus)

    log = vcat(data.log, "Fitted modulus($modulusname) in $timetaken, finished due to $ret with parans $minx")

    RheologyModel(modulus, minx, log)

end

"""
    modelpredict(data::RheologyData, model::RheologyModel)

Given data and model and parameters, predict new dataset based on both.
"""
function modelpredict(data::RheologyData, model::RheologyModel)::RheologyData

    # get modulus info from database
    (controlledvar, p0_default) = modeldatabase(model.modulus)

    # get singularity presence
    sing = singularitytest(model.modulus)

    if controlledvar == "σ"
        dcontrolled = deriv(data.σ, data.t)
    elseif controlledvar == "ϵ"
        dcontrolled = deriv(data.ϵ, data.t)
    end

    # generate time series difference array (for convolution)
    dt_series = deriv(data.t)

    # get convolution
    if !sing && data.sampling == "constant"
        convolved = boltzconvolve_nonsing(model.modulus, data.t, deriv(data.t), model.parameters, dcontrolled)

    elseif sing && data.sampling == "constant"
        convolved = boltzconvolve_sing(model.modulus, data.t, deriv(data.t), model.parameters, dcontrolled)

    elseif !sing && data.sampling == "variable"
        convolved = boltzintegral_nonsing(model.modulus, data.t, model.parameters, dcontrolled)

    elseif sing && data.sampling == "variable"
        convolved = boltzintegral_sing(model.modulus, data.t, model.parameters, dcontrolled)

    end

    if controlledvar == "σ"
        σ = data.σ
        ϵ = convolved
    elseif controlledvar == "ϵ"
        σ = convolved
        ϵ = data.ϵ
    end   

    # store operation
    log = vcat(data.log, "Predicted data from model:", model.log)

    RheologyData(σ, ϵ, data.t, data.sampling, log)

end

##############################
#~ Postprocessing Functions ~#
##############################

"""
    fiteval(self::RheologyData[, modelname::String]; singularity = false)

Show plot of data vs. fitted data for specified model. If no model specified,
shows plot with all fitted models.
"""
function fiteval(self::RheologyData, modelname::String)

    # params
    params = self.fittedmodels[modelname][2]

    # modulus function
    model = moduli(modelname, self.test_type)

    # get fit
    if !model.singularity && self.sampling == "constant"
        fitted = boltzconvolve_nonsing(model, self.t, deriv(self.t), params, self.dcontrolled)

    elseif model.singularity && self.sampling == "constant"
        fitted = boltzconvolve_sing(model, self.t, deriv(self.t), params, self.dcontrolled)

    elseif !model.singularity && self.sampling == "variable"
        fitted = boltzintegral_nonsing(model, self.t, params, self.dcontrolled)

    elseif model.singularity && self.sampling == "variable"
        fitted = boltzintegral_sing(model, self.t, params, self.dcontrolled)

    end

    # print params
    println(modelname, " fit: ", self.fittedmodels[modelname])

    # special for AFMData, Hertz sphere
    if typeof(self) == AFMData && self.test_type == "strlx"
        measured = self.measured.^(2/3)
        fitted = fitted.^(2/3)
    else
        measured = self.measured
    end

    if model.singularity
        plot(self.t[1:end], measured)
        plot(self.t[2:end], fitted, "--")
        show()
    else
        plot(self.t, measured)
        plot(self.t, fitted, "--")
        show()
    end

end

function fiteval(self::RheologyData)

    # special for AFMData, Hertz sphere
    if typeof(self) == AFMData && self.test_type == "strlx"
        measured = self.measured.^(2/3)
    else
        measured = self.measured
    end

    fig, ax = subplots()

    ax[:plot](self.t, measured, label="Original Data", alpha=0.65, linewidth=2)

    for (modelname, value) in self.fittedmodels

        # params
        params = self.fittedmodels[modelname][2]

        # modulus function
        model = moduli(modelname, self.test_type)

        # get fit
        if !model.singularity && self.sampling == "constant"
            fitted = boltzconvolve_nonsing(model, self.t, deriv(self.t), params, self.dcontrolled)

        elseif model.singularity && self.sampling == "constant"
            fitted = boltzconvolve_sing(model, self.t, deriv(self.t), params, self.dcontrolled)

        elseif !model.singularity && self.sampling == "variable"
            fitted = boltzintegral_nonsing(model, self.t, params, self.dcontrolled)

        elseif model.singularity && self.sampling == "variable"
            fitted = boltzintegral_sing(model, self.t, params, self.dcontrolled)

        end

        # print params
        println(modelname, " fit: ", self.fittedmodels[modelname])

        # special for AFMData, Hertz sphere
        if typeof(self) == AFMData && self.test_type == "strlx"
            fitted = fitted.^(2/3)
        end

        if model.singularity
            ax[:plot](self.t[2:end], fitted, label=modelname, alpha=0.85, linewidth=2)
        else
            ax[:plot](self.t, fitted, label=modelname, alpha=0.85, linewidth=2)
        end

    end

    ax[:legend](loc="best")

    show()

end

"""
    saveresult(self::RheologyData; include_data::Bool = false)

Save RheologyData object using JLD format. If include_data set to false
(by default) then :σ, :ϵ, :t, :dσ, :dϵ fields are set to [0.0] to save
space on disk. If include_data is set as true then these fields are saved
as is.
"""
function saveresult(self::RheologyData; include_data::Bool = false)

    # include original/processed numerical data σ, ϵ, t...etc.
    if include_data

        # save
        jldopen(string(self.filedir[1:end-4], "_RheologyData.jld"), "w") do file
            # register RHEOS module (with RheologyData type) to JLD
            addrequire(file, RHEOS)
            # write self to file
            write(file, "self", self)
        end
    # or not    
    elseif !include_data

        # get copy
        self_copy = deepcopy(self)

        # member variables to erase
        to_reset = self.numericdata
        for n in to_reset
            setfield!(self_copy, n, [0.0])
        end

        # save 
        jldopen(string(self_copy.filedir[1:end-4], "_RheologyMetadata.jld"), "w") do file
            # register RHEOS module (with RheologyData type) to JLD
            addrequire(file, RHEOS)
            # write self_copy to file
            write(file, "self", self_copy)
        end
    end
end

"""
    loadresult(filedir::String)

Convenience function loads result without having to call loadresult(filedir)["self"]
"""
function loadresult(filedir::String)
    # load in result
    loaded_result = load(filedir)
    # return
    loaded_result["self"]
end

