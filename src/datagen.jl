#!/usr/bin/env julia

#=
----------------------------------------------------------------------------
timeline generation function forms the basis of any load-generation workflow
----------------------------------------------------------------------------
=#

"""
    timeline(r)

Generate `RheoTimeData` struct with only the time data based on the range provided.

# Example:

`timeline(0:0.1:10)` returns a RheoTimeData with only time data, with the following series: [0.  0.1   0.2  ...  9.9  10].
"""
function timeline(r::R; savelog = true) where R <: AbstractRange
    log = loginit(savelog,:timeline, params=(r,), info=(comment="timeline created", type=time_only))
    RheoTimeData([], [], Vector{RheoFloat}(r), log)
end

"""
    timeline(;t_start::Real=0., t_end::Real, step::Real=(t_end - t_start)/250.)

Generate `RheoTimeData` struct with only the time data based on the keywords provided.

# Arguments

- `t_start`: Starting time, typically 0
- `t_end`: End time
- `step`: Time between sample. Default set in order to provide ≈ 250 time points.
"""
function timeline(;t_start::Real = 0., t_end::Real, step::Real = (t_end - t_start)/250., savelog = true)
    timeline(t_start:step:t_end, savelog = savelog)
end



"""
    frequencyspec(r::R; logscale=true)

Generate `RheoFreqData` struct with only the frequency data distributed on a linear or log scale.

# Arguments

- `r`: range, e.g. 1:0.1:100
- `logscale`: if true, set the range in logscale, e.g. -2:0.1:20  --> 10 values per decade between 10^-2 and 10^2
"""
function frequencyspec(r::R; logscale=true, savelog = true) where R <: AbstractRange
    if logscale
        ω=[10^o for o in r]
    else
        ω=Vector{RheoFloat}(r)
    end
    log = loginit(savelog, :frequencyspec, params=(r,), keywords=(logscale=logscale,), info=(comment="frequency range created", type=freq_only))
    RheoFreqData([], [], ω, log)
end


"""
    frequencyspec(;ω_start::Real, ω_end::Real, logstep::Real= 0.05)

Generate `RheoFreqData` struct with only the frequency data distributed on a log scale.

# Arguments

- `ω_start`: min frequency
- `ω_end`: max frequency
- `logstep`: step between frequencies in log scale, e.g. 0.1 --> 10 values per decade.
"""
function frequencyspec(;ω_start::Real, ω_end::Real, logstep::Real = 0.05, savelog = true)
    frequencyspec(log10(ω_start):logstep:log10(ω_end), logscale=true, savelog = savelog)
end

#=
------------------------------------------------------------------------
strainfunction and stressfunction can be used to impose strain or stress
respectively. The subsequent convenience functions are passed to them to
generate the appropriate loading.
------------------------------------------------------------------------
=#
"""
    strainfunction(data::RheoTimeData, f::T) where T<:Function
    strainfunction(f::T, data::RheoTimeData) where T<:Function

Accepts a `RheoTimeData` and outputs a new `RheoTimeData` with the strain imposed.
The strain signal is determined by the function provided, which should take
time as its only argument. The original data's time signal is used.

Normally used with a RheoTimeData generated using the `timeline` function.
"""
function strainfunction(data::RheoTimeData, f::T) where T<:Function
    log = logadd_process(data, :strainfunction, params=(f=f,), comment="Strain function applied to timeline" )
    return RheoTimeData(data.σ, rheoconvert(map(f, data.t)), data.t, log)
end

function strainfunction(f::T, data::RheoTimeData) where T<:Function
    strainfunction(data,f)
end

"""
    strainfunction!(data::RheoTimeData, f::T) where T<:Function
    strainfunction!(f::T, data::RheoTimeData) where T<:Function

In-place version of `strainfunction`.
"""
function strainfunction!(data::RheoTimeData, f::T) where T<:Function
    logadd_process!(data, :strainfunction, params=(f=f,), comment="Strain function applied to timeline" )
    _mapdata!(f,data.ϵ,data.t) 
    return data
end

function strainfunction!(f::T, data::RheoTimeData) where T<:Function
    strainfunction!(data,f)
end

"""
    stressfunction(data::RheoTimeData, f::T) where T<:Function
    stressfunction(f::T, data::RheoTimeData) where T<:Function

Accepts a `RheoTimeData` and outputs a new `RheoTimeData` with a stress imposed.
The stress signal is determined by the function provided, which should take
time as its only argument. The original data's time signal is used.

Normally used with a `RheoTimeData` generated using the `timeline` function.
"""
function stressfunction(data::RheoTimeData, f::T) where T<:Function
    log = logadd_process(data, :stressfunction, params=(f=f,), comment="Stress function applied to timeline" )
    return RheoTimeData(convert(Vector{RheoFloat}, map(f, data.t)), data.ϵ, data.t, log)
end

function stressfunction(f::T, data::RheoTimeData) where T<:Function
    stressfunction(data,f)
end

"""
    stressfunction!(data::RheoTimeData, f::T) where T<:Function
    stressfunction!(f::T, data::RheoTimeData) where T<:Function

In-place version of `stressfunction`.
"""
function stressfunction!(data::RheoTimeData, f::T) where T<:Function
    logadd_process!(data, :stressfunction, params=(f=f,), comment="Stress function applied to timeline" )
    _mapdata!(f,data.σ,data.t) 
    return data
end

function stressfunction!(f::T, data::RheoTimeData) where T<:Function
    stressfunction!(data,f)
end


#=
------------------------------------------------------------------------------
These functions provide convenient tools to design steps, ramps and
other input patterns for the stress and strain data

Functions may also return an anonymous function if the parameter t is omitted.
------------------------------------------------------------------------------
=#
"""
    hstep(t; offset=0., amp=1.)

Step generation function for use with `stressfunction` or `strainfunction`.
`offset` keyword arguent determines start of step. `amp` argument determines
amplitude (height) of step.
"""
function hstep(t; offset=0., amp=1.)
    return (t<offset) ? 0 : amp
end

function hstep(; kwargs...)
    return t -> hstep(t; kwargs...)
end

"""
    ramp(t; offset=0., gradient=1.)

Ramp signal generation function for use with `stressfunction` or `strainfunction`.
`offset` keyword argument determines start of ramp. `gradient` argument determines
the linear gradient of the ramp.
"""
function ramp(t; offset=0., gradient=1.)
    return (t<offset) ? 0 : (t-offset) * gradient
end

function ramp(; kwargs...)
    return t -> ramp(t; kwargs...)
end

"""
    stairs(t; offset=0., amp=1., width=1.)

Stairs signal generation function for use with `stressfunction` or `strainfunction`.
Equivalent to additional steps being added every `width` seconds. `offset` keyword
argument determines start of stairs signal. `amp` argument determines the height
of each additional step.
"""
function stairs(t; offset=0., amp=1., width=1.)
    return (t<offset) ? 0 : amp * floor((t - offset)/width + 1)
end

function stairs(; kwargs...)
    return t -> stairs(t; kwargs...)
end

"""
    square(t; offset=0., amp=1., period=1., width=0.5*period)

Square signal generation function for use with `stressfunction` or `strainfunction`.
`offset` keyword argument determines start of square signal. `amp` argument determines
the height of each square pulse. `period` determines the period of one off/on section
of the square wave signal. `width` determines the width of each square pulse.
"""
function square(t; offset=0., amp=1., period=1., width=0.5*period)
    return t<offset ? 0. : ( ((t-offset)%period) < width ? amp : 0.)
end

function square(; kwargs...)
    return t -> square(t; kwargs...)
end

"""
    sawtooth(t; offset=0., amp=1., period=1.)

Sawtooth signal generation function for use with `stressfunction` or `strainfunction`.
`offset` keyword argument determines start of sawtooth signal. `amp` argument determines
the height of each sawtooth pulse. `period` determines the period of the sawtooth wave signal.
"""
function sawtooth(t; offset=0., amp=1., period=1.)
    return t<offset ? 0. : amp*((t-offset)%period)/period
end

function sawtooth(; kwargs...)
    return t -> sawtooth(t; kwargs...)
end

"""
    triangle(t; offset=0., amp=1., period=1.)

Triangle signal generation function for use with `stressfunction` or `strainfunction`.
`offset` keyword argument determines start of triangle signal. `amp` argument determines
the height of each triangle pulse. `period` determines the period of the triangle wave signal.
`width` determines the width of the triangles.
"""
function triangle(t; offset=0., amp=1., period=1., width=0.5*period)
    if t<offset
        return 0.
    end
    t0=(t-offset)%period
    if t0 <width
        return amp*t0/width
    else
        return amp*(period-t0)/(period-width)
    end
end

function triangle(; kwargs...)
    return t -> triangle(t; kwargs...)
end








"""
    modulusfunction(data::RheoFreqData, Gp::T1, Gpp::T2) where {T1<:Function, T2<:Function}
    modulusfunction(Gp::T1, Gpp::T2, data::RheoFreqData) where {T1<:Function, T2<:Function}

Returns a new `RheoFreqData` with a modulus calculated using the function 'f' provided, 
applied to freqency values present in the input 'data'.

Normally used with a `RheoFreqData` generated using the `frequencyspec` function.
"""
function modulusfunction(data::RheoFreqData, Gp::T1, Gpp::T2) where {T1<:Function, T2<:Function}
    log = logadd_process(data, :modulusfunction, params=(Gp=Gp,Gpp=Gpp), comment="Modulus function applied to frequency range" )
    return RheoFreqData(convert(Vector{RheoFloat}, map(Gp, data.ω)), convert(Vector{RheoFloat}, map(Gpp, data.ω)),  data.ω, log)
end

function modulusfunction(Gp::T1, Gpp::T2, data::RheoFreqData) where {T1<:Function, T2<:Function}
    modulusfunction(data,Gp,Gpp)
end

"""
    modulusfunction!(data::RheoFreqData, f::T) where T<:Function
    modulusfunction!(f::T, data::RheoFreqData) where T<:Function

In-place version of `modulusfunction`.
"""
function modulusfunction!(data::RheoFreqData, Gp::T1, Gpp::T2) where {T1<:Function, T2<:Function}
    logadd_process!(data, :modulusfunction, params=(Gp=Gp,Gpp=Gpp), comment="Modulus function applied to frequency range" )
    _mapdata!(Gp,data.Gp,data.ω) 
    _mapdata!(Gpp,data.Gpp,data.ω) 
    return data
end

function modulusfunction!(Gp::T1, Gpp::T2, data::RheoFreqData) where {T1<:Function, T2<:Function}
    modulusfunction!(data,Gp,Gpp)
end


















#=
------------------------------------------------------------------------------------
Functionality not yet integrated with new data generation format, see GitHub issues:
-----------------------------------------------------------------------------------
"""
    stepgen(t_total::Real, t_on::Real; t_trans::Real = 0.0, stepsize::Real = 1.0)

Generate RheologyData struct with a step loading of height 1.0. If `t_trans` is 0.0 then
the step is instantaneous, otherwise the step is approximated by a logistic function
approximately centered at `t_on`.

# Arguments

- `t_total`: Total time length of data
- `t_on`: Step on time
- `t_trans`: Step transition time
- `stepsize`: Time sampling period
"""
function stepgen(t_total::Real,
                  t_on::Real;
                  t_trans::Real = 0.0,
                  stepsize::Real = 1.0)

    t = collect(0.0:stepsize:t_total)

    data = similar(t)
    # for smooth transition
    if t_trans>0.0
        k = 10.0/t_trans
        data = 1.0./(1 .+ exp.(-k*(t.-t_on)))
    # discrete jump
    elseif t_trans==0.0
        for (i, t_i) in enumerate(t)
            if t_i>=t_on
                data[i] = 1.0
            elseif t_i<t_on
                data[i] = 0.0
            end
        end
    end

    RheologyData(data, t, ["stepgen: t_total: $t_total, t_on: $t_on, t_trans: $t_trans, stepsize: $stepsize"])

end

"""
    noisegen(t_total::Real; seed::Union{Int, Nothing} = nothing, stepsize::Real = 1.0)

Generate uniform random noise of maximum amplitude +/- 1.0. If reproducibility is required,
always use the same number in the `seed` keyword argument with the same non-negative integer.

# Arguments

- `t_total`: Total time length of data
- `seed`: Seed used for random number generation
- `baseval`: Initial amplitude before oscillation started
- `stepsize`: Time sampling period
"""
function noisegen(t_total::Real; seed::Union{Int, Nothing} = nothing, stepsize::Real = 1.0)

    # conditional loading of random as this is the only RHEOS function which requires it
    @eval import Random: seed!, rand

    t = collect(0.0:stepsize:t_total)

    if seed!=nothing
        # use specified seed
        @assert seed>=0 "Seed integer must be non-negative"
        log = ["noisegen: seed: $seed, stepsize: $stepsize"]
        seed!(seed)
    else
        log = ["noisegen: seed: nothing, stepsize: $stepsize"]
    end

    data = (2*rand(eltype(t), length(t)) .- 1)



    RheologyData(data, t, log)

end

"""
    repeatdata(self::RheologyData, n::Integer)

Repeat a given RheologyData data set `n` times.
"""
function repeatdata(self::RheologyData, n::Integer; t_trans = 0.0)

    @assert self.σ==self.ϵ "Repeat data only works when σ==ϵ which is the state of RheologyData after being generated using the built-in RHEOS data generation functions"

    dataraw = self.σ

    @assert constantcheck(self.t) "Data sample-rate must be constant"

    step_size = self.t[2] - self.t[1]

    t = collect(0.0:step_size:(self.t[end]*n))

    # smooth transition between repeats
    if t_trans>0.0
        elbuffer = round(Int, (1/2)*(t_trans/step_size))
        selflength = length(self.t)

        # ensure smooth transition from end of data to new beginning so no discontinuities
        data_smooth_end = stepgen(self.t[end], self.t[end]; t_trans = t_trans, amplitude = dataraw[1] - dataraw[end], baseval = dataraw[selflength - elbuffer], stepsize = step_size)
        data_smooth_start = stepgen(self.t[end], 0.0; t_trans = t_trans, amplitude = dataraw[1] - dataraw[end], baseval = dataraw[selflength - elbuffer], stepsize = step_size)
        data_smoother = data_smooth_end + data_smooth_start

        data_single = self + data_smoother

        data = repeat(data_single.data[1:end] - dataraw[selflength - elbuffer], outer=[n])

        for i = 1:(n-1)
            deleteat!(data, i*length(self.t) + (2-i))
        end

        # fix first repeat
        for (i, v) in enumerate(dataraw)
            if i<(elbuffer*100) && i<round(Int, length(self.t)/2)
                data[i] = v
            end
        end

        log = vcat(self.log, ["repeated data $n times with transition time $t_trans"])

        return RheologyData(data, t, log)

    # discrete jump
    elseif t_trans==0.0

        data = repeat(dataraw, outer=[n])

        for i = 1:(n-1)
            deleteat!(data, i*length(self.t) + (2-i))
        end

        log = vcat(self.log, ["repeated data $n times with transition time $t_trans"])

        return RheologyData(data, t, log)

    end

end
=#
