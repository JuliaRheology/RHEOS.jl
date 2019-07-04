#!/usr/bin/env julia

#=
----------------------------------------------------------------------------
timeline generation function forms the basis of any load-generation workflow
----------------------------------------------------------------------------
=#
"""
    timeline(;t_start::Real=0., t_end::Real=10., step::Real=(t_end-t_start)/250.)

Generate RheoTimeData struct with only the time data abd kig,

# Arguments

- `t_start`: Starting time, typically 0
- `t_end`: End time
- `step`: Time between sample
"""
function timeline(;t_start::Real=0., t_end::Real=10., step::Real=(t_end-t_start)/250.)

    log = OrderedDict{Any,Any}(:n=>1,"activity"=>"timeline created - t_start: $t_start, t_end: $t_end, step: $step")

    RheoTimeData([], [], collect(t_start:step:t_end), log)
end

#=
------------------------------------------------------------------------
strainfunction and stressfunction can be used to impose strain of stress
respectively. The subsequent convenience functions are passed to them to 
generate the appropriate loading.
------------------------------------------------------------------------
=#
function strainfunction(data::RheoTimeData, f::T) where T<:Function
    log = copy(data.log)
    log[:n] = log[:n] + 1
    log[string("activity_", log[:n])] = "Strain generated"
    #log = OrderedDict{Any,Any}("activity"=>"strain function", "data_source"=>d.log)
    return RheoTimeData(data.σ, convert(Vector{RheoFloat}, map(f, data.t)), data.t, log)
end

function stressfunction(d::RheoTimeData, f::T) where T<:Function
    log = copy(data.log)
    log[:n] = log[:n] + 1
    log[string("activity_", log[:n])] = "Stress generated"
    #log = OrderedDict{Any,Any}("activity"=>"stress function", "data_source"=>data.log)
    return RheoTimeData(convert(Vector{RheoFloat}, map(f, data.t)), data.ϵ, data.t, log)
end

#=
------------------------------------------------------------------------------
These functions provide convenient tools to design steps, ramps and
other input patterns for the stress and strain data

Functions may also return an anonymous function if the parameter t is omitted.
------------------------------------------------------------------------------
=#
function hstep(t; offset=0., amp=1.)
    return (t<offset) ? 0 : amp
end

function hstep(; offset=0., amp=1.)
    return t->(t<offset) ? 0 : amp
end

function ramp(t;offset=0., amp=1.)
    return (t<offset) ? 0 : (t-offset) * amp
end

function ramp(;offset=0., amp=1.)
    return t -> (t<offset) ? 0 : (t-offset) * amp
end

function stairs(t;offset=0., amp=1., width=1.)
    return (t<offset) ? 0 : amp * ceil((t-offset)/width)
end

function stairs(;offset=0., amp=1., width=1.)
    return t -> (t<offset) ? 0 : amp * ceil((t-offset)/width)
end

function square(t;offset=0., amp=1., period=1., width=0.5*period)
    return t<offset ? 0. : ( ((t-offset)%period) <width ? amp : 0.)
end

function square(;offset=0., amp=1., period=1., width=0.5*period)
    return t -> t<offset ? 0. : ( ((t-offset)%period) <width ? amp : 0.)
end

function sawtooth(t;offset=0., amp=1., period=1.)
    return t<offset ? 0. : amp*((t-offset)%period)/period
end

function sawtooth(;offset=0., amp=1., period=1.)
    return t -> t<offset ? 0. : amp*((t-offset)%period)/period
end

function triangle(t;offset=0., amp=1., period=1., width=0.5*period)
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

function triangle(;offset=0., amp=1., period=1., width=0.5*period)
    return t -> triangle(t,offset=offset, amp=amp, period=period, width=width)
end

function frequency_spec(;ω_start::Real=1.0e-2, ω_end::Real=1.0e2, step::Real=(ω_end-ω_start)/1.0e5)
    log = OrderedDict{Any,Any}(:n=>1,"activity"=>"frequency_spec created: ω_start: $ω_start, ω_end: $ω_end, step: $step")
    RheoFreqData([],[],collect(ω_start:step:ω_end),log)
end

#=
---------------------------------------------------------------------------
Functionality not yet integrated with new data generation format, see TODO:
---------------------------------------------------------------------------
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
