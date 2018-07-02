#!/usr/bin/env julia

"""
    stepdata(t_total::Float64, t_on::Float64; t_trans::Float64 = (t_total - t_on)/100.0, amplitude::Float64 = 1.0, baseval::Float64 = 0.0, stepsize::Float64 = 0.5)

Generate RheologyArtificial struct with a step function approximated by a logisitic function.
"""
function stepdata(t_total::Float64, 
                  t_on::Float64;
                  t_trans::Float64 = (t_total - t_on)/100.0,
                  amplitude::Float64 = 1.0,
                  baseval::Float64 = 0.0,
                  stepsize::Float64 = 0.5, )

    t = collect(0.0:stepsize:t_total)

	k = 10.0/t_trans

	data = baseval + amplitude./(1 + exp.(-k*(t-t_on)))

    RheologyArtificial(data, t, stepsize, ["stepdata: t_total: $t_total, t_on: $t_on, t_trans: $t_trans, amplitude: $amplitude, baseval: $baseval, stepsize: $stepsize"])

end

"""
    rampdata(t_total::Float64, t_start::Float64, t_stop::Float64; amplitude::Float64 = 1.0, baseval::Float64 = 0.0, stepsize::Float64 = 0.5)

Generate RheologyArtificial struct with a ramp function.
"""
function rampdata(t_total::Float64,
                  t_start::Float64,
                  t_stop::Float64;
                  amplitude::Float64 = 1.0,
                  baseval::Float64 = 0.0,
                  stepsize::Float64 = 0.5 )

    t = collect(0.0:stepsize:t_total)

    m = amplitude/(t_stop - t_start)

    data = baseval + zeros(length(t))

    for (i, v) in enumerate(t)

        if t_start<=v<=t_stop
            data[i] += data[i-1] + stepsize*m
        elseif v>t_stop
            data[i] = data[i-1]
        end

    end

    RheologyArtificial(data, t, stepsize, ["rampdata: t_total: $t_total, t_start: $t_start, t_stop: $t_stop, amplitude: $amplitude, baseval: $baseval, stepsize: $stepsize"])

end

"""
    oscillatordata(t_total::Float64, frequency::Float64; t_start::Float64 = 0.0, phase::Float64 = 0.0, amplitude::Float64 = 1.0, baseval::Float64 = 0.0, stepsize::Float64 = 0.5)

Generate RheologyArtificial struct with a sinusoidal function.
"""
function oscillatordata(t_total::Float64,
                        frequency::Float64;
                        t_start::Float64 = 0.0,
                        phase::Float64 = 0.0,
                        amplitude::Float64 = 1.0,
                        baseval::Float64 = 0.0,
                        stepsize::Float64 = 0.5 )

    t = collect(0.0:stepsize:t_total)

    data = baseval + zeros(length(t))

    for (i, v) in enumerate(t)

        if v>=t_start
            data[i] += amplitude*sin(2*π*frequency*(v - t_start) + phase)
        end

    end

    RheologyArtificial(data, t, stepsize, ["oscillatordata: t_total: $t_total, frequency: $frequency, t_start: $t_start, phase: $phase, amplitude: $amplitude, baseval: $baseval, stepsize: $stepsize"])

end

"""
    repeatdata(self::RheologyArtificial, n::Integer)

Repeat a given RheologyArtificial generated data set n times.

Needs fix to increase accuracy over many repeats.
"""
function repeatdata(self::RheologyArtificial, n::Integer; t_trans = 0.5)

    t = collect(0.0:self.stepsize:(self.t[end]*n))
    t = t[1:(length(t)-1)]

    # ensure smooth transition from end of self.data to new beginning so no discontinuities
    data_smooth_single = stepdata(self.t[end], t_trans; t_trans = t_trans, amplitude = self.data[1] - self.data[end], baseval = self.data[end], stepsize = self.stepsize)

    data_single = data_smooth_single + self

    data = repeat(data_single.data[1:(length(data_single.data)-1)] - self.data[1], outer=[n])

    # fix first repeat
    for (i, v) in enumerate(self.data)
        data[i] = v
    end

    log = vcat(self.log, ["repeated data $n times"])

    RheologyArtificial(data, t, self.stepsize, log)

end

