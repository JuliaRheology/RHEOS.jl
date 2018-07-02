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

using PyPlot

"""
    repeatdata(self::RheologyArtificial, n::Integer)

Repeat a given RheologyArtificial generated data set n times.

Needs fix to increase accuracy over many repeats.
"""
function repeatdata(self::RheologyArtificial, n::Integer; t_trans = 0.1)

    t = collect(0.0:self.stepsize:(self.t[end]*n))

    elbuffer = round(Int, (1/2)*(t_trans/self.stepsize))
    selflength = length(self.t)

    # ensure smooth transition from end of self.data to new beginning so no discontinuities
    data_smooth_end = stepdata(self.t[end], self.t[end]; t_trans = t_trans, amplitude = self.data[1] - self.data[end], baseval = self.data[selflength - elbuffer], stepsize = self.stepsize)
    data_smooth_start = stepdata(self.t[end], 0.0; t_trans = t_trans, amplitude = self.data[1] - self.data[end], baseval = self.data[selflength - elbuffer], stepsize = self.stepsize)
    data_smoother = data_smooth_end + data_smooth_start

    data_single = self + data_smoother

    data = repeat(data_single.data[1:end] - self.data[selflength - elbuffer], outer=[n])

    for i = 1:(n-1)
        deleteat!(data, i*length(self.t) + (2-i))
    end

    # fix first repeat
    for (i, v) in enumerate(self.data)
        if i<(elbuffer*100) && i<round(Int, length(self.t)/2) 
            data[i] = v
        end
    end

    log = vcat(self.log, ["repeated data $n times"])

    RheologyArtificial(data, t, self.stepsize, log)

end

