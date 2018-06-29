#!/usr/bin/env julia
using PyPlot
import Base: +, -

###############################
#~ Definitions and Overloads ~#
###############################

struct RheologyArtificial

    # original data
    data::Array{Float64,1}
    t::Array{Float64,1}

    # constant sample rate step size
    stepsize::Float64

    # operations applied, stores history of which functions (including arguments)
    log::Array{String,1}

end

function +(self1::RheologyArtificial, self2::RheologyArtificial)

    @assert self1.stepsize==self2.stepsize "Step size must be same for both datasets"

    # get time array
    if length(self1.t) >= length(self2.t)
        t  = self1.t
    else
        t = self2.t
    end

    # init data array and fill by summing over each argument's indices
    data = zeros(length(t))

    for i in 1:length(self1.t)
        data[i] += self1.data[i]
    end

    for i in 1:length(self2.t)
        data[i] += self2.data[i]
    end

    # log
    log = vcat(self1.log, self2.log, ["previous two logs added"])

    RheologyArtificial(data, t, self1.stepsize, log)

end

#####################
#~ Data Generators ~#
#####################

"""
    function stepdata(t_total::Float64, t_on::Float64; t_trans::Float64 = (t_total - t_on)/100.0, amplitude::Float64 = 1.0, start_val::Float64 = 0.0, step_size::Float64 = 0.5)

Generate RheologyData struct with a step function approximated by a logisitic function.
Sends generated data to both ϵ and σ so as to be compatible with all types of tests.
"""
function stepdata(t_total::Float64, 
                  t_on::Float64;
                  t_trans::Float64 = (t_total - t_on)/100.0,
                  amplitude::Float64 = 1.0,
                  start_val::Float64 = 0.0,
                  step_size::Float64 = 0.5, )

    t = collect(0:step_size:t_total)

	k = 10.0/t_trans

	data = start_val + amplitude./(1 + exp.(-k*(t-t_on)))

    RheologyArtificial(data, t, step_size, ["stepdata: t_total: $t_total, t_on: $t_on, t_trans: $t_trans, amplitude: $amplitude, step_size: $step_size"])

end

foo = stepdata(1000.0, 100.0; start_val = 0.0, step_size = 0.5)
bar = stepdata(1200.0, 500.0; amplitude = -1.0)


baz = foo + bar;

plot(foo.t, foo.data)
plot(bar.t, bar.data)
plot(baz.t, baz.data)
show()