#!/usr/bin/env julia

"""
    function stepdata_generate(t_total::Float64, t_on::Float64, t_off::Float64, t_transition::Float64, amplitude::Float64, test_type::String; step_size::Union{String, Float64} = "auto")

Generate partially filled RheologyData struct with a step function approximation
(2 logisitc functions).
"""
function stepdata_generate(t_total::Float64, t_on::Float64, t_off::Float64, t_transition::Float64, amplitude::Float64, test_type::String; step_size::Union{String, Float64} = "auto" )

    if step_size == "auto"
        step_size = t_total/5000.0
    end

    t = collect(0:step_size:t_total)

	k = 10.0/t_transition

	data_controlled = amplitude./(1 + exp.(-k*(t-t_on))) - amplitude./(1 + exp.(-k*(t-t_off)))

    RheologyData(data_controlled, t, test_type, "INTERNAL STEP")

end
