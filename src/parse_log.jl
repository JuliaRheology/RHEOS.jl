# parse_log.jl

"""
    parse_log(log)

Parse a log of actions and extract model fitting entries into a dictionary.

# Arguments
- `log`: A list of tuples representing log entries.

# Returns
A dictionary where keys are model names and values are lists of dictionaries,
each containing model parameters, error, and index from the log.

# Example
```julia
data.log = [...]  # Define your log data
models_data = parse_log(data.log)
println(models_data)

"""
function parse_log(log)
    models_dict = Dict{String, Vector{Dict{String, Any}}}()

    for (idx, entry) in enumerate(log)
        if entry.action.funct == :modelfit
            model_name = entry.info.model_name
            model_params = entry.info.model_params
            error = entry.info.error

            model_entry = Dict("params" => model_params, "error" => error, "index" => idx)

            if haskey(models_dict, model_name)
                push!(models_dict[model_name], model_entry)
            else
                models_dict[model_name] = [model_entry]
            end
        end
    end

    return models_dict
end
