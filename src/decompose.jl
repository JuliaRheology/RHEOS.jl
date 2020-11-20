macro Name(arg)
    string(arg)
 end

"""
    decompose(model, data, output, r_count = 0)

Decomposes `model` into its components recursively, and stores the stress-strain data along the way.

# Arguments
- `model`: RheoModel to decompose
- `data`: RheoTimeData that is complete ie has no missing columns
- `output`: This is the returned NamedTuple. Built recursively.
- `r_count`: The recursion count. Helps keep track of how deep the recursion is, useful for printing a branched structure.
"""

function decompose(model::RheoModel, data::RheoTimeData, output::Dict=Dict())
    
    output["data"] = data

    # check if we have reached a spring or dashpot
    if model.description.type == "basic"
        # we can be sure that there will be no further division of this model
        # output["data"] = data  
        print("basic", typeof(output))
        return output
    end

    # The model can be split into simpler elements. 
    # Check if it is a series or parallel connection
    if model.description.type == "series"
        datae = extract(data, stress_only);
    elseif model.description.type == "parallel"
        datae = extract(data, strain_only);
    end
    print(model.description.type)

    # model_params = model.params # to store a separate copy of the model parameters - to use for param matching
    model_param_values = [mp for mp in model.params] # store sequence of parameter values separately
    print("type ", typeof(model_param_values))

    # start the recursion 
    for (index, c) in enumerate(model.description.components)
        # c will be like :Dashpot or :Maxwell or :Spring
        component = eval(c)
        print(component.name)

        """ 
        # --- Matching the actual parameter symbol ---

        # Get that component's parameter values
        pkeys = component.params # this is the list of Symbols
        print(pkeys)
        pval = Vector{Float64}()
        for p in pkeys # here we are iterating over the actual symbol
            print(model_params[p])
            push!(pval, model_params[p]) # want to push model parameter values to components
        end
        ptuple = (; zip(pkeys, pval)...)
        """

        # --- Matching the parameters in order ---
        pkeys = component.params
        pval = Vector{Float64}()
        for i = 1:length(pkeys)
            push!(pval, model_param_values[1])
            popfirst!(model_param_values)
        end
        ptuple = (; zip(pkeys, pval)...)
        

        # Create the component's RheoModel, and complete the dataset
        component_RM = RheoModel(component, ptuple);
        # NOTE: datae has one field missing, either stress or strain, depending on if it is a series or parallel branch
        component_data = modelpredict(datae, component_RM); 

        # -- Recurse with that
        # output["data"] = data
        # output[c.name] = Dict("data"=>c_data) 
        output[component.name] = Dict()
        # print(keys(output))

        output[component.name] = decompose(component_RM, component_data)

        # println("\t"^r_count, "> ",c.name, " = ", component_name, ", ", ptuple)

    end
    return output
end