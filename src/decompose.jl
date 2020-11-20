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
    # print(model.description.type)

    # start the recursion 
    for c in model.description.components
        # c will be like :Dashpot or :Maxwell or :Spring
        component = eval(c)
        print(component.name)

        # Get that component's parameter values
        pkeys = component.params # this is the list of Symbols
        print(pkeys)
        pval = Vector{Float64}()
        for p in pkeys
            print(model.params[p])
            push!(pval, model.params[p]) # want to push model parameter values to components
        end
        ptuple = (; zip(pkeys, pval)...)
        

        # Create the component's RheoModel, and complete the dataset
        component_RM = RheoModel(component, ptuple);
        # NOTE: datae has one field missing, either stress or strain, depending on if it is a series or parallel branch
        component_data = modelpredict(datae, component_RM); 

        # -- STEP 4.4: Recurse with that
        # output["data"] = data
        # output[c.name] = Dict("data"=>c_data) 
        output[component.name] = Dict()
        # print(keys(output))

        output[component.name] = decompose(component_RM, component_data)

        # println("\t"^r_count, "> ",c.name, " = ", component_name, ", ", ptuple)

    end
    return output
end