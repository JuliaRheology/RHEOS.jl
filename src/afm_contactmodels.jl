#!/usr/bin/env julia

"""
    contact_hertz()
"""
function contact_hertz

end

"""
    contact_threshold()
"""
function contact_threshold(f::Array{Float64,1}, δ::Array{Float64,1}; param::Float64 = NaN64)

    # get index of contact
    cp_index = minimum(find(f .>= param))

    # return 
    cp_index

end

"""
    contact_none()
"""
function contact_none(f::Array{Float64,1}, δ::Array{Float64,1}; param::Float64 = NaN64)

    # return index 1
    1
    
end
