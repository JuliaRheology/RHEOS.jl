#!/usr/bin/env julia

"""
    trapz(y::Array{Float64,1}, x::Array{Float64,1}[; init::Float64 = 0.0])

Array based trapezoidal integration of y with respect to x.

Limits of integration defined by the start and end points of the arrays. 'init'
keyword argument is used for setting an initial condition.
"""
function trapz(y::Array{Float64,1}, x::Array{Float64,1}; init::Float64=0.0)::Float64

    @assert length(x)==length(y) "X and Y array length must match."

    # init*2 to simplify final division by 2
    r = 2.0*init

    # trapezoidal rule
    for i in 2:length(x)
        @inbounds r += (y[i-1] + y[i])*(x[i] - x[i-1])
    end

    # return summation
    r/2.0
end


"""
    deriv(y::Array{Float64,1}[, x::Array{Float64,1}])

Given two arrays of data, x and y, calculate dy/dx using central difference
method and forward and backward difference for array boundaries.

If only one array passed in, y, then calculates derivative y with respect to
array element indices.

# Example

```jldoctest
x = collect(0.0:0.1:5.0)

y = x.^2

dy_dx = deriv(y, x) # ~= 2*x
```
"""
function deriv(y::Array{Float64,1}, x::Array{Float64,1})::Array{Float64,1}

    # assert y and x arrays are same length
    @assert length(y)==length(x) "X and Y Array lengths must match."

    # initialise zero array of length y
    ydot = zeros(Float64, length(y))
    for i in 1:length(y)
        if i==1
            # right handed difference for lower boundary
            @inbounds ydot[i] = (y[i+1] - y[i])/(x[i+1] - x[i])
        elseif i==length(y)
            # left handed difference for upper boundary
            @inbounds ydot[i] = (y[i] - y[i-1])/(x[i] - x[i-1])
        else
            # central difference with uneven spacing
            @inbounds Δx₁ = x[i] - x[i-1]
            @inbounds Δx₂ = x[i+1] - x[i]
            @inbounds ydot[i] = (y[i+1]*Δx₁^2 + (Δx₂^2 - Δx₁^2)*y[i] - y[i-1]*Δx₂^2)/(Δx₁*Δx₂*(Δx₁ + Δx₂))
        end
    end
    ydot
end

function deriv(y::Array{Float64,1})::Array{Float64,1}

    # initialise zero array of length y
    ydot = zeros(Float64, length(y))
    for i in 1:length(y)
        if i==1
            # right handed difference for lower boundary
            @inbounds ydot[i] = (y[i+1] - y[i])
        elseif i==length(y)
            # left handed difference for upper boundary
            @inbounds ydot[i] = (y[i] - y[i-1])
        else
            # central difference with even unit 1 spacing
            @inbounds ydot[i] = (y[i+1] - y[i-1])/2.0
        end
    end
    ydot
end

"""
    mittleff(α::Float64[, β::Float64], x::Array{Float64,1})

Threaded convience wrapper around MittagLeffler.mittleff for array arguments. X
array is limited to float type in accordance with the requirements of RHEOS.
"""
function mittleff(α::Float64, x_array::Array{Float64,1})::Array{Float64,1}

    # initialise array
    y = Array{Float64}(length(x_array))
    # do static scheduled threaded loop with no bounds checking
    @threads for i = 1:length(x_array)
        @inbounds y[i] = MittagLeffler.mittleff(α, x_array[i])
    end
    # return
    y
end

function mittleff(α::Float64, β::Float64, x_array::Array{Float64,1})::Array{Float64,1}

    # initialise array
    y = Array{Float64}(length(x_array))
    # do static scheduled threaded loop with no bounds checking
    @threads for i = 1:length(x_array)
        @inbounds y[i] = MittagLeffler.mittleff(α, β, x_array[i])
    end
    # return
    y
end
# NON THREADED VERSIONS
# function mittleff(α::Float64, xList::Array{Float64,1})::Array{Float64,1}
#     # call mittagleffler within an array comprehension
#     y = [MittagLeffler.mittleff(α, x) for x in xList]
# end
#
# function mittleff(α::Float64, β::Float64, xList::Array{Float64,1})::Array{Float64,1}
#     # call mittagleffler within an array comprehension
#     y = [MittagLeffler.mittleff(α, β, x) for x in xList]
# end


"""
    RheologyModel(form::Function, singularity::Bool, test_type::String)

Struct which contains the functional form of a chosen model, a bool stating whether
or not that functional form contains a singularity at t=0, and whether it is
considered a "strlx" (strain controlled) or "creep" (stress controlled) modulus.
"""
struct RheologyModel

    form::Function

    singularity::Bool

    test_type::String

end
