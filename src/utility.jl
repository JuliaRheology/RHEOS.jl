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
    yDot = zeros(Float64, length(y))
    for i in 1:length(y)
        if i==1
            # right handed difference for lower boundary
            @inbounds yDot[i] = (y[i+1] - y[i])/(x[i+1] - x[i])
        elseif i==length(y)
            # left handed difference for upper boundary
            @inbounds yDot[i] = (y[i] - y[i-1])/(x[i] - x[i-1])
        else
            # central difference with uneven spacing
            @inbounds Δx₁ = x[i] - x[i-1]
            @inbounds Δx₂ = x[i+1] - x[i]
            @inbounds yDot[i] = (y[i+1]*Δx₁^2 + (Δx₂^2 - Δx₁^2)*y[i] - y[i-1]*Δx₂^2)/(Δx₁*Δx₂*(Δx₁ + Δx₂))
        end
    end
    yDot
end

function deriv(y::Array{Float64,1})::Array{Float64,1}

    # initialise zero array of length y
    yDot = zeros(Float64, length(y))
    for i in 1:length(y)
        if i==1
            # right handed difference for lower boundary
            @inbounds yDot[i] = (y[i+1] - y[i])
        elseif i==length(y)
            # left handed difference for upper boundary
            @inbounds yDot[i] = (y[i] - y[i-1])
        else
            # central difference with even unit 1 spacing
            @inbounds yDot[i] = (y[i+1] - y[i-1])/2.0
        end
    end
    yDot
end
