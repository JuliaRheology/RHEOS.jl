VERSION >= v"0.4.0-dev+6521" && __precompile__()

module InverseLaplace

using Compat

export ILt, setNterms

export  Weeks, WeeksErr, optimize, opteval, setparameters

export Talbot, GWR, ILT

export ILtPair, abserr, iltpair_power

export ilt, talbot, gwr

@compat abstract type AbstractILt end

"""
   ft::Talbot = Talbot(func::Function, Nterms::Integer=32)

return `ft`, which estimates the inverse Laplace transform of `func` with
the fixed Talbot algorithm. `ft(t)` evaluates the transform at `t`.  You may
want to tune `Nterms` together with `setprecision(BigFloat,x)`.

# Example

Compute the inverse transform of the transform of `cos` at argument `pi/2`.
```julia-repl
julia> ft = Talbot(s -> s/(s^2+1), 80);

julia> ft(pi/2)
0.0
```

Note that given `Float64` input, the precision of the returned value may not be satisfactory.
```julia-repl
julia> ft(pi/2)
-3.5366510684573195e-5

julia> Float64(ft(big(pi)/2))
2.114425886215604e-49
```
"""
type Talbot{T<:Base.Callable} <: AbstractILt
    LSfunc::T  # Laplace space function
    Nterms::Int
end

Talbot(LSfunc::Base.Callable) = Talbot(LSfunc,32)

"""
   ft::GWR = GWR(func::Function, Nterms::Integer=16)

return `ft`, which estimates the inverse Laplace transform of `func` with
the GWR algorithm. `ft(t)` evaluates the transform at `t`.

# Example

Compute the inverse transform of the transform of `cos` at argument `pi/2`.
```
julia> ft = GWR(s -> s/(s^2+1), 16);

julia> ft(pi/2)
-0.001
```
"""
type GWR{T<:Base.Callable} <: AbstractILt
    LSfunc::T  # Laplace space function
    Nterms::Int
end

GWR(LSfunc::Base.Callable) = GWR(LSfunc,16)

"""
    ILT(function, Nterms=32)

This is an alias for the default `Talbot()` method.
"""
ILT(args...) = Talbot(args...)

function Base.show(io::IO, ailt::AbstractILt)
    print(io, string(typeof(ailt)), "(Nterms=", ailt.Nterms,')')
end

# TODO get rid of this in favor of above.
# Rely on broadcasting, as well.

type ILt{T<:Base.Callable, V<:Base.Callable} <: AbstractILt
    LSfunc::T
    iltfunc::V
    Nterms::Int
end

"""
    itrans = ILt(func, iltfunc=talbot, Nterms=32)

 * deprecated. Use, ILT, Talbot, GWR, or Weeks instead *

return an object that estimates the inverse Laplace transform of
the function `func` using the algorithm implemented by function `iltfunc`.
`itrans(t)` estimates the inverse transform for argument `t`.  The
accuracy of the estimates depends strongly on the choice of `iltfunc`, `t`, `Nterms`,
and the precision of the data type of the argument to `func`. The default value
of `32` may give extremely inaccurate estimates.

`iltfunc` may be either `talbot` or `gwr`.

## Example

```julia
julia> itr = ILt( s -> 1/(1+s^2), talbot);

julia> itr([ pi/4, pi/2, 3*pi/4, -pi])
4-element Array{Float64,1}:
  0.707107
  1.0
  0.707107
 -3.66676e-12
```
"""
ILt(func,iltfunc) = ILt(func, iltfunc, 32)

# Make this the default
"""
    ilt(func::Function, t::AbstractFloat, M::Integer=32)

`ilt` is an alias for the default inverse Laplace transform method `talbot`.
"""
ilt(args...) = talbot(args...)
ILt(func) = ILt(func, talbot, 32)


"""
    setNterms(ailt::AbstractILt, Nterms::Integer)

set the number of terms used in the inverse Laplace tranform `ailt`. If
`ailt` stores internal data, it will be recomputed, so that subsequent
calls `ailt(t)` reflect the new value of `Nterms`.
"""
setNterms(ailt::AbstractILt, N::Integer) = (ailt.Nterms = N)

(ailt::ILt)(t) = ailt.iltfunc(ailt.LSfunc, t, ailt.Nterms)
(ailt::ILt)(t,N) = ailt.iltfunc(ailt.LSfunc, t, N)

(ailt::Talbot)(t) = talbot(ailt.LSfunc, t, ailt.Nterms)
(ailt::Talbot)(t, n::Integer) = talbot(ailt.LSfunc, t, n)

(ailt::GWR)(t) = gwr(ailt.LSfunc, t, ailt.Nterms)
(ailt::GWR)(t, n::Integer) = gwr(ailt.LSfunc, t, n)

include("fixed_talbot.jl")
include("gwr.jl")
include("weeks.jl")
include("test.jl")

end # module
