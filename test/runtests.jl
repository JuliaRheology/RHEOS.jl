using RHEOS
@static if VERSION < v"0.7.0-DEV.2005"
    using Base.Test
else
    using Test
end

# run base tests
include("base.jl")

