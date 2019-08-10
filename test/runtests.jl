using RHEOS
using Test

const tol = (eps(RHEOS.RheoFloat))^(0.125)

include("definitions.jl")
include("base.jl")
include("datagen.jl")
include("IO.jl")
include("processing.jl")

include(joinpath(@__DIR__, "models", "elements.jl"))