using RHEOS
using Test

const tol = (eps(RHEOS.RheoFloat))^(0.125)

println("===============================================")
println("Testing RHEOS")
println("===============================================")
println("|")



include("symbols.jl")
include("definitions.jl")
include("base.jl")
include("datagen.jl")
include("IO.jl")
include("processing.jl")
include("interface.jl")

include(joinpath(@__DIR__, "models", "elements.jl"))
include(joinpath(@__DIR__, "models", "maxwell.jl"))
include(joinpath(@__DIR__, "models", "kelvinvoigt.jl"))
include(joinpath(@__DIR__, "models", "zener.jl"))
include(joinpath(@__DIR__, "models", "poynting-thomson.jl"))
