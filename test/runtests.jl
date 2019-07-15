using RHEOS
using Test

const tol = (eps(RHEOS.RheoFloat))^(0.125) 

# include("base.jl")
# include("datagen.jl")
include("IO.jl")
