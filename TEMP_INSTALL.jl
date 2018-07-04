#!/usr/bin/env julia

Pkg.update()

# From Julia Repository
Pkg.add("QuadGK")
Pkg.add("uCSV")
Pkg.add("ImageFiltering")
Pkg.add("Interpolations")
Pkg.add("NLopt")
Pkg.add("JLD")
Pkg.add("InverseLaplace")
Pkg.add("CSV")

# From rheos-cambridge forked repositories (so that versions are fixed and controlled)
try
    Pkg.clone("https://github.com/rheos-cambridge/MittagLeffler.jl")
catch
    # if already installed 
    LoadError
end

try
    Pkg.clone("https://github.com/rheos-cambridge/FastConv.jl")
catch
    # if already installed
    LoadError
end

#= temporarily add DataFrames, it is now included in the latest version of uCSV,
 but this is not yet on the official Julia repository =#
Pkg.add("DataFrames")

