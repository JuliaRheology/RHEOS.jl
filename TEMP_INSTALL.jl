#!/usr/bin/env julia

Pkg.add("QuadGK")
Pkg.add("uCSV")
Pkg.add("ImageFiltering")
Pkg.add("Interpolations")
Pkg.add("NLopt")
Pkg.add("JLD")
Pkg.add("InverseLaplace")

#= temporarily add DataFrames, it is now included in the latest version of uCSV,
 but this is not yet on  the official Julia repository =#
Pkg.add("DataFrames")

Pkg.update()