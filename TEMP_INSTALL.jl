#!/usr/bin/env julia

Pkg.add("QuadGK")
#= temporarily add DataFrames, it is now included in the latest version of uCSV,
 but this is not yet on  the official Julia repository =#
Pkg.add("DataFrames")
##
Pkg.add("uCSV")
Pkg.add("ImageFiltering")
Pkg.add("Interpolations")
Pkg.add("NLopt")
Pkg.add("PyPlot")
Pkg.add("JLD")

Pkg.update()