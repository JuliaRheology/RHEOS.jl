#!/usr/bin/env julia
push!(LOAD_PATH, "../deps")

module RHEOS
# local dependencies, stored in "/RHEOS/deps/"
using MittagLeffler
using FastConv
# install using Pkg.add from within Julia REPL
using uCSV
using ImageFiltering
using Interpolations
using NLopt
using PyPlot
# using Plots; gr() # add support for different plotting backends

######################################################
# utility.jl
export deriv, trapz

######################################################
# preproc.jl
export var_resample, downsample, fixed_resample

######################################################
# proc.jl
export leastsquares_init, objectivefunc, boltzconvolve, boltzintegral

######################################################
# models.jl
export mittleff, moduli

export G_SLS, G_SLS2, G_burgers,
        G_springpot, G_fractKV,
        G_fractmaxwell, G_fractzener,
        G_fractspecial

export J_SLS, J_SLS2, J_burgers,
        J_springpot, J_fractKV,
        J_fractmaxwell, J_fractzener

######################################################
# hilevel_init.jl
export RheologyData, fileload, RheologyModel

######################################################
# hilevel_preproc.jl
export var_resample, downsample, fixed_resample, smooth, mapbackdata

######################################################
# hilevel_proc.jl
export modelfit!, modelcomplete!

######################################################
# hilevel_postproc.jl
export fiteval

######################################################
# code
include("utility.jl")
include("hilevel_init.jl")
include("preproc.jl")
include("proc.jl")
include("models.jl")
include("hilevel_preproc.jl")
include("hilevel_proc.jl")
include("hilevel_postproc.jl")
######################################################
end
