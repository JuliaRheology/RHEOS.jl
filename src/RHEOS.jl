#!/usr/bin/env julia
__precompile__()

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
using JLD
# experimental dependency not automatically imported from base.jl
using Base.Threads
# using Plots; gr() # add support for different plotting backends

println("\n===========================")
println("Number of threads in use: ", nthreads())
println("===========================\n")

######################################################
# base_utility.jl
export deriv, trapz, mittleff, RheologyModel

######################################################
# base_preproc.jl
export var_resample, downsample, fixed_resample

######################################################
# base_proc.jl
export leastsquares_init, objectivefunc, boltzconvolve, boltzintegral

######################################################
# base_models.jl
export moduli

export G_SLS, G_SLS2, G_burgers,
        G_springpot, G_fractKV,
        G_fractmaxwell, G_fractzener,
        G_fractspecial

export J_SLS, J_SLS2, J_burgers,
        J_springpot, J_fractKV,
        J_fractmaxwell, J_fractzener

######################################################
# rheology_utility.jl
export RheologyData, fileload

######################################################
# rheology_preproc.jl
export var_resample, downsample, fixed_resample, smooth, mapbackdata

######################################################
# rheology_proc.jl
export modelfit!, modelcomplete!

######################################################
# rheology_postproc.jl
export fiteval, saveresult, loadresult

######################################################
# afm_contactmodels.jl
export contact_hertz, contact_threshold, contact_none

######################################################
# afm_utility.jl
export AFMfileload

######################################################
# Main functionality
include("base_utility.jl")
include("base_preproc.jl")
include("base_proc.jl")
include("base_models.jl")

# High level rheology interface
include("rheology_utility.jl")
include("rheology_preproc.jl")
include("rheology_proc.jl")
include("rheology_postproc.jl")

# High level AFM/contact mechanics
include("afm_contactmodels.jl")
include("afm_utility.jl")

######################################################
end
