#!/usr/bin/env julia
__precompile__()

module RHEOS
# installed using Pkg.clone from rheos-cambridge forked repos
using MittagLeffler
using FastConv
# install using Pkg.add from Julia central package repository
using InverseLaplace
using uCSV
using ImageFiltering
using Interpolations
using NLopt
using JLD
using DataFrames
# experimental dependency not automatically imported from base.jl
using Base.Threads
# for overloading in datagen.jl
import Base: +, -, *

# println("\n===========================")
# println("Number of threads in use: ", nthreads())
# println("===========================\n")

######################################################
# base.jl
export deriv, trapz, mittleff, closestindex

export var_resample, downsample, fixed_resample

export leastsquares_init, objectivefunc, boltzconvolve, boltzintegral

######################################################
# models.jl
export G_SLS, G_SLS2, G_burgers,
        G_springpot, G_fractKV,
        G_fractmaxwell, G_fractzener,
        G_fractspecial

export J_SLS, J_SLS2, J_burgers,
        J_springpot, J_fractKV,
        J_fractmaxwell, J_fractzener

export modeldatabase

######################################################
# definitions.jl
export RheologyArtificial, RheologyData, fileload, RheologyModel, RheologyModelTemp

######################################################
# datagen.jl
export stepdata, rampdata, oscillatordata, repeatdata

######################################################
# processing.jl
export var_resample, downsample, fixed_resample, smooth, mapbackdata

export modelfit, modelpredict

export savedata, loaddata, savemodel, loadmodel, exportdata

######################################################
# Main functionality
include("base.jl")
include("models.jl")

# High level rheology interface
include("definitions.jl")
include("datagen.jl")
include("processing.jl")
######################################################
end
