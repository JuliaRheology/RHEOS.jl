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
export deriv, trapz, mittleff, closestindex, closestindices

export singularitytest, var_resample, downsample, fixed_resample

export leastsquares_init, objectivefunc, boltzconvolve, boltzintegral

######################################################
# definitions.jl
export RheologyArtificial, RheologyData, fileload, RheologyModel, RheologyModelTemp

######################################################
# models.jl
export SLS, SpringPot, FractionalMaxwell, FractionalKelvinVoigt, FractionalSpecial, PowerLawPlateau

export null_modulus

######################################################
# datagen.jl
export stepgen, rampgen, singen, repeatdata, addnoise

######################################################
# processing.jl
export var_resample, downsample, fixed_resample, smooth, mapbackdata

export modelfit, modelpredict, modelstepfit, modelsteppredict

export savedata, loaddata, savemodel, loadmodel, exportdata

######################################################
include("base.jl")
include("definitions.jl")
include("models.jl")
include("datagen.jl")
include("processing.jl")
######################################################
end
