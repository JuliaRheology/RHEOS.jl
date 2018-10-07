#!/usr/bin/env julia
__precompile__(true)

module RHEOS

# bundled dependencies from rheos-cambridge forked repos
cd("..")
MittLeffLiteDir = joinpath(@__DIR__, "..", "deps", "MittLeffLite", "src", "MittLeffLite.jl")
include(MittLeffLiteDir)
using MittLeffLite
FastConvDir = joinpath(@__DIR__, "..", "deps", "FastConv", "src", "FastConv.jl")
include(FastConvDir)
using FastConv

# installed from Julia package repository
using InverseLaplace
using uCSV
using ImageFiltering
using Interpolations
using NLopt
using DataFrames

import SpecialFunctions.gamma

# Base and stdlib imports
import Base: +, -, *
import Random

######################################################
# debug
export trapz, G_springpot, G_spring
export boltzintegral_sing, boltzintegral_nonsing
export boltzconvolve

# utility
export closestindex, deriv
# definitions.jl
export RheologyData, RheologyModel, RheologyModelTemp, RheologyDynamic
# IO.jl
export fileload #, savedata, loaddata, savemodel, loadmodel, exportdata
# models.jl
export null_modulus
export SpringPot, Spring, Dashpot
export FractionalMaxwell, FractionalMaxwellSpring, FractionalMaxwellDashpot, Maxwell
export FractionalKelvinVoigt, FractionalKVspring, FractionalKVdashpot, KelvinVoigt
export FractionalZener, FractionalSLS, SLS
export FractionalSpecial
export JeffreysPT
export SLS2, PowerLawPlateau
# datagen.jl
export stepgen, rampgen, singen, repeatdata, addnoise
# processing.jl
export var_resample, downsample, fixed_resample, smooth, mapbackdata
export modelfit, modelpredict, modelstepfit, modelsteppredict
export dynamicmodelfit, dynamicmodelpredict

include("base.jl")
include("definitions.jl")
include("IO.jl")
include("modeldatabase.jl")
include("datagen.jl")
include("processing.jl")

######################################################

end
