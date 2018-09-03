#!/usr/bin/env julia
__precompile__()

module RHEOS
# installed using Pkg.clone from rheos-cambridge forked repos
using MittLeffLite
using FastConv
# install using Pkg.add from Julia central package repository
using InverseLaplace
using uCSV
using ImageFiltering
using Interpolations
using NLopt
using JLD
using DataFrames
import Base: +, -, *

######################################################
# debug
export trapz, G_springpot, G_spring
export boltzintegral_sing, boltzintegral_nonsing
export boltzconvolve_sing, boltzconvolve_nonsing

# utility
export closestindex, deriv
# definitions.jl
export RheologyData, RheologyModel, RheologyModelTemp
# IO.jl
export fileload, savedata, loaddata, savemodel, loadmodel, exportdata
# models.jl
export null_modulus
export SpringPot, Spring, Dashpot, SLS, FractionalSLS, FractionalKVspring, PowerLawPlateau
export FractionalMaxwell, FractionalMaxwellSpring, FractionalMaxwellDashpot, Maxwell
export FractionalSpecial
export SLS2
#  SLS, FractionalMaxwell, FractionalKelvinVoigt, FractionalSpecial, PowerLawPlateau
# datagen.jl
export stepgen, rampgen, singen, repeatdata, addnoise
# processing.jl
export var_resample, downsample, fixed_resample, smooth, mapbackdata
export modelfit, modelpredict, modelstepfit, modelsteppredict

include("base.jl")
include("definitions.jl")
include("IO.jl")
include("modeldatabase.jl")
include("datagen.jl")
include("processing.jl")

######################################################

end
