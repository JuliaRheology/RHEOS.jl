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
import Base: +, -, *

######################################################

# definitions.jl
export RheologyData, RheologyModel, RheologyModelTemp
# IO.jl
export fileload, savedata, loaddata, savemodel, loadmodel, exportdata
# models.jl
export null_modulus
export SLS, SpringPot, FractionalMaxwell, FractionalKelvinVoigt, FractionalSpecial, PowerLawPlateau
# datagen.jl
export stepgen, rampgen, singen, repeatdata, addnoise
# processing.jl
export var_resample, downsample, fixed_resample, smooth, mapbackdata
export modelfit, modelpredict, modelstepfit, modelsteppredict

include("base.jl")
include("definitions.jl")
include("IO.jl")
include("models.jl")
include("datagen.jl")
include("processing.jl")

######################################################

end
