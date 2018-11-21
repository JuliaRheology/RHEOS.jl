#!/usr/bin/env julia
__precompile__(true)

module RHEOS

# installed from Julia package repository
using InverseLaplace
using NLopt
using JLD2
import DSP.conv
# Base and stdlib imports
using Base.Cartesian
import SpecialFunctions.gamma
import Base: +, -, *
import DelimitedFiles: readdlm, writedlm

######################################################
# definitions.jl
export RheologyData, RheologyModel, RheologyDynamic
# IO.jl
export fileload, savedata, loaddata, savemodel, loadmodel, exportdata
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
export var_resample, downsample, fixed_resample, smooth, mapbackdata, zerotime
export modelfit, modelpredict, modelstepfit, modelsteppredict
export dynamicmodelfit, dynamicmodelpredict

######################################################
# bundled dependencies from rheos-cambridge forked repos
MittLeffLiteDir = joinpath(@__DIR__, "..", "deps", "MittLeffLite", "MittLeffLite.jl")
include(MittLeffLiteDir)
FastConvDir = joinpath(@__DIR__, "..", "deps", "FastConv", "FastConv.jl")
include(FastConvDir)

include("base.jl")
include("definitions.jl")
include("IO.jl")
include("modeldatabase.jl")
include("datagen.jl")
include("processing.jl")
######################################################
end
