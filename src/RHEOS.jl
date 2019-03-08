#!/usr/bin/env julia
__precompile__(true)

module RHEOS

RheoFloat = Float32
# installed from Julia package repository
using InverseLaplace
using NLopt
using JLD2
using Interpolations
import DSP.conv
# Base and stdlib imports
using Base.Cartesian
import Base.eltype
import SpecialFunctions.gamma
import Base: +, -, *
import DelimitedFiles: readdlm, writedlm
export constantcheck

######################################################
# definitions.jl
export RheologyData, RheologyModel, RheologyDynamic
# IO.jl
export importdata, exportdata, savedata, loaddata, savemodel, loadmodel
# models.jl
export null_modulus
export SpringPot, Spring, DashPot
export FractionalMaxwell, FractionalMaxwellSpring, FractionalMaxwellDashpot, Maxwell
export FractionalKelvinVoigt, FractionalKVspring, FractionalKVdashpot, KelvinVoigt
export FractionalZener, FractionalSLS, SLS
export FractionalSpecial
export JeffreysPT
export SLS2, PowerLawPlateau
# datagen.jl
export linegen, stepgen, rampgen, singen, noisegen, repeatdata
# processing.jl
export variableresample, downsample, fixedresample, smooth, zerotime
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
include("advance.jl")
######################################################
end
