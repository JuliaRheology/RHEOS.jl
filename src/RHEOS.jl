#!/usr/bin/env julia
__precompile__(true)

module RHEOS

# installed from Julia package repository
using InverseLaplace
using NLopt
using JLD2
using DataStructures
using FunctionWrappers: FunctionWrapper
import DSP.conv
import SpecialFunctions.gamma
# Base and stdlib imports
import Base: +, -, *
import DelimitedFiles: readdlm, writedlm

######################################################################

# This defines the data type for all arrays, parameters and processing
# it is defined as a const to avoid performance penalties.
# See julia docs for more info on this.
const RheoFloat = Float64

# convenience data types used as in many function parameters
const IntOrNone = Union{Integer, Nothing}
const RheovecOrNone = Union{Vector{RheoFloat}, Nothing}

######################################################################
export RheoFloat

# definitions.jl
export RheoLogItem, RheoLog, rheologrun
export RheoTimeData, RheoTimeDataType, RheoFreqData, RheoFreqDataType, check_time_data_consistency
export RheoModelClass, RheoModel, model_parameters
export LoadingType, strain_imposed, stress_imposed
export TimeDataType, time_only, strain_only, stress_only, strain_and_stress
export FreqDataType, invalid_freq_data, freq_only, with_modulus
export rheoconv,invLaplace
export freeze_params

# IO.jl
export importcsv, exportcsv#, savedata, loaddata, savemodel, loadmodel

# models.jl
export null_modulus
export Springpot, Spring, Dashpot
export FractionalMaxwell, FractionalMaxwellSpring, FractionalMaxwellDashpot, Maxwell
export FractionalKelvinVoigt, FractionalKVspring, FractionalKVdashpot, KelvinVoigt
export FractionalZener, FractionalSLS, SLS, FractionalJeffreys, Jeffreys
export FractionalSpecial
export Jeffreys_PT, SLS_PT
export SLS2, PowerLawPlateau
export BurgersLiquid

# datagen.jl
export timeline
export strainfunction, stressfunction
export hstep, ramp, stairs, square, sawtooth, triangle
export frequencyspec

# processing.jl
export resample, cutting, smooth, extract
export modelfit, modelpredict, modelstepfit, modelsteppredict
export dynamicmodelfit, dynamicmodelpredict

######################################################
# bundled dependencies from rheos-cambridge forked repos
MittLeffLiteDir = joinpath(@__DIR__, "..", "deps", "MittLeffLite", "MittLeffLite.jl")
include(MittLeffLiteDir)

include("base.jl")
include("definitions.jl")
include("IO.jl")
include("modeldatabase.jl")
include("datagen.jl")
include("processing.jl")
######################################################
end
