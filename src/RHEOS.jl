#!/usr/bin/env julia
__precompile__(true)

module RHEOS

# installed from Julia package repository
using NLopt
using JLD2
using DataStructures
using FunctionWrappers: FunctionWrapper
using Dierckx

# useful for the various model functions
using InverseLaplace
import MittagLeffler: mittleff as mittlefforiginal
import SpecialFunctions: gamma

import DSP.conv

# Base and stdlib imports
import Base: +, -, *,|
import DelimitedFiles: readdlm, writedlm

######################################################################

# This defines the data type for all arrays, parameters and processing
# it is defined as a const to avoid performance penalties.
# See julia docs for more info on this.
const RheoFloat = Float64

# Singularity approximation constant
const singularity_offset = 10.0

# convenience data types used as in many function parameters
const IntOrNone = Union{Integer, Nothing}
const RheovecOrNone = Union{Vector{RheoFloat}, Nothing}

######################################################################
export RheoFloat

# definitions.jl
export RheoLogItem, RheoLog, rheologrun, showlog
export RheoTimeData, RheoFreqData
export rheotimedatatype, rheofreqdatatype, check_time_data_consistency
export LoadingType, strain_imposed, stress_imposed
export hastime, hasstress, hasstrain
export TimeDataType, time_only, strain_only, stress_only, strain_and_stress
export FreqDataType, invalid_freq_data, freq_only, with_modulus
export RheoModelClass, RheoModel, model_parameters
export rheoconv, invLaplace
export freeze_params
export relaxmod, creepcomp, storagemod, lossmod, dynamicmod

# IO.jl
export importcsv, exportcsv#, savedata, loaddata, savemodel, loadmodel

# models.jl
export null_modulus
export Springpot, Spring, Dashpot
export Fract_Maxwell, FractS_Maxwell, FractD_Maxwell, Maxwell
export Fract_KelvinVoigt, FractS_KelvinVoigt, FractD_KelvinVoigt, KelvinVoigt
export Fract_Zener, FractSLS_Zener, SLS_Zener, FractJeffreys_Zener, Jeffreys_Zener, FractSolid
export Fract_PT, FractSLS_PT, SLS_PT, Jeffreys_PT, FractJeffreys_PT
export BurgersLiquid
export SLS2, PowerLawPlateau

# datagen.jl
export timeline
export strainfunction, stressfunction
export hstep, ramp, stairs, square, sawtooth, triangle
export frequencyspec

# processing.jl
export resample, indexweight, cutting, smooth, extract
export modelfit, modelpredict, modelstepfit, modelsteppredict
export dynamicmodelfit, dynamicmodelpredict

#interface.jl
export Interface
export AFM, Tweezers

######################################################
# bundled dependencies from rheos-cambridge forked repos
#MittLeffLiteDir = joinpath(@__DIR__, "..", "deps", "MittLeffLite", "MittLeffLite.jl")
#include(MittLeffLiteDir)

include("base.jl")
include("definitions.jl")
include("IO.jl")
include("modeldatabase.jl")
include("datagen.jl")
include("processing.jl")
include("interface.jl")

end
