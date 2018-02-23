#!/usr/bin/env julia
include("../src/RHEOS.jl")
using RHEOS

filedir = "../data/rheologyData1.csv"

data = fileload(filedir, ["time","stress","strain"])

data.insight = true

fixed_resample!(data, [1,450],[8],["up"])
# downsample!(data, [1,450], [2])
# var_resample!(data, :σᵦ, 0.1, _mapback = true)

# SLS plot
params = [843.83, 2037.72, 5.17872]
