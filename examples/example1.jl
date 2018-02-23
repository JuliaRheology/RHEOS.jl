#!/usr/bin/env julia
include("../src/RHEOS.jl")
using RHEOS

filedir = "../data/rheologyData1.csv"

data = fileload(filedir, ["time","stress","strain"])

# data.insight = true

fixed_resample!(data, [1,450],[8],["up"])
# downsample!(data, [1,450], [2])
# var_resample!(data, :σᵦ, 0.1, _mapback = true)

# SLS fit
lb = [0.0, 0.0, 0.0]
ub = [1e5, 1e5, 1e5]
p0 = [1000.0, 1000.0, 100.0]
modelfit!(data, "strlx", "SLS", p0, lb, ub)

# Spring-pot fit
lb = [0.0, 0.0]
ub = [1e5, 1.0]
p0 = [1000.0, 0.5]
modelfit!(data, "strlx", "springpot", p0, lb, ub; singularity = true)

fiteval(data, "SLS", "strlx")
fiteval(data, "springpot", "strlx"; singularity = true)
