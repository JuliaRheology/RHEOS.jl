#!/usr/bin/env julia
include("../src/RHEOS.jl")
using RHEOS
using PyPlot

fdir = "../data/afmData1.txt"

data = AFMfileload(fdir, "strlx"; cpfind = "hertz", R = 25e-6)

data.insight = true

p0 = [12363.7, 0.5]
lb = [0.0, 0.01]
ub = [20000.0, 0.99]
modelfit!(data, "springpot", p0, lb, ub)
fiteval(data, "springpot")