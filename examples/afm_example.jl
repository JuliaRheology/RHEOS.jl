#!/usr/bin/env julia
include("../src/RHEOS.jl")
using RHEOS
using PyPlot

filedir = "../data/afmData1.txt"

data = AFMfileload(filedir, "strlx"; cpfind = "threshold", cp_param = 1e-8)

plot(data.t, data.f)
plot(data.t, data.δ, "--")
show()