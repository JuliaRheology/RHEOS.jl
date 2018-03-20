#!/usr/bin/env julia
include("../src/RHEOS.jl")
using RHEOS

filedir = "../data/afmData1.txt"

data = AFMfileload(filedir, "strlx")

