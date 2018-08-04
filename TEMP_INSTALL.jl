#!/usr/bin/env julia

Pkg.update()

# From rheos-cambridge forked repositories (so that versions are fixed and controlled)
try
    Pkg.clone("https://github.com/rheos-cambridge/MittagLeffler.jl")
catch
    # if already installed 
    LoadError
end

try
    Pkg.clone("https://github.com/rheos-cambridge/FastConv.jl")
catch
    # if already installed
    LoadError
end
