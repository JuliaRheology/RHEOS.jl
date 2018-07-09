using Documenter, RHEOS

makedocs(modules=[RHEOS],
         doctest=false, clean=true,
         format =:html,
         sitename="RHEOS.jl",
         authors="Louis Kaplan, Alessandra Bonfanti, Alexandre Kabla"
         )

# deploydocs(repo = "github.com:rheos-cambridge/RHEOS.jl.git",
#            target = "build",
#            osname = "windows",
#            julia = "0.6.3",
#            deps = nothing,
#            make = nothing
#            )