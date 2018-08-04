# using Documenter, RHEOS

makedocs(modules=[RHEOS],
         doctest=false, clean=true,
         format =:html,
         sitename="RHEOS.jl",
         authors="Louis Kaplan",
         pages = [
             "Home" => "index.md"
         ]
         )

deploydocs(
    repo = "github.com:rheos-cambridge/RHEOS.jl.git",
    julia = "0.6.3",
    deps = nothing,
    make = nothing,
    target = "build"
)