using Documenter, RHEOS

makedocs(modules=[RHEOS],
         doctest=false, clean=true,
         format =:html,
         sitename=" ",
         authors="Louis Kaplan",
         pages = [
             "Home" => "index.md",
            #  "Fitting Data" => "fittingdata.md",
            #  "Predicting Responses" => "predictingresponse.md",
            #  "Generating Data" => "generatingdata.md",
            #  "Sampling and Filtering" => "samplingandfiltering.md",
            #  "File I/O" => "fileIO.md",
            #  "Models" => "models.md",
             "API" => "API.md"
         ]
         )

deploydocs(
    repo = "github.com/JuliaRheology/RHEOS.jl.git",
    deps = nothing,
    make = nothing,
    target = "build"
)
