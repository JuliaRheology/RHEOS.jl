using Documenter, RHEOS

makedocs(modules=[RHEOS],
         doctest=false, clean=true,
         format =:html,
         sitename="RHEOS.jl", # try making this a blank string!
         authors="Louis Kaplan",
         pages = [
             "Home" => "index.md",
             "Fitting Data" => "fittingdata.md",
             "Predicting Responses" => "predictingresponse.md",
             "Generating Loading - Under Construction" => "generatingdata.md",
             "Preprocessing Tools - Under Construction" => "preprocessing.md",
             "File I/O" => "fileIO.md",
             "Models" => "models.md",
             "API" => "API.md"
         ]
         )

deploydocs(
    repo = "github.com/JuliaRheology/RHEOS.jl.git",
    deps = nothing,
    make = nothing,
    target = "build"
)
