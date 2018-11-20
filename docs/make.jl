using Documenter, RHEOS

makedocs(modules=[RHEOS],
         doctest=false, clean=true,
         format =:html,
         sitename="RHEOS.jl", # try making this a blank string!
         authors="Louis Kaplan",
         pages = [
             "Home" => "01_index.md",
             "Fitting Data" => "02_fittingdata.md",
             "Predicting Responses" => "03_predictingresponse.md",
             "Generating Loading" => "04_generatingdata.md",
             "Pre/Post-Processing" => "05_pre_post_processing.md"
             "Models" => "06_models.md"
             "API" => "07_API.md"
         ]
         )

deploydocs(
    repo = "github.com/JuliaRheology/RHEOS.jl.git",
    deps = nothing,
    make = nothing,
    target = "build"
)
