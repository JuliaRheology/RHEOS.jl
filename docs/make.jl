using Documenter
using Literate
using RHEOS

"""
    docprepare()

Convert all source documentation (.jl and .md files)
into markdown, ready for documentation build, and 
relevant .jl files into Jupyter notebooks ready for
deployment to notebooks branch.
"""
function docprepare()
    # remove remnants of previous build and recreate staging dir
    rm("docs/staging-docs", force=true, recursive=true)
    mkdir("docs/staging-docs")

    # copy readme to staging-docs, remove Logo image
    write("docs/staging-docs/index.md", 
            open("README.md") do input
                readuntil(input, "<!-- delim -->", keep = true)
                read(input)
            end)
    
    # copy assets to staging directory
    cp("docs/src/assets", "docs/staging-docs/assets")

    # iterate through src and convert/copy as appropriate
    for file in readdir("docs/src")
        if endswith(file, "md")
            cp("docs/src/$file", "docs/staging-docs/$file")
        elseif endswith(file, "jl")
            Literate.markdown("docs/src/$file", "docs/staging-docs/"; documenter=true)
        end
    end
end

function notebookprepare()
    # create notebook staging dir
    mkdir("docs/staging-docs/notebooks")

    # copy assets to notebooks staging directory
    cp("docs/src/assets", "docs/staging-docs/notebooks/assets")

    for file in readdir("docs/src")
        if endswith(file, "jl")
            Literate.notebook("docs/src/$file", "docs/staging-docs/notebooks/")
        end
    end
end

function maindocbuilder()
    # prepare doc files from source
    docprepare()

    # prepare notebook files from source
    notebookprepare()

    # build docs from staging area
    makedocs(modules=[RHEOS],
            doctest = false, clean=true,
            format = Documenter.HTML(),
            sitename ="RHEOS.jl",
            source = "staging-docs",
            authors = "J Louis Kaplan, Alessandra Bonfati, Alexandre Kabla",
            pages = ["Home" => "index.md",
                     "Architecture" => "architecture.md",
                     "File I/O" => "fileio.md",
                     "Preprocessing" => "preprocessing.md",
                     "Generating Data" => "gendata.md",
                     "Fit and Predict" => ["Time data" => "fitpredictTime.md",
                                           "Frequency data" => "fitpredictFreq.md"],
                     "Numerical Aspects" => "numerical.md",
                     "Models" => ["Basic Elements" => "model_elements.md",
                                  "Maxwell" => "model_maxwell.md",
                                  "Kelvin-Voigt" => "model_kv.md",
                                  "Zener" => "model_zener.md",
                                  "Poynting-Thomson" => "model_pt.md",
                                  "Burgers" => "model_burgers.md",
                                  "Create Your Model" => "model_create.md"],
                     "Additional Examples" => "examples.md",
                     "API" => "API.md"])

    deploydocs(repo = "github.com/JuliaRheology/RHEOS.jl.git",
               deps = nothing,
               make = nothing,
               target = "build")

end

maindocbuilder()
