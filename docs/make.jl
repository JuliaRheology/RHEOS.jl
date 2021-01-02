using Documenter, RHEOS
###
# convert Jupyter notebooks to markdown and do some simple preprocessing before
# feeding into Documenter.jl
function convert_to_markdown(file)
    run(`jupyter nbconvert examples/$file --to markdown --template docs/documenter.tpl --output-dir docs/src-staging`)
    return "docs/src-staging/$(replace(file, "ipynb"=>"md"))"
end

function convert_equations(file)
    contents = read(file, String)
    contents = replace(contents, r"\$\$(.*?)\$\$"s => s"""```math
    \g<1>
    ```""")
    contents = replace(contents, r"\* \$(.*?)\$" => s"* ``\g<1>``") # starting a line with inline math screws up tex2jax for some reason
    write(file, contents)
    return file
end

rm("docs/src-staging", force=true, recursive=true)
mkdir("docs/src-staging")
# println(readdir("docs"))
# println(pwd())
for file in readdir("examples")
    if endswith(file, "ipynb")
        # println(file)
        file |> convert_to_markdown |> convert_equations
    elseif !startswith(file, ".")
        cp("examples/$file", "docs/src-staging/$file")
    end
end
cp("docs/src/index.md","docs/src-staging/index.md")
cp("docs/src/API.md","docs/src-staging/API.md")
cp("docs/src/architecture.md","docs/src-staging/architecture.md")

mkdir("docs/src-staging/assets")
cp("docs/src/assets/logo.png", "docs/src-staging/assets/logo.png")
cp("docs/src/assets/diagram.svg", "docs/src-staging/assets/diagram.svg")

# highlight output cells (i.e. anything withouout a language specified) white
@eval Documenter.Writers.HTMLWriter function mdconvert(c::Markdown.Code, parent::MDBlockContext; kwargs...)
    @tags pre code
    language = isempty(c.language) ? "none" : c.language
    pre[".language-$(language)"](code[".language-$(language)"](c.code))
end

##
makedocs(modules=[RHEOS],
         doctest=false, clean=true,
         format =Documenter.HTML(),
         sitename=" ",
         source = "src-staging",
         authors="Louis Kaplan",
         pages = [
         "Home" => "index.md",
	     "Architecture" => "architecture.md",
	     "File I/O" => "fileIO.md",
	     "Preprocessing" => "preprocessing.md",
	     "Generating Data" => "generatingdata.md",
	     "Fit and Predict" => ["Time data" => "fitpredictTime.md",
				               "Frequency data" => "fitpredictFreq.md"],
         "Models" => ["Basic Elements" => "elements.md",
                      "Maxwell" => "fractionalMaxwell.md",
			          "Kelvin-Voigt" => "fractionalKelvinVoigt.md",
			          "Zener" => "fractionalZener.md",
			          "Poynting-Thomson" => "fractionalPT.md",
		              "Burgers" => "burgers.md",
			          "Create Your Model" => "createmodel.md"],
	    "Examples" => "examples.md",
         "API" => "API.md"
         ]
         )

deploydocs(
    repo = "github.com/JuliaRheology/RHEOS.jl.git",
    deps = nothing,
    make = nothing,
    target = "build"
)
