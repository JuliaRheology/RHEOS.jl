using Documenter, InverseLaplace

makedocs(
    format = :html,
    sitename = "InverseLaplace.jl",
    modules = [InverseLaplace],
    pages = [
        "index.md"
    ]
)

deploydocs(
    repo = "github.com/jlapeyre/InverseLaplace.jl.git",
    target = "build",
    julia  = "0.6",
    deps = nothing,
    make = nothing
)

# using Documenter, InverseLaplace

# makedocs(
#     modules = InverseLaplace,
# #    clean   = false,
# )

# # deploydocs(
# #     deps = Deps.pip("pygments", "mkdocs", "mkdocs-material"),
# #     repo = "github.com/MichaelHatherly/PrivateModules.jl.git",
# # )
