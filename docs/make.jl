using Documenter
using SargassumBOMB

makedocs(
    sitename = "SargassumBOMB.jl",
    authors = "Gage Bonner",
    format = Documenter.HTML(),
    modules = [SargassumBOMB],
    checkdocs = :export,
    warnonly = true,
    pages = [
        "Introduction" => "index.md",
        "First Steps" => "first-steps.md",
        "Interpolants" => "interpolants.md",
        "Theory" => "theory.md",
        "API" => "api.md",
        "Syntax" => "syntax.md"
    ]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
