using VortexDynamics2D
using Documenter

DocMeta.setdocmeta!(VortexDynamics2D, :DocTestSetup, :(using VortexDynamics2D); recursive=true)

makedocs(;
    modules=[VortexDynamics2D],
    authors="Marko <mbudisic@clarkson.edu> and contributors",
    repo="https://github.com/mbudisic/VortexDynamics2D.jl/blob/{commit}{path}#{line}",
    sitename="VortexDynamics2D.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://mbudisic.github.io/VortexDynamics2D.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/mbudisic/VortexDynamics2D.jl",
    devbranch="master",
)
