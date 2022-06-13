using VortexDance
using Documenter

DocMeta.setdocmeta!(VortexDance, :DocTestSetup, :(using VortexDance); recursive=true)

makedocs(;
    modules=[VortexDance],
    authors="Marko <mbudisic@clarkson.edu> and contributors",
    repo="https://github.com/mbudisic/VortexDance.jl/blob/{commit}{path}#{line}",
    sitename="VortexDance.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://mbudisic.github.io/VortexDance.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/mbudisic/VortexDance.jl",
    devbranch="master",
)
