using NetworkCascadesMCMC
using Documenter

DocMeta.setdocmeta!(NetworkCascadesMCMC, :DocTestSetup, :(using NetworkCascadesMCMC); recursive=true)

makedocs(;
    modules=[NetworkCascadesMCMC],
    authors="BrandnerN",
    repo="https://github.com/kafu16/NetworkCascadesMCMC.jl/blob/{commit}{path}#{line}",
    sitename="NetworkCascadesMCMC.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://kafu16.github.io/NetworkCascadesMCMC.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/kafu16/NetworkCascadesMCMC.jl",
)
