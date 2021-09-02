using BetaDispersion
using Documenter

DocMeta.setdocmeta!(BetaDispersion, :DocTestSetup, :(using BetaDispersion); recursive=true)

makedocs(;
    modules=[BetaDispersion],
    authors="Arthur Newbury",
    repo="https://github.com/EvoArt/BetaDispersion.jl/blob/{commit}{path}#{line}",
    sitename="BetaDispersion.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://EvoArt.github.io/BetaDispersion.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/EvoArt/BetaDispersion.jl",
)
