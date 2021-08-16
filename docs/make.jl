using BetaDisp
using Documenter

DocMeta.setdocmeta!(BetaDisp, :DocTestSetup, :(using BetaDisp); recursive=true)

makedocs(;
    modules=[BetaDisp],
    authors="Arthur Newbury",
    repo="https://github.com/EvoArt/BetaDisp.jl/blob/{commit}{path}#{line}",
    sitename="BetaDisp.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://EvoArt.github.io/BetaDisp.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/EvoArt/BetaDisp.jl",
)
