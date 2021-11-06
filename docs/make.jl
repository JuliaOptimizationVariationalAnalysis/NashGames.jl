using NashGames
using Documenter

DocMeta.setdocmeta!(NashGames, :DocTestSetup, :(using NashGames); recursive=true)

makedocs(;
    modules=[NashGames],
    authors="Tangi Migot tangi.migot@gmail.com",
    repo="https://github.com/JuliaOptimizationVariationalAnalysis/NashGames.jl/blob/{commit}{path}#{line}",
    sitename="NashGames.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://JuliaOptimizationVariationalAnalysis.github.io/NashGames.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/JuliaOptimizationVariationalAnalysis/NashGames.jl",
    devbranch="main",
)
