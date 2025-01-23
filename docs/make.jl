using LUSE_ENGR705_FiniteElementAnalysis
using Documenter

DocMeta.setdocmeta!(LUSE_ENGR705_FiniteElementAnalysis, :DocTestSetup, :(using LUSE_ENGR705_FiniteElementAnalysis); recursive=true)

makedocs(;
    modules=[LUSE_ENGR705_FiniteElementAnalysis],
    authors="Joby M. Anthony III",
    repo="https://github.com/jmanthony3/LUSE_ENGR705_FiniteElementAnalysis.jl/blob/{commit}{path}#{line}",
    sitename="LUSE_ENGR705_FiniteElementAnalysis.jl",
    doctest=false,
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://jmanthony3.github.io/LUSE_ENGR705_FiniteElementAnalysis.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/jmanthony3/LUSE_ENGR705_FiniteElementAnalysis.jl",
    devbranch="main",
)
