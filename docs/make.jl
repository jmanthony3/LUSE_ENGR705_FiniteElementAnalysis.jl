using FiniteElementAnalysis
using Documenter

DocMeta.setdocmeta!(FiniteElementAnalysis, :DocTestSetup, :(using FiniteElementAnalysis); recursive=true)

makedocs(;
    modules=[FiniteElementAnalysis],
    authors="Joby M. Anthony III",
    repo="https://github.com/jmanthony3/FiniteElementAnalysis.jl/blob/{commit}{path}#{line}",
    sitename="FiniteElementAnalysis.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://jmanthony3.github.io/FiniteElementAnalysis.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/jmanthony3/FiniteElementAnalysis.jl",
    devbranch="main",
)
