using GLTimeseries
using Documenter

DocMeta.setdocmeta!(GLTimeseries, :DocTestSetup, :(using GLTimeseries); recursive=true)

makedocs(;
    modules=[GLTimeseries],
    authors="Galen Lynch <galen@galenlynch.com>",
    sitename="GLTimeseries.jl",
    format=Documenter.HTML(;
        canonical="https://galenlynch.github.io/GLTimeseries.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/galenlynch/GLTimeseries.jl",
    devbranch="main",
)
