using DynamicTimeseries
using Documenter

DocMeta.setdocmeta!(DynamicTimeseries, :DocTestSetup, :(using DynamicTimeseries); recursive=true)

makedocs(;
    modules=[DynamicTimeseries],
    authors="Galen Lynch <galen@galenlynch.com>",
    sitename="DynamicTimeseries.jl",
    format=Documenter.HTML(;
        canonical="https://galenlynch.github.io/DynamicTimeseries.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/galenlynch/DynamicTimeseries.jl",
    devbranch="main",
)
