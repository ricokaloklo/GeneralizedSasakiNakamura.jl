pushfirst!(LOAD_PATH, joinpath(@__DIR__, ".."))

using Documenter, GeneralizedSasakiNakamura

makedocs(
    sitename = "GeneralizedSasakiNakamura.jl",
    modules = [GeneralizedSasakiNakamura],
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        edit_link = nothing,
    ),
    pages = [
        "Home" => "index.md",
        "Examples" => "examples.md",
        "API Reference" => "APIs.md",
    ],
    checkdocs = :none,
)

deploydocs(
    repo = "github.com/ricokaloklo/GeneralizedSasakiNakamura.jl.git",
    devbranch = "master",
)
