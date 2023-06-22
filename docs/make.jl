using Pkg; Pkg.add("Documenter")
using Documenter, GeneralizedSasakiNakamura

makedocs(
    sitename="GeneralizedSasakiNakamura.jl",
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true"
    )
)

deploydocs(
    repo = "github.com/ricokaloklo/GeneralizedSasakiNakamura.jl.git",
)
