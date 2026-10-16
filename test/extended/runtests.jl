# Optional suite; the registered test/runtests.jl does not include this file.
using Test

const TEST_GROUPS = (
    "api", "homogeneous", "isem", "qnm", "pointparticle_mode", "flux", "regression", "threads",
)

selection = get(ENV, "GSN_TEST_GROUPS", haskey(ENV, "GSN_TEST_LEVEL") ? "full" : "")
if isempty(selection)
    println("Set GSN_TEST_GROUPS=quick, full, list, or a comma-separated group list.")
elseif selection == "list"
    println(join(TEST_GROUPS, "\n"))
else
    selected = selection in ("full", "quick") ? TEST_GROUPS :
        Tuple(String.(strip.(split(selection, ','))))
    all(group -> group in TEST_GROUPS, selected) ||
        throw(ArgumentError("Unknown test selection: $selection"))
    length(unique(selected)) == length(selected) ||
        throw(ArgumentError("Repeated test group: $selection"))
    const GSN_TEST_LEVEL_SELECTED = Symbol(get(ENV, "GSN_TEST_LEVEL",
        selection == "quick" ? "quick" : "full"))
    GSN_TEST_LEVEL_SELECTED in (:quick, :full, :extended) ||
        throw(ArgumentError("GSN_TEST_LEVEL must be quick, full or extended."))

    include(joinpath(@__DIR__, "load_package.jl"))
    println("test selection: ", selection, "; level: ", GSN_TEST_LEVEL_SELECTED)
    const _GSN_USING = Expr(:using, Expr(:., fullname(GSN)...))
    @testset "GSN extended suite" begin
        for group in selected
            scope = Module(gensym(Symbol("GSNTest_", group)))
            Core.eval(scope, :(using Test))
            Core.eval(scope, _GSN_USING)
            Core.eval(scope, :(const GSN = $GSN))
            Core.eval(scope, :(const GSN_TEST_LEVEL = $(QuoteNode(GSN_TEST_LEVEL_SELECTED))))
            Core.eval(scope, :(include(path) = Base.include($scope, path)))
            println("Running ", group)
            flush(stdout)
            @time Base.include(scope, joinpath(@__DIR__, group, "runtests.jl"))
            flush(stdout)
        end
    end
end
