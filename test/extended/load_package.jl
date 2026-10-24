# Load GeneralizedSasakiNakamura for the extended suite and bind it to `GSN`.
#
# * Without GSN_SOURCE: `using GeneralizedSasakiNakamura` from the active environment (Pkg.test or any
#   environment that has the package).
# * With GSN_SOURCE=<source checkout>: include "$GSN_SOURCE/src/GeneralizedSasakiNakamura.jl" directly
#   (GSN_SOURCE_FILE overrides the file). GSN_SWSH_PATH optionally selects a dependency checkout.
#
# Included at top level of Main by runtests.jl, and by each group file when that file is run on its own.

using Test

const _GSN_SOURCE_DIR = get(ENV, "GSN_SOURCE", "")

if isempty(_GSN_SOURCE_DIR)
    using GeneralizedSasakiNakamura
else
    using LinearAlgebra
    # Dependencies used by include-style loading.
    using ApproxFun, AbstractTrees, LogarithmicNumbers, SpecialFunctions
    let swsh = get(ENV, "GSN_SWSH_PATH", "")
        if !isempty(swsh) && isdir(swsh)
            pushfirst!(LOAD_PATH, swsh)
        end
    end
    BLAS.set_num_threads(1)
    include(get(ENV, "GSN_SOURCE_FILE",
        joinpath(_GSN_SOURCE_DIR, "src", "GeneralizedSasakiNakamura.jl")))
    using .GeneralizedSasakiNakamura
end

const GSN = GeneralizedSasakiNakamura

println("GSN extended tests: package = ", pathof(GSN) === nothing ?
    (isempty(_GSN_SOURCE_DIR) ? "?" : _GSN_SOURCE_DIR) : pathof(GSN),
    "; Julia ", VERSION, "; threads ", Threads.nthreads())
for dependency in (GSN.KerrGeodesics, GSN.SpinWeightedSpheroidalHarmonics)
    println("  ", nameof(dependency), " ", pkgversion(dependency))
end
flush(stdout)
