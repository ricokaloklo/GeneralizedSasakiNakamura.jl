# Arbitrary-precision quasinormal-mode acceptance tests against the independent reference roots in
# references/qnm_reference.txt (s = -2, l = m = 2, a = 2/5 and 9/10, n = 0..9 plus the second n = 8
# root "8_1"; each root was obtained twice on another machine, by the MST incidence-amplitude zero of the
# QNM paclet at 120 digits and by a Leaver continued fraction at 100 digits, and the two agree to 80-89
# digits; the frequency derivative of the Teukolsky incidence amplitude was taken there by fourth-order
# central differences at 120 digits with two step sizes, which agree to 38-45 digits). Opt-in, not part
# of test/runtests.jl.
#
# For every root and every working precision in GSN_AP_BITS, `QNM.qnm_arbitrary_precision` polishes the
# reference's binary64 seed (so that the same root is found, including 8_1) and the frequency and the
# Teukolsky-normalized incidence derivative are compared with the reference:
#     digits agreed >= min(N log10(2) - 6, certified reference digits - 1).
#
#     GSN_AP_BITS=113,256 julia --project=<env> test/arbitrary_precision/qnm_runtests.jl
# GSN_AP_QNM_CASES selects "quick" (default: a = 2/5, n = 0 and 8_1), "all", or a comma-separated list
# of tag substrings such as "a=2/5;l=2;m=2;n=5,a=9/10;l=2;m=2;n=8".

using Test
using Printf

const _AP_SOURCE = get(ENV, "GSN_SOURCE", "")
if isempty(_AP_SOURCE)
    using GeneralizedSasakiNakamura
else
    include(joinpath(_AP_SOURCE, "src", "GeneralizedSasakiNakamura.jl"))
    using .GeneralizedSasakiNakamura
end
const GSN = GeneralizedSasakiNakamura
include(joinpath(@__DIR__, "bhpt_references.jl"))

# GSN_AP_QNM_REFERENCE=<file> uses another reference file of the same format (e.g. a larger benchmark set).
const QNM_REFERENCE = get(ENV, "GSN_AP_QNM_REFERENCE", joinpath(@__DIR__, "references", "qnm_reference.txt"))
const ACCEPTANCE_MARGIN_DIGITS = 6

# Working precisions: integers are BigFloat bits; "double64" polishes in Double64 (DoubleFloats.jl in the environment).
bits_list() = [lowercase(strip(s)) == "double64" ? :double64 : parse(Int, s) for s in split(get(ENV, "GSN_AP_BITS", "113,256"), ",")]
if :double64 in bits_list()
    using DoubleFloats
end
_precision_bits(bits) = bits === :double64 ? 106 : bits

# tag "a=2/5;l=2;m=2;n=8_1" -> (a::Rational, l, m, label)
function parse_qnm_tag(tag::AbstractString)
    fields = Dict(split(part, "=", limit=2) for part in split(tag, ";"))
    num, den = split(fields["a"], "/")
    return parse(Int, num) // parse(Int, den), parse(Int, fields["l"]), parse(Int, fields["m"]), String(fields["n"])
end

# The reference file holds strings and numbers; only the numeric entries are parsed.
function load_qnm_reference(path::AbstractString)
    refs = Dict{String,Dict{String,Any}}()
    for line in eachline(path)
        startswith(line, "#") && continue
        parts = split(line, "|")
        length(parts) == 3 || continue
        tag, key, value = parts
        d = get!(refs, String(tag), Dict{String,Any}())
        d[String(key)] = startswith(strip(value), "\"") || startswith(strip(value), "{") ? String(value) :
            _parse_mathematica_number(value)
    end
    return refs
end

function run_root(tag, ref, bits)
    a, l, m, label = parse_qnm_tag(tag)
    seed = ComplexF64(ref["seed"])
    omega_ref = ref["omega.official_QNM_MST120_cold"]
    derivative_ref = ref["dA_inc.h1e-15"]
    certified_root = Float64(real(ref["digits.independent_cold_MST120_CF100"]))
    certified_derivative = Float64(real(ref["digits.derivative_steps"]))
    expected = _precision_bits(bits) * log10(2) - ACCEPTANCE_MARGIN_DIGITS
    t0 = time()
    result = bits === :double64 ? GSN.QNM.qnm_arbitrary_precision(a, -2, l, m, 0; type=Double64, guess=seed) :
        GSN.QNM.qnm_arbitrary_precision(a, -2, l, m, 0; bits=bits, guess=seed)
    dt = time() - t0
    agreed_root = digits_agreed(Complex{BigFloat}(result.omega), omega_ref)
    agreed_derivative = digits_agreed(Complex{BigFloat}(result.incidence_derivative_teukolsky), derivative_ref)
    @printf("%-22s %s bits  omega agreed %6.1f  expected %6.1f  reference %6.1f | dA_inc/domega agreed %6.1f  reference %6.1f | %d secant steps  %6.1f s\n",
        tag, bits === :double64 ? " d64" : @sprintf("%4d", bits), agreed_root, expected, certified_root, agreed_derivative, certified_derivative, result.iterations, dt)
    @test agreed_root >= min(expected, certified_root - 1)
    @test agreed_derivative >= min(expected, certified_derivative - 1)
end

@testset "arbitrary-precision QNM roots vs independent reference" begin
    selection = get(ENV, "GSN_AP_QNM_CASES", "quick")
    refs = load_qnm_reference(QNM_REFERENCE)
    for tag in sort(collect(keys(refs)))
        a, l, m, label = parse_qnm_tag(tag)
        if selection == "quick"
            a == 2 // 5 && label in ("0", "8_1") || continue
        elseif selection != "all"
            any(occursin(strip(part), tag) for part in split(selection, ",")) || continue
        end
        ref = refs[tag]
        occursin("completed", string(get(ref, "status", ""))) || continue
        for bits in bits_list()
            run_root(tag, ref, bits)
        end
    end
end
