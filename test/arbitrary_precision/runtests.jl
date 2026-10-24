# Arbitrary-precision acceptance tests of the ISEM route against the Black Hole Perturbation Toolkit
# Teukolsky paclet (MST), opt-in and not part of test/runtests.jl.
#
# For every stored reference case (a, l, m, omega) and every working precision in GSN_AP_BITS the IN
# and UP solutions are computed with the public `Teukolsky_radial`, and the Teukolsky amplitude ratios
# and the radial values at r' = 3.1126842514806725 are compared with the paclet. The pass condition for
# each quantity is
#     digits agreed >= min(N log10(2) - 6 - input conditioning, certified reference digits - 1),
# where the certified reference digits are the agreement between the two stored paclet working
# precisions of the same case (120 and 200 digits; 200 and 400 for omega = 10), and the input
# conditioning is the loss that rounding a and omega to N bits costs: the paclet works with the exact
# rationals, the route with their N-bit roundings, and a strongly damped amplitude (B_inc near a
# quasinormal frequency, a small reflection) is sensitive to that at the level of its condition number.
# It is measured only when the plain requirement fails (one-ulp perturbations of a and of omega, two
# more builds), and printed. The table printed for every case gives the working bits, the digits
# agreed, the digits expected and the time.
#
# Running (the dependency SpinWeightedSpheroidalHarmonics must support BigFloat at complex c, i.e. the
# ArbitraryPrecision branch, v1.4.1 or later):
#     GSN_AP_BITS=113,256 julia --project=<env> test/arbitrary_precision/runtests.jl
#     GSN_AP_BITS=53,113,192,256 GSN_AP_CASES=all julia --project=<env> test/arbitrary_precision/runtests.jl
# GSN_AP_CASES selects "quick" (default: a = 2/5, l = m = 2, omega = 1/2 and 1/2 - 3i/10), "real",
# "complex", "axis", "all", or a comma-separated list of tag substrings. GSN_SOURCE=<checkout> loads the
# package by include as in test/extended.

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

const REFERENCE_DIR = joinpath(@__DIR__, "references")
const RPRIME = 3.1126842514806725      # the binary double, as the paclet's SetPrecision pads it
const ACCEPTANCE_MARGIN_DIGITS = 6     # N log10(2) - 6, as in the brief
const KEYS = ("Incidence/Transmission", "Reflection/Transmission", "R(rprime)/Transmission", "dR(rprime)/Transmission")

# Reference pairs (lower, higher working precision), in order of priority: a case is tested against
# the first pair that holds it. omega = 10 needs the paclet at 200 digits, so it is taken from the
# 200/400 pair; the regular cases from the 120/200 pair; the axis cases (omega = +-delta -+ i sigma,
# delta from 1e-6 to 1e-30, sigma = 1/10, 1/2, 1, 3, and the exact point i/10) from the 120/200 axis
# files, else from 60/120 where the 200-digit run did not finish.
const REFERENCE_PAIRS = (
    ("bhpt_200_a04l2.txt", "bhpt_400_omega10.txt"),
    ("bhpt_200_a09l2.txt", "bhpt_400_omega10.txt"),
    ("bhpt_200_l3m1.txt", "bhpt_400_omega10.txt"),
    ("bhpt_120_a04l2.txt", "bhpt_200_a04l2.txt"),
    ("bhpt_120_a09l2.txt", "bhpt_200_a09l2.txt"),
    ("bhpt_120_l3m1.txt", "bhpt_200_l3m1.txt"),
    ("bhpt_120_axis.txt", "bhpt_200_axis.txt"),
    ("bhpt_60_axis.txt", "bhpt_120_axis.txt"),
)

# Working precisions: integers are BigFloat bits; "double64" runs the route in Double64 (DoubleFloats.jl
# must be in the environment), a fixed type of 106 bits.
bits_list() = [lowercase(strip(s)) == "double64" ? :double64 : parse(Int, s) for s in split(get(ENV, "GSN_AP_BITS", "113,256"), ",")]
# DoubleFloats.jl must be in the environment: the apparent-point test below always runs in Double64.
using DoubleFloats
_precision_bits(bits) = bits === :double64 ? 106 : bits
_precision_label(bits) = bits === :double64 ? "   d64" : @sprintf("%6d", bits)

# Parse the case tag written by bhpt_acceptance_reference.wls into (a, l, m, omega) as exact numbers.
function parse_tag(tag::AbstractString)
    fields = Dict(split(part, "=", limit=2) for part in split(tag, ";"))
    a = _parse_exact(fields["a"])
    l = parse(Int, fields["l"])
    m = parse(Int, fields["m"])
    omega = _parse_exact(fields["omega"])
    return a, l, m, omega
end
function _parse_exact(text::AbstractString)
    # Mathematica InputForm of an exact number: rationals, integers, "x - (y*I)/z", "I/2", "2 - 2*I"
    t = replace(strip(text), " " => "")
    occursin("I", t) || return _parse_rational(t)
    # split at the sign before the imaginary part (the real part may be absent)
    idx = 0
    for i in length(t):-1:2
        if (t[i] == '+' || t[i] == '-') && t[i-1] != '('
            idx = i
            break
        end
    end
    real_text = idx == 0 ? "0" : t[1:idx-1]
    imag_text = idx == 0 ? t : t[idx:end]
    imag_text = replace(imag_text, "*I" => "", "I" => "1")
    imag_text = replace(imag_text, "(" => "", ")" => "")
    return _parse_rational(real_text) + im * _parse_rational(imag_text)
end
function _parse_rational(t::AbstractString)
    t = replace(t, "(" => "", ")" => "")
    occursin("/", t) || return parse(BigInt, t) // 1
    num, den = split(t, "/")
    return parse(BigInt, num) // parse(BigInt, den)      # 10^-30 appears in the axis cases
end

function selected(tag, selection)
    a, l, m, omega = parse_tag(tag)
    # The paclet's principal branch is not on the physical sheet for Re omega < 0 (the -delta side of
    # the negative imaginary axis): those references are kept in the files for the record only.
    real(omega) < 0 && return false
    # The exact points of the negative imaginary axis are on the branch cut: the paclet's principal
    # value there matches its +delta side in the amplitudes but not in the Up radial values (another
    # branch of the Up normalization), so those records are diagnostics, not references; the route gives
    # the +delta limit there, which is tested through the +delta cases.
    iszero(real(omega)) && imag(omega) < 0 && return false
    selection == "all" && return true
    # a comma-separated list of tag substrings selects those cases
    occursin(";", selection) && return any(occursin(strip(part), tag) for part in split(selection, ","))
    selection == "real" && return isreal(omega)
    selection == "complex" && return !isreal(omega)
    selection == "axis" && return !isreal(omega) && abs(real(omega)) <= 1 // 10^6
    # quick
    return a == 2 // 5 && l == 2 && m == 2 && (omega == 1 // 2 || omega == 1 // 2 - 3im // 10)
end

certified_digits(low, high) = (low === nothing || high === nothing) ? NaN : digits_agreed(low, high)

# Digits lost to the rounding of the inputs: the relative change of each quantity under a one-ulp
# change of a and of Re omega, divided by the relative change of the input (log10 of the larger).
function input_conditioning(l, m, ab, wb, bc, rprime, values)
    loss = Dict{String,Float64}()
    perturbed = (
        (nextfloat(ab), wb, Float64((nextfloat(ab) - ab) / ab)),
        (ab, wb isa Complex ? Complex{BigFloat}(nextfloat(real(wb)), imag(wb)) : nextfloat(wb),
            Float64((nextfloat(real(wb)) - real(wb)) / abs(wb))),
    )
    for (a_p, w_p, relative) in perturbed
        R = GSN.Teukolsky_radial(-2, l, m, a_p, w_p, bc)
        Rv = R.Teukolsky_solution(rprime)
        for (quantity, value, new_value) in zip(KEYS, values, (R.incidence_amplitude, R.reflection_amplitude, Rv[1], Rv[2]))
            sensitivity = Float64(abs(Complex{BigFloat}(new_value) - Complex{BigFloat}(value)) / abs(Complex{BigFloat}(value))) / relative
            loss[quantity] = max(get(loss, quantity, 0.0), log10(max(sensitivity, 1.0)))
        end
    end
    return loss
end

function run_case(tag, ref_high, ref_low, bits)
    a, l, m, omega = parse_tag(tag)
    T = bits === :double64 ? Double64 : BigFloat
    setprecision(BigFloat, _precision_bits(bits)) do
        ab = T(numerator(a)) / T(denominator(a))
        wb = isreal(omega) ? T(numerator(real(omega))) / T(denominator(real(omega))) :
            Complex{T}(T(numerator(real(omega))) / T(denominator(real(omega))),
                T(numerator(imag(omega))) / T(denominator(imag(omega))))
        rprime = T(RPRIME)
        expected = _precision_bits(bits) * log10(2) - ACCEPTANCE_MARGIN_DIGITS
        for (bc, key) in ((GSN.IN, "In"), (GSN.UP, "Up"))
            t0 = time()
            R = GSN.Teukolsky_radial(-2, l, m, ab, wb, bc)
            Rv = R.Teukolsky_solution(rprime)
            dt = time() - t0
            values = (R.incidence_amplitude, R.reflection_amplitude, Rv[1], Rv[2])
            conditioning = nothing
            for (quantity, value) in zip(KEYS, values)
                refkey = "$key.$quantity"
                ref = get(ref_high, refkey, nothing)
                ref === nothing && continue
                certified = certified_digits(get(ref_low, refkey, nothing), ref)
                agreed = digits_agreed(Complex{BigFloat}(value), ref)
                required = min(expected, certified - 1)
                loss = 0.0
                if agreed < required
                    conditioning === nothing && (conditioning = input_conditioning(l, m, ab, wb, bc, rprime, values))
                    loss = conditioning[quantity]
                    required = min(expected - loss, certified - 1)
                end
                @printf("%-42s %3s %-26s %s bits  agreed %6.1f  expected %6.1f  conditioning %5.1f  reference %6.1f  %6.2f s\n",
                    tag, key, quantity, _precision_label(bits), agreed, expected, loss, certified, dt)
                # A reference without a lower-precision counterpart is uncertified: reported, not tested.
                isnan(certified) ? (@test_skip agreed >= expected) : (@test agreed >= required)
            end
        end
    end
end

@testset "arbitrary-precision ISEM vs BHPT" begin
    selection = get(ENV, "GSN_AP_CASES", "quick")
    tested = Set{String}()
    for (low_file, high_file) in REFERENCE_PAIRS
        low = load_bhpt_refs(joinpath(REFERENCE_DIR, low_file))
        high = load_bhpt_refs(joinpath(REFERENCE_DIR, high_file))
        for tag in sort(collect(keys(high)))
            selected(tag, selection) || continue
            haskey(low, tag) || continue
            tag in tested && continue
            push!(tested, tag)
            for bits in bits_list()
                run_case(tag, high[tag], low[tag], bits)
            end
        end
    end
end

# A root of eta 5.4e-7 off the real axis inside the infinity march, (s, l, m, a, omega) = (-2, 2, 2, 0.99, 0.75885)
# with a and omega the binary doubles: the Taylor patches about real points collapsed onto it and Double64 failed
# with "noninvertible q0"; the route crosses it with the local series at the apparent point. References: R at
# unit transmission from an independent 256-bit Frobenius and Taylor solution of the Teukolsky equation. Double64 is
# tested whatever GSN_AP_BITS selects.
@testset "apparent point next to the real axis" begin
    references = ((1.6, 0.06194167414236735 + 0.02117128350690939im), (2.5, -0.22435973040728777 + 0.8634188314161765im),
        (4.0, -3.779344698099683 - 2.9519979793236772im), (10.0, -74.90514587608912 - 41.48465882312736im))
    for bits in unique([:double64; bits_list()])
        T = bits === :double64 ? Double64 : BigFloat
        setprecision(BigFloat, _precision_bits(bits)) do
            R = GSN.Teukolsky_radial(-2, 2, 2, T(0.99), T(0.75885), GSN.IN)
            for (r, value) in references
                @test ComplexF64(R.Teukolsky_solution(T(r))[1]) ≈ value rtol=1e-12
            end
        end
    end
end
