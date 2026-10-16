"""
    QNMMode(s, l, m, n, branch=:positive_real)

Label one QNM, accepted by [`qnm_frequency`](@ref) and [`qnm_sequence`](@ref)
as an alternative to the plain `(s, l, m, n, branch)` arguments. In this form
`branch=:negative_real` requires an explicit `guess`, from which the root is
polished on the negative-real branch. The plain-argument form with `mirror`
instead returns the reflection `-conj(ω)` of the ordinary `(s, l, -m, n)` root
at the same spin.
"""
struct QNMMode
    s::Int
    l::Int
    m::Int
    n::Int
    branch::Symbol

    function QNMMode(s::Int, l::Int, m::Int, n::Int,
            branch::Symbol=:positive_real)
        l >= max(abs(s), abs(m)) ||
            throw(ArgumentError("l must satisfy l >= max(abs(s), abs(m))."))
        n >= 0 || throw(ArgumentError("n must be nonnegative."))
        branch in (:positive_real, :negative_real) || throw(ArgumentError(
            "branch must be :positive_real or :negative_real."))
        return new(s, l, m, n, branch)
    end
end

function _qnm_convention(value)
    convention = Symbol(lowercase(String(value)))
    convention in (:overtone, :complete_spectrum, :leaver) ||
        throw(ArgumentError(
            "convention must be :overtone, :complete_spectrum, or :leaver."))
    return convention
end

"""
    QNMBranch

Branch selector for [`qnm`](@ref). `ordinary` selects the QNM with positive
real frequency and `mirror` the one with negative real frequency.
`Ordinary` and `Mirror` are capitalised aliases of `ordinary` and `mirror`.
"""
@enum QNMBranch begin
    ordinary
    mirror
end

const Ordinary = ordinary
const Mirror = mirror

"""
    QNMResult

Result of [`qnm`](@ref) when every check passes. The display shows the main
observables; all other fields are available as properties.
"""
struct QNMResult{D}
    data::D
end

"""
    QNMEstimate

Result of [`qnm`](@ref) when all observables are finite but one accuracy check
did not reach its target. The check is named in `limiting_gate`, and
`achieved_relative_precision` gives the precision it did reach (an estimate,
not a bound). Typical cases are a frequency derivative that did not converge to
`1e-8`, or a strongly damped mode whose incidence amplitude does not vanish to
`1e-8` at the root. The values are usable as estimates.
"""
struct QNMEstimate{D}
    data::D
end

"""
    QNMFailure

Result of [`qnm`](@ref) when the root, the radial evaluation or the excitation
factor fails. `stage` and `stop_reason` say which step failed and why; no
amplitudes are returned.
"""
struct QNMFailure{D}
    data::D
end

"""
    QNMEndpointResult

Result of [`qnm`](@ref) for a zero-damping mode (ZDM) at exact extremality,
`a = ±1`. These co-rotating families end at the synchronous frequency
`ω = m a/2`, which is a branch point of the Green function rather than a simple
pole, so the amplitude fields are `missing`. Damped modes (DM) at `a = ±1`
return a `QNMResult` or `QNMEstimate`.
"""
struct QNMEndpointResult{D}
    data::D
end

"""
    QNMPairResult

Compact public wrapper returned by [`qnm_pair`](@ref).
"""
struct QNMPairResult{D}
    data::D
end

for ResultType in (QNMResult, QNMEstimate, QNMFailure, QNMEndpointResult,
        QNMPairResult)
    @eval begin
        function Base.getproperty(result::$ResultType, name::Symbol)
            name === :data && return getfield(result, :data)
            return getproperty(getfield(result, :data), name)
        end

        Base.propertynames(result::$ResultType, private::Bool=false) =
            propertynames(getfield(result, :data), private)
    end
end

function Base.show(io::IO, result::QNMEstimate)
    print(io, "QNMEstimate(mode=", result.mode,
        ", frequency=", result.frequency,
        ", achieved_relative_precision=", result.achieved_relative_precision,
        ", limiting_gate=", result.limiting_gate, ")")
end

function Base.show(io::IO, result::QNMEndpointResult)
    print(io, "QNMEndpointResult(mode=", result.mode,
        ", frequency=", result.frequency,
        ", spectrum_object=", result.spectrum_object, ")")
end

function Base.show(io::IO, ::MIME"text/plain", result::QNMEndpointResult)
    println(io, "ExtremalQuasiNormalMode(")
    println(io, "    parameters = (a = ", result.a,
        ", s = ", result.mode.s, ", l = ", result.mode.l,
        ", m = ", result.mode.m, ", n = ", result.mode.n,
        _qnm_segment_suffix(result), ")")
    println(io, "    omega = ", _qnm_summary_number(result.frequency))
    println(io, "    lambda = ", _qnm_summary_number(result.lambda))
    println(io, "    family = ZDM")
    println(io, "    simple_pole = false")
    print(io, "    formalism = GSN)")
end

function Base.show(io::IO, result::QNMFailure)
    print(io, "QNMFailure(mode=", result.mode,
        ", stage=", result.stage,
        ", stop_reason=", result.stop_reason, ")")
end

function Base.show(io::IO, ::MIME"text/plain", result::QNMFailure)
    println(io, "QNMFailure(")
    println(io, "    parameters = (a = ", result.a,
        ", s = ", result.mode.s, ", l = ", result.mode.l,
        ", m = ", result.mode.m, ", n = ", result.mode.n,
        _qnm_segment_suffix(result), ")")
    println(io, "    stage = ", result.stage)
    if result.frequency isa Number && isfinite(result.frequency)
        println(io, "    omega = ", _qnm_summary_number(result.frequency))
    else
        println(io, "    omega = unavailable")
    end
    println(io, "    status = failed")
    println(io, "    stop_reason = ", result.stop_reason)
    print(io, "    formalism = GSN)")
end

function Base.show(io::IO, result::QNMResult)
    print(io, "QNMResult(mode=", result.mode,
        ", frequency=", result.frequency,
        ", status=", result.status, ")")
end

_qnm_summary_number(value::Real) = @sprintf("%.17g", value)

# Segment 0 adds nothing to the display; later segments are named inline.
function _qnm_segment_suffix(result)
    segment = :segment in propertynames(result) ? result.segment : 0
    return segment == 0 ? "" : string(", segment = ", segment)
end
function _qnm_summary_number(value::Complex)
    real_text = _qnm_summary_number(real(value))
    imaginary = imag(value)
    separator = signbit(imaginary) ? " - " : " + "
    imaginary_text = _qnm_summary_number(abs(imaginary))
    return string(real_text, separator, imaginary_text, "im")
end

# Exact-extremal (a = +-1) damped modes are labelled as such; zero-damping ones are
# QNMEndpointResult.
_qnm_extremal(result) = abs(result.a) == 1

function Base.show(io::IO, ::MIME"text/plain", result::QNMResult)
    println(io, _qnm_extremal(result) ? "ExtremalQuasiNormalMode(" : "QuasiNormalMode(")
    println(io, "    parameters = (a = ", result.a,
        ", s = ", result.mode.s, ", l = ", result.mode.l,
        ", m = ", result.mode.m, ", n = ", result.mode.n,
        _qnm_segment_suffix(result), ")")
    println(io, "    omega = ", _qnm_summary_number(result.frequency))
    println(io, "    lambda = ", _qnm_summary_number(result.lambda))
    println(io, "    incidence_amplitude = ",
        _qnm_summary_number(result.incidence_amplitude))
    println(io, "    reflection_amplitude = ",
        _qnm_summary_number(result.reflection_amplitude))
    println(io, "    incidence_derivative = ",
        _qnm_summary_number(result.incidence_derivative))
    println(io, "    excitation_factor = ",
        _qnm_summary_number(result.excitation_factor))
    _qnm_extremal(result) && println(io, "    family = DM")
    print(io, "    formalism = GSN)")
end

function Base.show(io::IO, ::MIME"text/plain", result::QNMEstimate)
    println(io, _qnm_extremal(result) ?
        "ExtremalQuasiNormalModeEstimate(" : "QuasiNormalModeEstimate(")
    println(io, "    parameters = (a = ", result.a,
        ", s = ", result.mode.s, ", l = ", result.mode.l,
        ", m = ", result.mode.m, ", n = ", result.mode.n,
        _qnm_segment_suffix(result), ")")
    println(io, "    omega = ", _qnm_summary_number(result.frequency))
    println(io, "    lambda = ", _qnm_summary_number(result.lambda))
    println(io, "    incidence_amplitude = ",
        _qnm_summary_number(result.incidence_amplitude))
    println(io, "    reflection_amplitude = ",
        _qnm_summary_number(result.reflection_amplitude))
    println(io, "    incidence_derivative = ",
        _qnm_summary_number(result.incidence_derivative))
    println(io, "    excitation_factor = ",
        _qnm_summary_number(result.excitation_factor))
    println(io, "    achieved_relative_precision = ",
        _qnm_summary_number(result.achieved_relative_precision))
    println(io, "    target_relative_precision = ",
        _qnm_summary_number(result.target_relative_precision))
    if result.inferred_root_shift !== missing
        println(io, "    scaled_incidence_residual = ",
            _qnm_summary_number(result.scaled_incidence_residual))
        println(io, "    incidence_target = ",
            _qnm_summary_number(result.incidence_target))
        println(io, "    inferred_root_shift = ",
            _qnm_summary_number(result.inferred_root_shift))
    end
    println(io, "    limiting_gate = ", result.limiting_gate)
    _qnm_extremal(result) && println(io, "    family = DM")
    println(io, "    status = estimated")
    print(io, "    formalism = GSN)")
end

function Base.show(io::IO, result::QNMPairResult)
    print(io, "QNMPairResult(a=", result.a,
        ", s=", result.s, ", l=", result.l,
        ", m=", result.m, ", n=", result.n,
        ", status=", result.status, ")")
end

function Base.show(io::IO, ::MIME"text/plain", result::QNMPairResult)
    println(io, "QNMPairResult")
    println(io, "  labels:")
    println(io, "    a = ", _qnm_summary_number(result.a))
    println(io, "    s = ", result.s, ", l = ", result.l,
        ", m = ", result.m, ", n = ", result.n)
    println(io, "  ordinary:")
    println(io, "    frequency = ",
        _qnm_summary_number(result.ordinary.frequency))
    println(io, "    status = ", result.ordinary.status)
    println(io, "  mirror:")
    println(io, "    frequency = ",
        _qnm_summary_number(result.mirror.frequency))
    println(io, "    status = ", result.mirror.status)
    println(io, "  status = ", result.status)
    print(io, "  stop_reason = ", result.stop_reason)
end

"""
    LeaverResult

Leaver continued-fraction root for one QNM: frequency, angular separation constant and eigenvalue, residuals and status.
`status` is `:accepted`, `:estimated` or `:failed`, with the reason in `stop_reason`; `overtone_index` is `missing`
for `convention=:leaver`, and `provenance` records the seed, continuation path and gate values.
"""
struct LeaverResult{T,C,A,M,P}
    mode::QNMMode
    convention::Symbol
    overtone_index::Union{Int,Missing}
    inversion_index::Union{Int,Missing}
    a::T
    omega::C
    angular_A::A
    lambda::A
    mixing::M
    cf_value::C
    cf_error::T
    cf_iterations::Int
    root_residual::T
    angular_residual::T
    precision_bits::Int
    status::Symbol
    stop_reason::Symbol
    provenance::P
end

"""
    ISEMValidationResult

Independent ISEM check of a QNM frequency: radial incidence, reflection and transmission amplitudes and the scaled incidence residual.
"""
struct ISEMValidationResult{O,C,T,M}
    omega::O
    incidence::C
    reflection::C
    transmission::C
    absolute_incidence::T
    scaled_incidence::T
    selected_method::Symbol
    status::Symbol
    stop_reason::Symbol
    metadata::M
end

"""
    ExcitationFactorResult

QNM excitation factor with the amplitudes it is built from and their error estimates.
"""
struct ExcitationFactorResult{O,A,AT,R,RT,B,BT,T,M}
    omega::O
    alpha::A
    alpha_teukolsky::AT
    reflection::R
    reflection_teukolsky::RT
    B_gsn::B
    B_teukolsky::BT
    direction_drift::T
    step_drift::T
    bridge_residual::T
    status::Symbol
    stop_reason::Symbol
    metadata::M
end

# Thrown by `leaver_cf_inversion` when the continued fraction does not converge
# within `maximum_iterations` and convergence is required.
struct LentzConvergenceError <: Exception
    error_estimate
    iterations::Int
    maximum_iterations::Int
end

function Base.showerror(io::IO, error::LentzConvergenceError)
    print(io, "Leaver continued fraction did not converge: estimate=",
        error.error_estimate, ", iterations=", error.iterations,
        ", maximum_iterations=", error.maximum_iterations)
end
