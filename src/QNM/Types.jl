"""
    QNMMode(s, l, m, n, branch=:positive_real)

Label one directly solved QNM branch. `branch=:negative_real` requires an
explicit directly sourced frequency guess; the solver never fills it by
conjugation.
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

Public branch selector for [`qnm`](@ref). `ordinary` selects the positive-real
QNM branch and `mirror` selects the directly solved negative-real branch.
"""
@enum QNMBranch begin
    ordinary
    mirror
end

const Ordinary = ordinary
const Mirror = mirror

"""
    QNMResult

Compact public wrapper returned by [`qnm`](@ref). All detailed fields remain
available through property access, while notebook display omits the internal
radial-route state.
"""
struct QNMResult{D}
    data::D
end

"""
    QNMEstimate

Finite QNM observables that miss one explicitly reported production gate. This
includes a Richardson derivative that has not reached its convergence target,
or a high-damping complex-frequency route whose raw GSN incidence does not
close at an independently identified Leaver root. The result preserves the
computed amplitudes and convergence indicators; it is not an accepted
production result or a rigorous error bound.
"""
struct QNMEstimate{D}
    data::D
end

"""
    QNMFailure

Structured public result returned when a QNM root, independent GSN-ISEM
validation, or excitation-factor gate fails.  It preserves the failing
scientific object and first failed gate without presenting nonfinite amplitudes
as a `QuasiNormalMode` success.
"""
struct QNMFailure{D}
    data::D
end

"""
    QNMEndpointResult

Supported exact-extremal spectral endpoint data. A corotating zero-damping
family ends at a Green-function branch point rather than an isolated simple
pole. This type records that coalesced synchronous endpoint and exposes its
ordinary simple-pole amplitude fields as `missing`. Discrete damped extremal
limits use `QNMResult` after their exact-extremal radial gates pass, or
`QNMEstimate` when only the Richardson production target remains unmet.
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
    println(io, "ExtremalQuasiNormalModeEndpoint(")
    println(io, "    parameters = (a = ", result.a,
        ", s = ", result.mode.s, ", l = ", result.mode.l,
        ", m = ", result.mode.m, ", n = ", result.mode.n,
        _qnm_segment_suffix(result), ")")
    println(io, "    omega = ", _qnm_summary_number(result.frequency))
    println(io, "    lambda = ", _qnm_summary_number(result.lambda))
    println(io, "    spectrum_object = ", result.spectrum_object)
    println(io, "    branch = ", result.mode.branch)
    println(io, "    simple_pole = ", result.simple_pole)
    println(io, "    requested_overtone_label = ",
        result.requested_overtone_label)
    println(io, "    coalesced_overtone_labels = ",
        result.coalesced_overtone_labels)
    println(io, "    overtone_interpretation = ",
        result.overtone_interpretation)
    println(io, "    radial_observables = ", result.radial_observables)
    println(io, "    endpoint_method = ", result.endpoint_method)
    print(io, "    formalism = GSN/Teukolsky endpoint)")
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

# Segment 0 keeps the historical display; later segments are named inline.
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

function Base.show(io::IO, ::MIME"text/plain", result::QNMResult)
    println(io, "QuasiNormalMode(")
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
    print(io, "    formalism = GSN)")
end

function Base.show(io::IO, ::MIME"text/plain", result::QNMEstimate)
    println(io, "QuasiNormalModeEstimate(")
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

Independent GSN-ISEM check of a QNM frequency: radial incidence, reflection and transmission amplitudes and the scaled incidence residual.
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
