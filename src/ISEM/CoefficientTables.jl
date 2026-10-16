module CoefficientTables

using ..Parameters: GSNParameters

export RationalCoefficients, EndpointCoefficientSet
export OrdinaryCoefficientEvaluator, CoefficientSet
export evaluate_ordinary_coefficients!
export evaluate_ordinary_variable!
export register_coefficient_backend!
export gsn_coefficients

# Numerator and denominator polynomial coefficients of one rational ODE coefficient.
struct RationalCoefficients{N<:AbstractVector{ComplexF64},D<:AbstractVector{ComplexF64}}
    numerator::N
    denominator::D
end

# Rational A(x) and B(x) of X'' + A X' + B X = 0 in one endpoint basis (:hc or :ic).
struct EndpointCoefficientSet{A<:RationalCoefficients,B<:RationalCoefficients}
    A::A
    B::B
    basis::Symbol
end

# Generated evaluator that returns the packed P/Q Taylor coefficients of A and B
# about an ordinary point; `eval!` takes x0 and `eval_variable!` takes the
# table's own expansion variable.
struct OrdinaryCoefficientEvaluator{C,F,G}
    basis::Symbol
    variable::Symbol
    coefficients::C
    value_count::Int
    eval!::F
    eval_variable!::G
    requires_latest_world::Bool
end

# Horner evaluation at t of each coefficient block in `ranges`, writing one
# value per block into `out` (zero for an empty block). Used by the generated tables.
@noinline function _evaluate_polynomial_blocks!(out, coefficients, t, ranges)
    @inbounds for j in eachindex(ranges)
        block = ranges[j]
        if isempty(block)
            out[j] = 0.0 + 0.0im
        else
            value = coefficients[last(block)]
            for k in (last(block)-1):-1:first(block)
                value = value*t + coefficients[k]
            end
            out[j] = value
        end
    end
    return out
end

# Fill the first `e.value_count` entries of `out` with the ordinary-point coefficients
# at expansion point `x0`; returns `out`.
function evaluate_ordinary_coefficients!(out::Vector{ComplexF64}, e::OrdinaryCoefficientEvaluator, x0)
    length(out) >= e.value_count || throw(ArgumentError("ordinary output buffer is too short."))
    return e.eval!(out, e.coefficients, x0)
end

# As `evaluate_ordinary_coefficients!`, at the already transformed expansion variable `value`.
function evaluate_ordinary_variable!(out::Vector{ComplexF64}, e::OrdinaryCoefficientEvaluator, value)
    length(out) >= e.value_count || throw(ArgumentError("ordinary output buffer is too short."))
    return e.eval_variable!(out, e.coefficients, value)
end

function (e::OrdinaryCoefficientEvaluator)(x0)
    out = Vector{ComplexF64}(undef, e.value_count)
    evaluate_ordinary_coefficients!(out, e, x0)
    return out
end

# All ODE coefficient data for one mode: horizon and infinity endpoint sets,
# the ordinary-point evaluator, and metadata about how they were built.
struct CoefficientSet{P,H,I,O,M}
    params::P
    horizon::H
    infinity::I
    ordinary::O
    metadata::M
end

const _BACKEND = Ref{Union{Nothing,Function}}(nothing)

# Install the function that builds a CoefficientSet; returns the previous backend.
function register_coefficient_backend!(backend::Function)
    old = _BACKEND[]
    _BACKEND[] = backend
    return old
end

function _missing_backend()
    error(
        "ISEM GSN coefficient tables have not been migrated into the package yet. " *
        "Register a generated backend with register_coefficient_backend! " *
        "or add the stable generated tables to CoefficientTables.jl.",
    )
end

# The CoefficientSet for `params` (or for `s, lambda, m, nu, omega`) from the registered
# backend; `controls` selects the endpoint and ordinary-point bases.
function gsn_coefficients(params::GSNParameters; controls=nothing)
    backend = _BACKEND[]
    backend === nothing && _missing_backend()
    return backend(params; controls=controls)
end

function gsn_coefficients(s::Integer, lambda, m, nu, omega; controls=nothing)
    backend = _BACKEND[]
    backend === nothing && _missing_backend()
    return backend(s, lambda, m, nu, omega; controls=controls)
end

end
