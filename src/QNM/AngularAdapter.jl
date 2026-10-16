# Angular eigen-branch at one c = a*omega: separation constant `angular_A`, `lambda`, the spherical-harmonic
# `mixing` vector, residual, truncation order, precision bits, sheet id, status and route-specific `metadata`.
struct AngularBranchResult{T,C,V,M}
    c::C
    angular_A::C
    lambda::C
    mixing::V
    residual::T
    truncation_order::Int
    precision_bits::Int
    sheet_id::Symbol
    status::Symbol
    metadata::M
end

const _QNM_REQUIRED_ANGULAR_APIS = (
    :continue_angular_mode,
)

# The BigFloat angular spectral matrix is `_big_spectral_matrix` in SpinWeightedSpheroidalHarmonics 1.4.0 and
# the precision-generic `construct_spectral_matrix` in later versions; both give the same matrix.
_has_big_spectral_matrix() = isdefined(SpinWeightedSpheroidalHarmonics, :_big_spectral_matrix) ||
    isdefined(SpinWeightedSpheroidalHarmonics, :construct_spectral_matrix)

_big_angular_spectral_matrix(c, s::Int, m::Int, size::Int) =
    isdefined(SpinWeightedSpheroidalHarmonics, :_big_spectral_matrix) ?
        SpinWeightedSpheroidalHarmonics._big_spectral_matrix(c, s, m, size) :
        SpinWeightedSpheroidalHarmonics.construct_spectral_matrix(c, s, m, size)

function _require_qnm_angular_backend()
    missing_apis = Any[filter(
        name -> !isdefined(SpinWeightedSpheroidalHarmonics, name),
        _QNM_REQUIRED_ANGULAR_APIS)...]
    _has_big_spectral_matrix() || push!(missing_apis, :construct_spectral_matrix)
    isempty(missing_apis) && return nothing
    backend_path = try
        pathof(SpinWeightedSpheroidalHarmonics)
    catch
        nothing
    end
    throw(ArgumentError(
        "qnm requires SpinWeightedSpheroidalHarmonics >= 1.4. " *
        "The loaded backend at $(repr(backend_path)) is missing " *
        "$(join(string.(missing_apis), ", ")). Update SpinWeightedSpheroidalHarmonics " *
        "and restart Julia before calling qnm."))
end

"""
    angular_A_to_lambda(angular_A, c, m)

Convert the angular eigenvalue `A` (Leaver's convention, `c = a*omega`) to the separation constant of the radial Teukolsky equation, `lambda = A + c^2 - 2 m c`.
"""
@inline angular_A_to_lambda(angular_A, c, m::Int) =
    angular_A + c^2 - 2m * c

"""
    lambda_to_angular_A(lambda, c, m)

Inverse of `angular_A_to_lambda`: `A = lambda - c^2 + 2 m c`.
"""
@inline lambda_to_angular_A(lambda, c, m::Int) =
    lambda - c^2 + 2m * c

# Float64 angular eigenpair at c, continued on `sheet_id` by SpinWeightedSpheroidalHarmonics.
function _angular_branch_float64(mode::QNMMode, c;
        truncation_order::Int=32, sheet_id::Symbol=:qnm_straight_from_spherical)
    pair = SpinWeightedSpheroidalHarmonics.continue_angular_mode(
        mode.s, mode.l, mode.m, ComplexF64(c);
        sheet_id, truncation_order, cache=nothing)
    return AngularBranchResult(
        pair.c,
        pair.angular_sep,
        pair.lambda,
        pair.coefficients,
        pair.residual,
        pair.truncation_order,
        pair.precision_bits,
        pair.sheet_id,
        pair.status,
        (
            route=:fast_eigenvalue_continuation,
            previous_overlap=pair.previous_overlap,
            overlap_margin=pair.overlap_margin,
            spectral_gap=pair.spectral_gap,
            eigenvector_condition=pair.eigenvector_condition,
            matrix_size=pair.matrix_size,
        ),
    )
end

# BigFloat angular eigenpair: Newton refinement of the Float64 eigenpair on the
# BigFloat spectral matrix.
function _angular_branch_bigfloat(mode::QNMMode, c::Complex{BigFloat};
        truncation_order::Int=32, sheet_id::Symbol=:qnm_straight_from_spherical,
        maximum_iterations::Int=24)
    seed = SpinWeightedSpheroidalHarmonics.continue_angular_mode(
        mode.s, mode.l, mode.m, ComplexF64(c);
        sheet_id, truncation_order, cache=nothing)
    matrix = _big_angular_spectral_matrix(
        c, mode.s, mode.m, seed.matrix_size)
    coefficients = Complex{BigFloat}.(seed.coefficients)
    pivot = argmax(abs.(coefficients))
    coefficients ./= coefficients[pivot]
    value = Complex{BigFloat}(
        BigFloat(real(seed.angular_sep)), BigFloat(imag(seed.angular_sep)))
    target = BigFloat(2)^(-BigFloat(precision(BigFloat) - 16))
    iterations = 0

    for iteration in 1:maximum_iterations
        residual_vector = matrix * coefficients - value * coefficients
        normalization_error = coefficients[pivot] - one(eltype(coefficients))
        residual = max(norm(residual_vector), abs(normalization_error))
        iterations = iteration
        residual <= target && break
        size = seed.matrix_size
        jacobian = zeros(Complex{BigFloat}, size + 1, size + 1)
        jacobian[1:size, 1:size] .= matrix
        for index in 1:size
            jacobian[index, index] -= value
            jacobian[index, size + 1] = -coefficients[index]
        end
        jacobian[size + 1, pivot] = one(eltype(coefficients))
        correction = jacobian \
            (-vcat(residual_vector, normalization_error))
        coefficients .+= view(correction, 1:size)
        value += correction[size + 1]
    end

    coefficients ./= norm(coefficients)
    residual = BigFloat(norm(matrix * coefficients - value * coefficients) /
        max(one(BigFloat), norm(matrix) * norm(coefficients),
            abs(value) * norm(coefficients)))
    lambda = angular_A_to_lambda(value, c, mode.m)
    overlap = abs(dot(Complex{BigFloat}.(seed.coefficients), coefficients)) /
        max(eps(BigFloat), norm(seed.coefficients) * norm(coefficients))
    status = residual <= sqrt(eps(BigFloat)) ?
        :high_precision_refined : :high_precision_stagnated
    return AngularBranchResult(
        c, value, lambda, coefficients, residual, truncation_order,
        precision(BigFloat), sheet_id, status,
        (
            route=:fast_eigenvalue_bigfloat_refinement,
            seed_c=seed.c,
            seed_residual=seed.residual,
            seed_overlap=overlap,
            matrix_size=seed.matrix_size,
            refinement_iterations=iterations,
        ),
    )
end

# Angular separation constant, eigenvalue and coefficients at c = a*omega on
# `sheet_id`, in the working precision of the inputs.
function angular_branch(mode::QNMMode, a, omega;
        truncation_order::Int=32, sheet_id::Symbol=:qnm_straight_from_spherical)
    _require_qnm_angular_backend()
    c = a * omega
    if c isa Complex{BigFloat}
        return _angular_branch_bigfloat(mode, c;
            truncation_order, sheet_id)
    elseif real(c) isa Union{Float16,Float32,Float64}
        return _angular_branch_float64(mode, c;
            truncation_order, sheet_id)
    end
    throw(ArgumentError(
        "QNM angular continuation supports Float16/32/64 and BigFloat inputs."))
end
