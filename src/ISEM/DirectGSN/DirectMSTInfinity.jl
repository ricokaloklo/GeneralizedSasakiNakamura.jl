module DirectMSTInfinity

using LinearAlgebra: cond

using HypergeometricFunctions
using SpecialFunctions: digamma, expint, gamma, loggamma, polygamma
using StaticArrays: MVector, SVector

using ....Transformation: eta_coefficient
using ..DirectCoefficientTables: DirectCoefficientSet
using ..DirectParameters: DirectGSNParameters
using ..DirectIteration:
    DirectBasis,
    DirectLogScaledState,
    direct_endpoint_scale,
    direct_iterate_from_state,
    direct_iterate_pair_from_state,
    direct_iterate_scaled_y,
    direct_logscaled_state,
    direct_materialize_logscaled_state
using ..DirectLFE:
    DDComplex,
    dc_add,
    dc_div,
    dc_imul,
    dc_mul,
    dc_neg,
    dc_sub,
    dc_value
export direct_abel_denominator, direct_mst_infinity_basis, direct_mst_infinity_pair
export direct_mst_plan, direct_mst_eval_plan, direct_mst_state, direct_mst_pin_state
export direct_mst_logscaled_seed, direct_mst_logscaled_state
export direct_mst_scaled_state
export direct_mst_principal_plan, direct_mst_principal_state, direct_mst_anchor
export MSTPhysicalPlan, direct_mst_physical_plan
export direct_mst_physical_state, direct_mst_physical_amplitudes
export direct_mst_monodromy_state
export direct_mst_nia_branch_strength, direct_mst_nia_jump
export mst_nia_amplitudes, mst_principal_amplitudes
export nia_strength_fast
export MSTParams, MSTSeriesData, CFResult
export MSTCertificateError

"""
    MSTCertificateError(kind, message)

Raised when an MST quantity fails its accuracy certificate (series budget, residual,
truncation or conditioning check). `kind` names the failed stage.
"""
struct MSTCertificateError <: Exception
    kind::Symbol
    message::String
end

Base.showerror(io::IO, err::MSTCertificateError) = print(io, err.message)
export alpha_mst, beta_mst, gamma_mst
export rn_cf, ln_cf, mst_guess, mst_nu, mst_nu_complex

const I = 1im
const CF_MAX = 200
const CF_SMALL = 1.0e-50
const MST_N_INIT = 5
const MST_N_MAX = 200
const MST_REL_ERROR = 1.0e-14
const MST_MILLER_MARGIN = 32
const MST_MILLER_STEP = 32
const MST_MILLER_MAX = 512
const MST_MILLER_TOL = MST_REL_ERROR
const MST_NEGATIVE_SEGMENT_TOL = 2.0e-13
const MST_NEGATIVE_SEGMENT_LO = -8
const MST_STATE_REANCHOR_TOL = 1.0e-10
const MST_NEGATIVE_REANCHOR_TOL = 1.0e-13
const MST_NIA_TAIL = 560
const MST_NIA_CHECK_TAIL = 480
const MST_NIA_SUM_MAX = 200
const MST_NIA_CHECK_SUM_MAX = 160
const MST_NIA_CERT_TOL = 1.0e-8
const MST_NIA_FAST_TAIL = 224
const MST_NIA_FAST_CHECK_TAIL = 200
const MST_NIA_FAST_SUM_MAX = 96
const MST_NIA_FAST_CHECK_SUM_MAX = 80
const LOG_FLOAT64_MAX = log(floatmax(Float64))
const LOG_FLOAT64_MIN = log(nextfloat(0.0))
const HGF_2F1 = getproperty(HypergeometricFunctions, Symbol("_\u2082F\u2081"))
const HGF_1F1 = getproperty(HypergeometricFunctions, Symbol("_\u2081F\u2081"))

_finite_complex(z) = isfinite(real(z)) && isfinite(imag(z))

function _isem_module()
    return parentmodule(parentmodule(@__MODULE__))
end

function _root_module()
    return parentmodule(_isem_module())
end

function _matching_module()
    isem = _isem_module()
    isdefined(isem, :Matching) ||
        error("ISEM.Matching is not loaded; direct GSN sfe=true requires the legacy SFE backend.")
    return getfield(isem, :Matching)
end

function _teukolsky_transformation_module()
    matching = _matching_module()
    isdefined(matching, :TeukolskyTransformation) ||
        error("ISEM.Matching.TeukolskyTransformation is not loaded.")
    return getfield(matching, :TeukolskyTransformation)
end

function _conversion_module()
    root = _root_module()
    isdefined(root, :ConversionFactors) ||
        error("GeneralizedSasakiNakamura.ConversionFactors is not loaded.")
    return getfield(root, :ConversionFactors)
end

# Closure r -> matrix converting the GSN function (X, dX/dr*) to the Teukolsky
# radial function.
function _teukolsky_from_gsn_matrix(params)
    solutions = getfield(_root_module(), :Solutions)
    matrix_function = getfield(solutions, :Teukolsky_radial_function_from_Sasaki_Nakamura_function_conversion_matrix)
    return r -> matrix_function(
        params.s,
        params.m,
        params.a,
        params.omega,
        params.lambda,
        r,
    )
end

# r and dx/dr* for the compact coordinate x of the direct GSN solver.
function _direct_x_to_r(params, x)
    rp = 1.0 + params.kappa
    rm = 1.0 - params.kappa
    return (rp - rm * x) / (1.0 - x)
end

function _direct_dx_drstar(params, x)
    kappa = params.kappa
    return kappa * (x - 1.0)^2 * x /
        (kappa + (x - 1.0)^2 + 2 * kappa^2 * x - kappa * x^2)
end

# Factors f3 (in) and f4 (out) between the raw MST series and the unit-normalized
# solutions, used by _raw_unit_logfactor.
function _legacy_factors(params)
    epsilon = 2.0 * params.omega
    kappa = params.kappa
    r_plus = 1.0 + kappa
    f3 = 2.0 * kappa^(2im * params.omega + 1) * exp(1im * params.omega * r_plus)
    f4 = 2.0^(1 + 2 * params.s) * kappa^(2 * params.s - 2im * params.omega + 1) *
        exp(-1im * params.omega * r_plus)
    return ComplexF64(f3), ComplexF64(f4)
end

# Scale from the physical-in MST normalization to the GSN one: GSN transmission
# amplitude Btrans over the transmission amplitude p_trans of P.
function _pin_unit_scale(params)
    epsilon = 2 * params.omega
    tau = (epsilon - params.m * params.a) / params.kappa
    r_plus = 1 + params.kappa
    horizon_p = params.kappa * (epsilon + tau) / (2 * r_plus)
    p_trans = ComplexF64(
        4.0^params.s *
        params.kappa^(2im * horizon_p + 2 * params.s) *
        exp(1im * horizon_p * r_plus),
    )
    conversion = _conversion_module()
    gsn_trans = ComplexF64(getfield(conversion, :Btrans)(
        params.s,
        params.m,
        params.a,
        params.omega,
        params.lambda,
    ))
    scale = gsn_trans / p_trans
    _finite_complex(scale) && !iszero(scale) ||
        error("direct GSN physical-in MST unit scale is nonfinite.")
    return scale
end

"""
    direct_mst_infinity_basis(coefficients, branch, match_x; controls, scratch=nothing, plan=nothing)

Direct GSN infinity solution of `branch` (`:in` or `:out`) at `match_x`: the MST series is
evaluated at the plan's seed point (a new `direct_mst_plan` when `plan` is `nothing`) and
the state is integrated to `match_x` with `direct_iterate_from_state`, or in the scaled
y variable when `controls.basis == :low_frequency_scaled_y`.
"""
function direct_mst_infinity_basis(
    coefficients::DirectCoefficientSet,
    branch::Symbol,
    match_x::Float64;
    controls,
    scratch=nothing,
    plan=nothing,
)
    normalized =
        branch in (:in, :ingoing, :IN, :down, :DOWN) ? :in :
        branch in (:out, :outgoing, :UP, :up) ? :out :
        throw(ArgumentError("MST infinity branch must be :in or :out."))
    kind = normalized == :in ? :infinity_in : :infinity_out
    mst = plan === nothing ?
        direct_mst_plan(coefficients, Float64(match_x); branches=(normalized,)) : plan
    data = mst.data
    seed_x = mst.seed_x
    r = _direct_x_to_r(coefficients.params, seed_x)
    transform = mst.converter(r)
    sequence = normalized == :in ? mst.in_sequence : mst.out_sequence
    budget = normalized == :in ? mst.in_budget : mst.out_budget
    factor = normalized == :in ? mst.in_factor : mst.out_factor
    seed_value, seed_derivative, _ =
        _mst_state_at(coefficients, data, normalized, seed_x, transform,
            sequence, budget, factor)
    if controls.basis == :low_frequency_scaled_y && seed_x != match_x
        return direct_iterate_scaled_y(coefficients, kind, seed_x,
            seed_value, seed_derivative;
            target_y=(1.0-match_x)/abs(coefficients.params.omega),
            order=controls.ordinary_order, tolerance=controls.tolerance)
    end
    return direct_iterate_from_state(
        coefficients,
        kind,
        seed_x,
        seed_value,
        seed_derivative,
        match_x;
        controls=controls,
        scratch=scratch,
    )
end

"""
    direct_mst_infinity_pair(coefficients, branch1, branch2, match_x; controls, scratch=nothing, plan=nothing)

Two infinity solutions as in `direct_mst_infinity_basis`, sharing one MST plan and
seed point.
"""
function direct_mst_infinity_pair(
    coefficients::DirectCoefficientSet,
    branch1::Symbol,
    branch2::Symbol,
    match_x::Float64;
    controls,
    scratch=nothing,
    plan=nothing,
)
    normalized1 =
        branch1 in (:in, :ingoing, :IN, :down, :DOWN) ? :in :
        branch1 in (:out, :outgoing, :UP, :up) ? :out :
        throw(ArgumentError("MST infinity branch must be :in or :out."))
    normalized2 =
        branch2 in (:in, :ingoing, :IN, :down, :DOWN) ? :in :
        branch2 in (:out, :outgoing, :UP, :up) ? :out :
        throw(ArgumentError("MST infinity branch must be :in or :out."))
    kind1 = normalized1 == :in ? :infinity_in : :infinity_out
    kind2 = normalized2 == :in ? :infinity_in : :infinity_out
    mst = plan === nothing ? direct_mst_plan(
        coefficients,
        Float64(match_x);
        branches=(normalized1, normalized2),
    ) : plan
    data = mst.data
    seed_x = mst.seed_x
    r = _direct_x_to_r(coefficients.params, seed_x)
    transform = mst.converter(r)
    sequence1 = normalized1 == :in ? mst.in_sequence : mst.out_sequence
    sequence2 = normalized2 == :in ? mst.in_sequence : mst.out_sequence
    budget1 = normalized1 == :in ? mst.in_budget : mst.out_budget
    budget2 = normalized2 == :in ? mst.in_budget : mst.out_budget
    factor1 = normalized1 == :in ? mst.in_factor : mst.out_factor
    factor2 = normalized2 == :in ? mst.in_factor : mst.out_factor
    seed_value1, seed_derivative1, _ =
        _mst_state_at(coefficients, data, normalized1, seed_x, transform,
            sequence1, budget1, factor1)
    seed_value2, seed_derivative2, _ =
        _mst_state_at(coefficients, data, normalized2, seed_x, transform,
            sequence2, budget2, factor2)
    if controls.basis == :low_frequency_scaled_y && seed_x != match_x
        target_y = (1.0-match_x)/abs(coefficients.params.omega)
        first = direct_iterate_scaled_y(coefficients, kind1, seed_x,
            seed_value1, seed_derivative1; target_y,
            order=controls.ordinary_order, tolerance=controls.tolerance)
        second = direct_iterate_scaled_y(coefficients, kind2, seed_x,
            seed_value2, seed_derivative2; target_y,
            order=controls.ordinary_order, tolerance=controls.tolerance)
        return first, second
    end
    return direct_iterate_pair_from_state(
        coefficients,
        kind1,
        kind2,
        seed_x,
        seed_value1,
        seed_derivative1,
        seed_value2,
        seed_derivative2,
        match_x;
        controls=controls,
        scratch=scratch,
    )
end

"""
    MSTParams(s, l, m, a, omega, lambda, nu[, nu_offset])

Parameters of an MST series for |a| < 1: the mode, `epsilon = 2 omega`,
`kappa = sqrt(1 - a^2)`, `tau = (epsilon - m a) / kappa`, `nu` and `nu_offset = nu - l`
(kept separately so that its digits survive when nu is close to l).
"""
struct MSTParams
    s::Int
    l::Int
    m::Int
    a::Float64
    omega::ComplexF64
    lambda::ComplexF64
    epsilon::ComplexF64
    kappa::Float64
    tau::ComplexF64
    nu::ComplexF64
    nu_offset::ComplexF64
end

function MSTParams(s::Integer, l::Integer, m::Integer, a::Real, omega::Number, lambda, nu)
    nu_value = ComplexF64(nu)
    return MSTParams(s, l, m, a, omega, lambda, nu_value,
        nu_value - ComplexF64(l))
end

function MSTParams(
    s::Integer,
    l::Integer,
    m::Integer,
    a::Real,
    omega::Number,
    lambda,
    nu,
    nu_offset,
)
    af = Float64(a)
    omegac = ComplexF64(omega)
    kappa2 = 1.0 - af * af
    kappa2 > 0.0 || error("MSTParams requires non-extremal |a| < 1")
    epsilon = 2.0 * omegac
    kappa = sqrt(kappa2)
    tau = (epsilon - Int(m) * af) / kappa
    return MSTParams(Int(s), Int(l), Int(m), af, omegac, ComplexF64(lambda),
        epsilon, kappa, tau, ComplexF64(nu), ComplexF64(nu_offset))
end

@inline _nu_shift(p::MSTParams, shift::Integer=0) =
    ComplexF64(p.l + shift) + p.nu_offset

@inline _npnu(p::MSTParams, n::Integer) =
    ComplexF64(p.l + n) + p.nu_offset

@inline _nu_factor(p::MSTParams, n::Integer, shift::Integer=0) =
    ComplexF64(p.l + n + shift) + p.nu_offset

@inline _twice_nu_factor(p::MSTParams, n::Integer, shift::Integer=0) =
    ComplexF64(2 * (p.l + n) + shift) + 2 * p.nu_offset

"""
    alpha_mst(n, p), beta_mst(n, p), gamma_mst(n, p)

Coefficients of the MST three-term recurrence
alpha_n a_(n+1) + beta_n a_n + gamma_n a_(n-1) = 0
(Sasaki & Tagoshi, Living Rev. Relativ. 6, 6 (2003)).
`beta_mst` switches to double-double when its terms cancel.
"""
function alpha_mst(n::Integer, p::MSTParams)
    eps = ComplexF64(p.epsilon)
    kappa = ComplexF64(p.kappa)
    tau = ComplexF64(p.tau)
    npnu1 = _nu_factor(p, n, 1)
    spin_factor = _nu_factor(p, n, 1 + p.s)
    return I * eps * kappa *
        (spin_factor + I * eps) *
        (spin_factor - I * eps) *
        (npnu1 + I * tau) /
        (npnu1 * _twice_nu_factor(p, n, 3))
end

"""
    beta_mst(n, p)

Middle coefficient of the MST recurrence; see `alpha_mst`.
"""
function beta_mst(n::Integer, p::MSTParams)
    s = ComplexF64(p.s)
    eps = ComplexF64(p.epsilon)
    kappa = ComplexF64(p.kappa)
    tau = ComplexF64(p.tau)
    product = _nu_factor(p, n) * _nu_factor(p, n, 1)
    angular = -p.lambda - s * (s + 1) + product
    eps2 = eps^2
    linear = eps * kappa * tau
    correction = linear * (s^2 + eps2) / product
    value = angular + eps2 + linear + correction
    scale = abs(p.lambda) + abs(s*(s+1)) + abs(product) +
        abs(eps2) + abs(linear) + abs(correction)
    # Both the angular terms and the frequency terms can cancel; under severe
    # cancellation beta_n is recomputed in double-double (_dd_beta).
    if Base.eps(Float64)*scale > MST_COEFF_ERROR_MAX*abs(value)
        return dc_value(_dd_beta(Int(n), p, DDComplex(p.nu_offset)))
    end
    return value
end

"""
    gamma_mst(n, p)

Lower coefficient of the MST recurrence; see `alpha_mst`.
"""
function gamma_mst(n::Integer, p::MSTParams)
    eps = ComplexF64(p.epsilon)
    kappa = ComplexF64(p.kappa)
    tau = ComplexF64(p.tau)
    npnu = _nu_factor(p, n)
    spin_factor = _nu_factor(p, n, -p.s)
    return -I * eps * kappa *
        (spin_factor + I * eps) *
        (spin_factor - I * eps) *
        (npnu - I * tau) /
        (npnu * _twice_nu_factor(p, n, -1))
end

"""
    CFResult

Continued-fraction result: `value`, last relative change `relerr`, `iterations` and
`converged`.
"""
struct CFResult
    value::ComplexF64
    relerr::Float64
    iterations::Int
    converged::Bool
end

_cf_safe(v::ComplexF64) = iszero(v) ? ComplexF64(CF_SMALL) : v

# Modified Lentz evaluation of a_0 / (b_0 + a_1 / (b_1 + ...)) for coefficient functions of
# j = 0, 1, ...; a zero denominator is replaced by CF_SMALL. Returns a CFResult.
function lentz_cf(a_coeff, b_coeff; tol::Float64=eps(Float64), maxiter::Int=CF_MAX)
    c_prev = ComplexF64(CF_SMALL)
    d_prev = ComplexF64(0)
    value_prev = ComplexF64(CF_SMALL)
    relerr = Inf
    j = -1
    for iter in 1:maxiter
        j += 1
        an = ComplexF64(a_coeff(j))
        bn = ComplexF64(b_coeff(j))
        d_next = bn + an * d_prev
        c_next = bn + an / _cf_safe(c_prev)

        modified = false
        if iszero(c_next) || (!iszero(an) && abs(1 + bn * c_prev / an) <= 1.0e-14)
            c_next = ComplexF64(CF_SMALL)
            modified = true
        end
        if iszero(d_next) || (!iszero(an) && abs(1 + bn / (an * _cf_safe(d_prev))) <= 1.0e-14)
            d_next = ComplexF64(1 / CF_SMALL)
            modified = true
        else
            d_next = inv(d_next)
        end

        value_next = value_prev * d_next * c_next
        relerr = abs(1 - value_prev / value_next)
        c_prev = c_next
        d_prev = d_next
        value_prev = value_next
        if relerr <= tol && !modified
            return CFResult(value_next, relerr, iter, true)
        end
    end
    return CFResult(value_prev, relerr, maxiter, false)
end

"""
    rn_cf(n, p)

Continued fraction for the MST coefficient ratio R_n = a_n / a_(n-1), as a `CFResult`.
"""
function rn_cf(n::Integer, p::MSTParams)
    numerator = -alpha_mst(n - 1, p) * gamma_mst(n, p)
    normalize = !iszero(numerator) &&
        CF_SMALL*abs(beta_mst(n,p)) > MST_COEFF_ERROR_MAX*abs(numerator)
    a_coeff(j) = normalize && j == 0 ? one(numerator) :
        -alpha_mst(j + n - 1, p) * gamma_mst(j + n, p)
    b_coeff(j) = beta_mst(j + n, p)
    cf = lentz_cf(a_coeff, b_coeff)
    value = normalize ? -gamma_mst(n,p)*cf.value : cf.value / alpha_mst(n-1,p)
    return CFResult(value, cf.relerr, cf.iterations, cf.converged)
end

# L_n = a_n / a_(n+1) from the continued fraction alone.
function _ln_cf_unpaired(n::Integer, p::MSTParams)
    numerator = -alpha_mst(n, p) * gamma_mst(n + 1, p)
    # Keep a small, nonzero leading numerator outside Lentz's zero guard.
    normalize = !iszero(numerator) &&
        CF_SMALL*abs(beta_mst(n,p)) > MST_COEFF_ERROR_MAX*abs(numerator)
    a_coeff(j) = normalize && j == 0 ? one(numerator) :
        -alpha_mst(-j + n, p) * gamma_mst(-j + n + 1, p)
    b_coeff(j) = beta_mst(-j + n, p)
    cf = lentz_cf(a_coeff, b_coeff)
    value = normalize ? -alpha_mst(n,p)*cf.value : cf.value / gamma_mst(n+1,p)
    return CFResult(value, cf.relerr, cf.iterations, cf.converged)
end

"""
    ln_cf(n, p)

Continued fraction for the MST coefficient ratio L_n = a_n / a_(n+1), as a `CFResult`.
For -l <= n < 0 the rows at nu + n = 0 and -1 are eliminated together when their
cancellation would otherwise lose accuracy.
"""
function ln_cf(n::Integer, p::MSTParams)
    if -p.l <= n < 0
        # Eliminate the adjacent nu+n=0,-1 rows together. Their 1/(nu-l)^2
        # terms cancel algebraically, before floating-point evaluation.
        j = -p.l
        t = p.nu_offset
        q = t*t
        e = p.epsilon
        ek = e*p.kappa
        tau2 = p.tau*p.tau
        spin2 = p.s*p.s
        A = spin2 + e*e
        L = p.lambda + p.s*(p.s+1) - e*e - ek*p.tau
        D = ek*p.tau*A
        polynomial = ((-q - 2A - tau2 + 4spin2 + 1)*q +
            (-A*A - 2A*tau2 + 2A + 4tau2*spin2 + tau2 - 4spin2))*q +
            (3A*A*tau2 + A*A + 2A*tau2 - 4tau2*spin2)
        d1 = (q-L)^2 - q
        d2 = 2D*(q+1-L)/(q-1)
        d3 = ek*ek*polynomial/((q-1)*(4q-1))
        determinant = d1+d2+d3
        roundoff = eps(Float64)*(abs(beta_mst(j,p)*beta_mst(j-1,p)) +
            abs(gamma_mst(j,p)*alpha_mst(j-1,p)))/abs(determinant)
        roundoff > MST_COEFF_ERROR_MAX || return _ln_cf_unpaired(n,p)
        tail = _ln_cf_unpaired(j-2,p)
        coupled = gamma_mst(j-1,p)*tail.value
        W = beta_mst(j-1,p) + coupled
        d4 = beta_mst(j,p)*coupled
        denominator = d1+d2+d3+d4
        value = -alpha_mst(j,p)*W/denominator
        error = (eps(Float64)*(abs(d1)+abs(d2)+abs(d3)+abs(d4)) +
            abs(d4)*tail.relerr)/abs(denominator) +
            (eps(Float64)*(abs(W-coupled)+abs(coupled)) +
                abs(coupled)*tail.relerr)/abs(W)
        for k in (j+1):n
            beta = beta_mst(k,p)
            coupled = gamma_mst(k,p)*value
            denominator = beta+coupled
            error = (eps(Float64)*(abs(beta)+abs(coupled)) +
                abs(coupled)*error)/abs(denominator) + eps(Float64)
            value = -alpha_mst(k,p)/denominator
        end
        return CFResult(value,error,tail.iterations+n-j+2,tail.converged)
    end
    return _ln_cf_unpaired(n,p)
end

"""
    MSTSeriesData(params)

Lazily filled MST series coefficients: `log a_n` with a_0 = 1, the relative error of each
ratio, its convergence flag and the accumulated coefficient errors, guarded by a lock.
"""
mutable struct MSTSeriesData
    params::MSTParams
    log_coeffs::Dict{Int,ComplexF64}
    cf_relerrs::Dict{Int,Float64}
    cf_ok::Dict{Int,Bool}
    coeff_errors::Dict{Int,Float64}
    coefficient_lock::ReentrantLock
end

MSTSeriesData(params::MSTParams) =
    MSTSeriesData(
        params,
        Dict(0 => ComplexF64(0)),
        Dict{Int,Float64}(),
        Dict{Int,Bool}(),
        Dict(0 => 0.0),
        ReentrantLock(),
    )

const MST_COEFF_CONDITION_MAX = 4.0
const MST_COEFF_ERROR_MAX = 2.5e-14

# Next ratio a_key / a_(key-1) (_next_rn) or a_key / a_(key+1) (_next_ln) from the
# recurrence and the two previous coefficients, with its error; `valid` is false when the
# step is ill conditioned and the continued fraction is used instead.
function _next_rn(data::MSTSeriesData, key::Int)
    n = key - 1
    previous = exp(data.log_coeffs[n] - data.log_coeffs[n - 1])
    beta = beta_mst(n, data.params)
    coupled = gamma_mst(n, data.params) / previous
    numerator = beta + coupled
    denominator = alpha_mst(n, data.params)
    ratio = -numerator / denominator
    condition = (abs(beta) + abs(coupled)) /
        max(abs(numerator), floatmin(Float64))
    error = condition * eps(Float64) +
        abs(coupled) * data.cf_relerrs[n] /
        max(abs(numerator), floatmin(Float64))
    valid = _finite_complex(ratio) && !iszero(ratio) &&
        isfinite(condition) && condition <= MST_COEFF_CONDITION_MAX &&
        isfinite(error) && error <= MST_COEFF_ERROR_MAX
    return ratio, Float64(error + eps(Float64)), valid
end

function _next_ln(data::MSTSeriesData, key::Int)
    n = key + 1
    previous = exp(data.log_coeffs[n] - data.log_coeffs[n + 1])
    beta = beta_mst(n, data.params)
    coupled = alpha_mst(n, data.params) / previous
    numerator = beta + coupled
    denominator = gamma_mst(n, data.params)
    ratio = -numerator / denominator
    condition = (abs(beta) + abs(coupled)) /
        max(abs(numerator), floatmin(Float64))
    error = condition * eps(Float64) +
        abs(coupled) * data.cf_relerrs[n] /
        max(abs(numerator), floatmin(Float64))
    valid = _finite_complex(ratio) && !iszero(ratio) &&
        isfinite(condition) && condition <= MST_COEFF_CONDITION_MAX &&
        isfinite(error) && error <= MST_COEFF_ERROR_MAX
    return ratio, Float64(error + eps(Float64)), valid
end

# log a_n, filled on demand outward from n = 0 by recurrence when well conditioned,
# otherwise by continued fraction.
function log_series_coefficient!(data::MSTSeriesData, n::Integer)
    return lock(data.coefficient_lock) do
        _log_series_coefficient!(data, Int(n))
    end
end

function _log_series_coefficient!(data::MSTSeriesData, key::Int)
    if haskey(data.log_coeffs, key)
        return data.log_coeffs[key]
    elseif key > 0
        prev = _log_series_coefficient!(data, key - 1)
        ratio, ratio_error, recurrence_ok = key > 1 ?
            _next_rn(data, key) : (ComplexF64(NaN), Inf, false)
        if recurrence_ok
            data.log_coeffs[key] = prev + log(ratio)
            data.cf_relerrs[key] = ratio_error
            data.cf_ok[key] = true
            step_error = ratio_error
        else
            cf = rn_cf(key, data.params)
            data.log_coeffs[key] = prev + log(cf.value)
            data.cf_relerrs[key] = cf.relerr
            data.cf_ok[key] = cf.converged
            step_error = cf.converged ? cf.relerr : Inf
        end
        data.coeff_errors[key] = min(Inf,
            data.coeff_errors[key - 1] + step_error + eps(Float64))
        return data.log_coeffs[key]
    else
        prev = _log_series_coefficient!(data, key + 1)
        ratio, ratio_error, recurrence_ok = key < -1 ?
            _next_ln(data, key) : (ComplexF64(NaN), Inf, false)
        if recurrence_ok
            data.log_coeffs[key] = prev + log(ratio)
            data.cf_relerrs[key] = ratio_error
            data.cf_ok[key] = true
            step_error = ratio_error
        else
            cf = ln_cf(key, data.params)
            data.log_coeffs[key] = prev + log(cf.value)
            data.cf_relerrs[key] = cf.relerr
            data.cf_ok[key] = cf.converged
            step_error = cf.converged ? cf.relerr : Inf
        end
        data.coeff_errors[key] = min(Inf,
            data.coeff_errors[key + 1] + step_error + eps(Float64))
        return data.log_coeffs[key]
    end
end

# MST coordinates x = (r_+ - r) / (r_+ - r_-) and z = epsilon kappa (1 - x).
mst_x(p::MSTParams, r::Number) =
    (1.0 + p.kappa - r) / (2.0 * p.kappa)
mst_z(p::MSTParams, r::Number) = p.epsilon * p.kappa * (1.0 - mst_x(p, r))

struct MSTSum
    value::ComplexF64
    derivative::ComplexF64
    estimated_relerr::Float64
    nmin::Int
    nmax::Int
    max_term_abs::Float64
    log_scale::Float64
end

struct MSTTriplet
    value::ComplexF64
    deriv::ComplexF64
    second::ComplexF64
end

struct LogNormSum
    phase::ComplexF64
    logabs::Float64
    estimated_relerr::Float64
    nmin::Int
    nmax::Int
end

function log_pochhammer(a, n::Integer)
    ac = ComplexF64(a)
    n == 0 && return ComplexF64(0)
    total = ComplexF64(0)
    if n > 0
        for k in 0:(n - 1)
            total += log(ac + k)
        end
        return total
    end
    for k in 1:(-n)
        total += log(ac - k)
    end
    return -total
end

# Tolerances and limits of the confluent hypergeometric U(a, b, z) evaluators.
struct HyperUBackend
    reltol::Float64
    series_max::Int
    asymptotic_max::Int
    asymptotic_min_abs_z::Float64
    u11_limit_tol::Float64
    integer_b_tol::Float64
    integer_b_priority_max_abs_z::Float64
    near_integer_max_delta::Float64
    near_integer_max_abs_z::Float64
end

const HYPER_U = HyperUBackend(
    1.0e-14,
    800,
    120,
    5.0,
    1.0e-6,
    1.0e-6,
    12.0,
    1.0e-3,
    12.0,
)

struct HyperEval
    value::ComplexF64
    route::Symbol
    estimated_relerr::Float64
    iterations::Int
    status::Symbol
    message::String
end

function _near_zero(z; atol=1.0e-14, rtol=1.0e-14)
    return abs(ComplexF64(z)) <= atol + rtol * max(1.0, abs(z))
end

function _heval(value, route, relerr, iterations; status=:OK, message="")
    z = ComplexF64(value)
    st = _finite_complex(z) ? status : :NONFINITE
    return HyperEval(z, route, Float64(relerr), Int(iterations), st, String(message))
end

_hfail(route, message) =
    HyperEval(ComplexF64(NaN + NaN * im), route, Inf, 0, :FAILED, String(message))

_usable(status::Symbol) =
    status in (:OK, :WARN_LIMIT_APPROX, :WARN_CONDITIONED, :WARN_ASYMPTOTIC_LEASTTERM, :WARN_MAXITER)

function _hrelerr(term, total)
    return abs(term) / max(abs(total), floatmin(Float64))
end

# Truncated Taylor series ("jets") in delta, stored as SVector coefficients.
const U_JET_DEGREE = 6
const U_JET_SIZE = U_JET_DEGREE + 1
const UJet = SVector{U_JET_SIZE,ComplexF64}
const U_SMALL_JET_SIZE = 4
const U_SMALL_DELTA = 1.0e-8

@inline _jet_zero() = zero(UJet)

@inline _jet_zero(::Val{N}) where {N} = zero(SVector{N,ComplexF64})

@inline function _jet_const(value)
    return Base.setindex(_jet_zero(), ComplexF64(value), 1)
end

@inline function _jet_const(value, size::Val{N}) where {N}
    return Base.setindex(_jet_zero(size), ComplexF64(value), 1)
end

@inline function _jet_linear(value, slope)
    jet = Base.setindex(_jet_zero(), ComplexF64(value), 1)
    return Base.setindex(jet, ComplexF64(slope), 2)
end

@inline function _jet_linear(value, slope, size::Val{N}) where {N}
    jet = Base.setindex(_jet_zero(size), ComplexF64(value), 1)
    return Base.setindex(jet, ComplexF64(slope), 2)
end

@inline function _jet_mul(
    a::SVector{N,ComplexF64},
    b::SVector{N,ComplexF64},
) where {N}
    out = MVector{N,ComplexF64}(undef)
    @inbounds for n in 0:(N - 1)
        value = ComplexF64(0)
        for k in 0:n
            value += a[k + 1] * b[n - k + 1]
        end
        out[n + 1] = value
    end
    return SVector{N,ComplexF64}(out)
end

@inline function _jet_inv(a::SVector{N,ComplexF64}) where {N}
    iszero(a[1]) && return SVector{N,ComplexF64}(
        ntuple(_ -> ComplexF64(NaN + NaN * im), N),
    )
    out = MVector{N,ComplexF64}(undef)
    out[1] = inv(a[1])
    @inbounds for n in 1:(N - 1)
        value = ComplexF64(0)
        for k in 1:n
            value += a[k + 1] * out[n - k + 1]
        end
        out[n + 1] = -value / a[1]
    end
    return SVector{N,ComplexF64}(out)
end

@inline _jet_div(
    a::SVector{N,ComplexF64},
    b::SVector{N,ComplexF64},
) where {N} = _jet_mul(a, _jet_inv(b))

@inline function _jet_exp(a::SVector{N,ComplexF64}) where {N}
    out = MVector{N,ComplexF64}(undef)
    out[1] = exp(a[1])
    @inbounds for n in 1:(N - 1)
        value = ComplexF64(0)
        for k in 1:n
            value += k * a[k + 1] * out[n - k + 1]
        end
        out[n + 1] = value / n
    end
    return SVector{N,ComplexF64}(out)
end

@inline function _jet_eval(a::SVector{N,ComplexF64}, delta) where {N}
    value = a[end]
    @inbounds for n in (N - 1):-1:1
        value = muladd(value, delta, a[n])
    end
    return value
end

@inline function _jet_norm(a::SVector{N,ComplexF64}) where {N}
    value = 0.0
    @inbounds for coefficient in a
        value = max(value, abs(coefficient))
    end
    return value
end

# Jet of 1/Gamma(z0 + slope delta); the argument is shifted to Re >= 2 before loggamma
# and polygamma are used.
_rgamma_jet(z0, slope) = _rgamma_jet(z0, slope, Val(U_JET_SIZE))

function _rgamma_jet(z0, slope, size::Val{N}) where {N}
    zc = ComplexF64(z0)
    qc = ComplexF64(slope)
    shift = max(0, ceil(Int, 2.0 - real(zc)))
    base = zc + shift
    h = MVector{N,ComplexF64}(undef)
    h[1] = -loggamma(base)
    @inbounds for n in 1:(N - 1)
        h[n + 1] = -(qc^n) * polygamma(n - 1, base) / factorial(n)
    end
    result = _jet_exp(SVector{N,ComplexF64}(h))
    @inbounds for j in 0:(shift - 1)
        result = _jet_mul(result, _jet_linear(zc + j, qc, size))
    end
    return result
end

# Jets of 1/Gamma(beta0 + k + slope delta) for k = 0..max(0, 1 - beta0); cached below for
# |beta0| <= 256 and slope +-1.
function _gamma_prefix(beta0::Int, beta_slope, size::Val{N}) where {N}
    k0 = max(0, 1 - beta0)
    prefix = Vector{SVector{N,ComplexF64}}(undef, k0 + 1)
    prefix[k0 + 1] = _rgamma_jet(beta0 + k0, beta_slope, size)
    @inbounds for k in (k0 - 1):-1:0
        prefix[k + 1] = _jet_mul(
            _jet_linear(beta0 + k, beta_slope, size),
            prefix[k + 2],
        )
    end
    return prefix
end

const _GAMMA_PREFIX7 = let
    cache = Dict{Tuple{Int,Int},Vector{SVector{7,ComplexF64}}}()
    for beta0 in -256:256, slope in (-1, 1)
        cache[(beta0, slope)] = _gamma_prefix(beta0, slope, Val(7))
    end
    cache
end

const _GAMMA_PREFIX4 = let
    cache = Dict{Tuple{Int,Int},Vector{SVector{4,ComplexF64}}}()
    for beta0 in -256:256, slope in (-1, 1)
        cache[(beta0, slope)] = _gamma_prefix(beta0, slope, Val(4))
    end
    cache
end

function _cached_prefix(beta0::Int, beta_slope, size::Val{N}) where {N}
    if beta_slope isa Integer && beta_slope in (-1, 1)
        cache = N == 7 ? _GAMMA_PREFIX7 : N == 4 ? _GAMMA_PREFIX4 : nothing
        if cache !== nothing
            cached = get(cache, (beta0, Int(beta_slope)), nothing)
            cached === nothing || return cached
        end
    end
    return _gamma_prefix(beta0, beta_slope, size)
end

# Jet in delta of the regularized Kummer series M(a, b, z) / Gamma(b) with
# a = alpha0 + alpha_slope delta and b = beta0 + beta_slope delta (compensated summation),
# with its error estimate and condition. _mreg_triplet adds the first two z-derivatives.
function _mreg_jet(alpha0, alpha_slope, beta0::Int, beta_slope, z, backend::HyperUBackend)
    k0 = max(0, 1 - beta0)
    k0 + 4 < backend.series_max ||
        return (jet=_jet_zero(), relerr=Inf, iterations=0,
            status=:FAILED, condition=Inf)

    gamma_prefix = _cached_prefix(beta0, beta_slope, Val(U_JET_SIZE))

    total = _jet_zero()
    compensation = _jet_zero()
    pochhammer = _jet_const(1)
    gamma_term = gamma_prefix[1]
    zterm = ComplexF64(1)
    maxterm = 0.0
    last_rel = Inf
    iterations = 0
    converged = false

    @inbounds for k in 0:(backend.series_max - 1)
        if k <= k0
            gamma_term = gamma_prefix[k + 1]
        else
            gamma_term = _jet_div(
                gamma_term,
                _jet_linear(beta0 + k - 1, beta_slope),
            )
        end
        term = zterm * _jet_mul(pochhammer, gamma_term)
        if !all(_finite_complex, term)
            return (jet=total, relerr=Inf, iterations=k,
                status=:NONFINITE, condition=Inf)
        end
        corrected = term - compensation
        updated = total + corrected
        compensation = (updated - total) - corrected
        total = updated
        term_norm = _jet_norm(term)
        total_norm = _jet_norm(total)
        maxterm = max(maxterm, term_norm)
        last_rel = term_norm / max(total_norm, floatmin(Float64))
        iterations = k + 1
        if k >= max(k0 + 3, U_JET_DEGREE + 2) && last_rel <= backend.reltol
            converged = true
            break
        end
        pochhammer = _jet_mul(
            pochhammer,
            _jet_linear(ComplexF64(alpha0) + k, alpha_slope),
        )
        zterm *= ComplexF64(z) / (k + 1)
    end

    total_norm = _jet_norm(total)
    condition = maxterm / max(total_norm, floatmin(Float64))
    roundoff = eps(Float64) * condition
    status = converged ? :OK : :WARN_MAXITER
    return (jet=total, relerr=max(last_rel, roundoff), iterations=iterations,
        status=status, condition=condition)
end

function _mreg_triplet(alpha0, alpha_slope, beta0::Int, beta_slope, z,
        backend::HyperUBackend)
    return _mreg_triplet(alpha0, alpha_slope, beta0, beta_slope, z,
        backend, Val(U_JET_SIZE))
end

function _mreg_triplet(alpha0, alpha_slope, beta0::Int, beta_slope, z,
        backend::HyperUBackend, size::Val{N}) where {N}
    degree = N - 1
    k0 = max(0, 1 - beta0)
    k0 + 4 < backend.series_max ||
        return (jet=_jet_zero(size), dz=_jet_zero(size), dzz=_jet_zero(size),
            relerr=Inf, iterations=0, status=:FAILED, condition=Inf)

    gamma_prefix = _cached_prefix(beta0, beta_slope, size)

    total = _jet_zero(size)
    compensation = _jet_zero(size)
    total_dz = _jet_zero(size)
    compensation_dz = _jet_zero(size)
    total_dzz = _jet_zero(size)
    compensation_dzz = _jet_zero(size)
    pochhammer = _jet_const(1, size)
    gamma_term = gamma_prefix[1]
    zterm = ComplexF64(1)
    maxterm = 0.0
    last_rel = Inf
    iterations = 0
    converged = false

    @inbounds for k in 0:(backend.series_max - 1)
        if k <= k0
            gamma_term = gamma_prefix[k + 1]
        else
            gamma_term = _jet_div(
                gamma_term,
                _jet_linear(beta0 + k - 1, beta_slope, size),
            )
        end
        term = zterm * _jet_mul(pochhammer, gamma_term)
        if !all(_finite_complex, term)
            return (jet=total, dz=total_dz, dzz=total_dzz,
                relerr=Inf, iterations=k, status=:NONFINITE, condition=Inf)
        end
        term_dz = k == 0 ? _jet_zero(size) : (k / ComplexF64(z)) * term
        term_dzz = k <= 1 ? _jet_zero(size) :
            (k * (k - 1) / ComplexF64(z)^2) * term
        corrected = term - compensation
        updated = total + corrected
        compensation = (updated - total) - corrected
        total = updated
        corrected_dz = term_dz - compensation_dz
        updated_dz = total_dz + corrected_dz
        compensation_dz = (updated_dz - total_dz) - corrected_dz
        total_dz = updated_dz
        corrected_dzz = term_dzz - compensation_dzz
        updated_dzz = total_dzz + corrected_dzz
        compensation_dzz = (updated_dzz - total_dzz) - corrected_dzz
        total_dzz = updated_dzz
        term_norm = max(_jet_norm(term), _jet_norm(term_dz), _jet_norm(term_dzz))
        total_norm = max(_jet_norm(total), _jet_norm(total_dz), _jet_norm(total_dzz))
        maxterm = max(maxterm, term_norm)
        last_rel = term_norm / max(total_norm, floatmin(Float64))
        iterations = k + 1
        if k >= max(k0 + 3, degree + 2) && last_rel <= backend.reltol
            converged = true
            break
        end
        pochhammer = _jet_mul(
            pochhammer,
            _jet_linear(ComplexF64(alpha0) + k, alpha_slope, size),
        )
        zterm *= ComplexF64(z) / (k + 1)
    end

    total_norm = max(_jet_norm(total), _jet_norm(total_dz), _jet_norm(total_dzz))
    condition = maxterm / max(total_norm, floatmin(Float64))
    roundoff = eps(Float64) * condition
    status = converged ? :OK : :WARN_MAXITER
    return (jet=total, dz=total_dz, dzz=total_dzz,
        relerr=max(last_rel, roundoff), iterations=iterations,
        status=status, condition=condition)
end

# pi delta / sin(pi delta), using its series for small delta; _pi_delta_jet gives the jet.
@inline function _pi_delta_ratio(delta)
    dc = ComplexF64(delta)
    if abs(dc) <= 1.0e-3
        x2 = (pi * dc)^2
        return 1 + x2 / 6 + 7 * x2^2 / 360 + 31 * x2^3 / 15120
    end
    return pi * dc / sin(pi * dc)
end

# Divides a jet by delta, dropping its constant term.
@inline function _jet_drop_constant(a::SVector{N,ComplexF64}) where {N}
    out = MVector{N,ComplexF64}(undef)
    @inbounds for n in 1:(N - 1)
        out[n] = a[n + 1]
    end
    out[N] = 0
    return SVector{N,ComplexF64}(out)
end

@inline _pi_delta_jet() = _pi_delta_jet(Val(U_JET_SIZE))

@inline function _pi_delta_jet(::Val{N}) where {N}
    out = MVector{N,ComplexF64}(ntuple(_ -> ComplexF64(0), N))
    out[1] = 1
    N >= 3 && (out[3] = pi^2 / 6)
    N >= 5 && (out[5] = 7pi^4 / 360)
    N >= 7 && (out[7] = 31pi^6 / 15120)
    return SVector{N,ComplexF64}(out)
end

# Condition number of the divided difference (component_a - component_b) / delta.
@inline function _divided_condition(
        component_a::SVector{N,ComplexF64},
        component_b::SVector{N,ComplexF64},
        difference::SVector{N,ComplexF64}, delta) where {N}
    value = ComplexF64(0)
    power = ComplexF64(1)
    weighted = 0.0
    @inbounds for j in 1:(N - 1)
        value += difference[j + 1] * power
        weighted += (abs(component_a[j + 1]) + abs(component_b[j + 1])) * abs(power)
        power *= delta
    end
    return weighted / max(abs(value), floatmin(Float64))
end

# U(a, b, z) and its first two z-derivatives for b = integer_b + delta near an integer
# and a = a0 + delta/2, from the delta-jets of the two regularized M terms of the
# connection formula, whose difference vanishes at delta = 0. Returns nothing when
# delta or z is out of range or a jet fails.
function _near_integer_triplet(a0, integer_b::Int, delta, z,
        backend::HyperUBackend)
    size = abs(ComplexF64(delta)) <= U_SMALL_DELTA ?
        Val(U_SMALL_JET_SIZE) : Val(U_JET_SIZE)
    return _near_integer_triplet(a0, integer_b, delta, z, backend, size)
end

function _near_integer_triplet(a0, integer_b::Int, delta, z,
        backend::HyperUBackend, size::Val{N}) where {N}
    degree = N - 1
    dc = ComplexF64(delta)
    zc = ComplexF64(z)
    abs(dc) <= backend.near_integer_max_delta || return nothing
    abs(zc) <= backend.near_integer_max_abs_z || return nothing
    iszero(zc) && return nothing

    rho = ComplexF64(0.5)
    a0c = ComplexF64(a0)
    m_a = _mreg_triplet(a0c, rho, integer_b, 1, zc, backend, size)
    m_b = _mreg_triplet(1 + a0c - integer_b, rho - 1, 2 - integer_b, -1, zc,
        backend, size)
    m_a.status == :OK || return nothing
    m_b.status == :OK || return nothing

    gamma_a = _rgamma_jet(1 + a0c - integer_b, rho - 1, size)
    A = _jet_mul(m_a.jet, gamma_a)
    Az = _jet_mul(m_a.dz, gamma_a)
    Azz = _jet_mul(m_a.dzz, gamma_a)

    gamma_b = _rgamma_jet(a0c, rho, size)
    C = _jet_mul(m_b.jet, gamma_b)
    Cz = _jet_mul(m_b.dz, gamma_b)
    Czz = _jet_mul(m_b.dzz, gamma_b)
    delta_jet = _jet_linear(0, 1, size)
    delta_plus_one = _jet_linear(1, 1, size)
    z_delta = _jet_exp(_jet_linear(0, -log(zc), size))
    C0 = _jet_mul(z_delta, C)
    C1 = _jet_mul(z_delta,
        Cz - _jet_mul(delta_jet, C) / zc)
    C2 = _jet_mul(z_delta,
        Czz - 2 * _jet_mul(delta_jet, Cz) / zc +
        _jet_mul(_jet_mul(delta_jet, delta_plus_one), C) / zc^2)

    p = 1 - integer_b
    B0 = C0
    B1 = C1 + p * C0 / zc
    B2 = C2 + 2p * C1 / zc + p * (p - 1) * C0 / zc^2
    log_scale = p * log(zc)
    if real(log_scale) > 0
        inverse_scale = exp(-log_scale)
        component_a = inverse_scale * A
        component_az = inverse_scale * Az
        component_azz = inverse_scale * Azz
        component_b = B0
        component_bz = B1
        component_bzz = B2
    else
        scale = exp(log_scale)
        component_a = A
        component_az = Az
        component_azz = Azz
        component_b = scale * B0
        component_bz = scale * B1
        component_bzz = scale * B2
        log_scale = ComplexF64(0)
    end

    difference = component_a - component_b
    difference_z = component_az - component_bz
    difference_zz = component_azz - component_bzz

    ratio = _pi_delta_jet(size)
    sign = isodd(integer_b) ? -1.0 : 1.0
    u_jet = sign * _jet_mul(ratio, _jet_drop_constant(difference))
    uz_jet = sign * _jet_mul(ratio, _jet_drop_constant(difference_z))
    uzz_jet = sign * _jet_mul(ratio, _jet_drop_constant(difference_zz))
    phase = cis(imag(log_scale))
    value = phase * _jet_eval(u_jet, dc)
    deriv = phase * _jet_eval(uz_jet, dc)
    second = phase * _jet_eval(uzz_jet, dc)
    if !(_finite_complex(value) && _finite_complex(deriv) &&
            _finite_complex(second) && !iszero(value))
        return nothing
    end
    scale = max(abs(value), abs(deriv), abs(second), floatmin(Float64))
    tail = max(
        abs(u_jet[N - 1] * dc^(degree - 1)),
        abs(uz_jet[N - 1] * dc^(degree - 1)),
        abs(uzz_jet[N - 1] * dc^(degree - 1)),
    ) / scale
    condition = max(
        m_a.condition,
        m_b.condition,
        _divided_condition(component_a, component_b, difference, dc),
        _divided_condition(component_az, component_bz, difference_z, dc),
        _divided_condition(component_azz, component_bzz, difference_zz, dc),
    )
    actual_a = a0c + rho * dc
    actual_b = ComplexF64(integer_b) + dc
    ode_terms = (
        zc * second,
        (actual_b - zc) * deriv,
        -actual_a * value,
    )
    ode_residual = abs(sum(ode_terms)) /
        max(sum(abs, ode_terms), floatmin(Float64))
    relerr = max(
        m_a.relerr,
        m_b.relerr,
        tail,
        eps(Float64) * condition,
        ode_residual,
    )
    status = relerr <= backend.reltol ? :OK : :WARN_CONDITIONED
    return (
        value=ComplexF64(value),
        deriv=ComplexF64(deriv),
        second=ComplexF64(second),
        log_scale=Float64(real(log_scale)),
        estimated_relerr=Float64(relerr),
        iterations=max(m_a.iterations, m_b.iterations),
        condition=Float64(condition),
        ode_residual=Float64(ode_residual),
        status=status,
    )
end

# log U(a, b, z) by the same near-integer jet construction (value only).
function _near_integer_u(a0, integer_b::Int, delta, z, backend::HyperUBackend)
    dc = ComplexF64(delta)
    zc = ComplexF64(z)
    abs(dc) <= backend.near_integer_max_delta ||
        return _hfail(:log_u_near_integer_jet, "outside near-integer delta range")
    abs(zc) <= backend.near_integer_max_abs_z ||
        return _hfail(:log_u_near_integer_jet, "outside near-integer argument range")
    iszero(zc) && return _hfail(:log_u_near_integer_jet, "zero argument")

    rho = ComplexF64(0.5)
    a0c = ComplexF64(a0)
    m_a = _mreg_jet(a0c, rho, integer_b, 1, zc, backend)
    m_b = _mreg_jet(1 + a0c - integer_b, rho - 1, 2 - integer_b, -1, zc, backend)
    m_a.status == :OK ||
        return _hfail(:log_u_near_integer_jet, "first regularized M jet failed")
    m_b.status == :OK ||
        return _hfail(:log_u_near_integer_jet, "second regularized M jet failed")

    a_jet = _jet_mul(
        m_a.jet,
        _rgamma_jet(1 + a0c - integer_b, rho - 1),
    )
    b_core = _jet_mul(
        m_b.jet,
        _rgamma_jet(a0c, rho),
    )
    z_delta = _jet_exp(_jet_linear(0, -log(zc)))
    b_core = _jet_mul(z_delta, b_core)

    log_scale = (1 - integer_b) * log(zc)
    if real(log_scale) > 0
        a_scaled = exp(-log_scale) * a_jet
        difference = a_scaled - b_core
        component_a = a_scaled
        component_b = b_core
    else
        b_scaled = exp(log_scale) * b_core
        difference = a_jet - b_scaled
        component_a = a_jet
        component_b = b_scaled
        log_scale = ComplexF64(0)
    end

    polynomial = ComplexF64(0)
    power = ComplexF64(1)
    weighted_components = 0.0
    last_term = ComplexF64(0)
    @inbounds for j in 1:U_JET_DEGREE
        contribution = difference[j + 1] * power
        polynomial += contribution
        weighted_components +=
            (abs(component_a[j + 1]) + abs(component_b[j + 1])) * abs(power)
        last_term = contribution
        power *= dc
    end
    if !_finite_complex(polynomial) || iszero(polynomial)
        return _hfail(:log_u_near_integer_jet,
            "near-integer jet returned zero or nonfinite")
    end

    sign = isodd(integer_b) ? -1.0 : 1.0
    finite_part = sign * _pi_delta_ratio(dc) * polynomial
    if !_finite_complex(finite_part) || iszero(finite_part)
        return _hfail(:log_u_near_integer_jet,
            "near-integer finite part is zero or nonfinite")
    end
    condition = weighted_components / max(abs(polynomial), floatmin(Float64))
    tail = abs(last_term) / max(abs(polynomial), floatmin(Float64))
    relerr = max(m_a.relerr, m_b.relerr, tail, eps(Float64) * condition)
    status = relerr <= 100 * backend.reltol ? :OK : :WARN_CONDITIONED
    return HyperEval(
        ComplexF64(log_scale + log(finite_part)),
        :log_u_near_integer_jet,
        relerr,
        max(m_a.iterations, m_b.iterations),
        status,
        "integer_b=$(integer_b), delta=$(dc), condition=$(condition)",
    )
end

function _factorial_float(n::Integer)
    n < 0 && throw(ArgumentError("factorial needs nonnegative n"))
    value = 1.0
    for k in 2:n
        value *= k
    end
    return value
end

function _pochhammer_product(a, n::Integer)
    n < 0 && throw(ArgumentError("pochhammer needs nonnegative n"))
    value = ComplexF64(1)
    ac = ComplexF64(a)
    for k in 0:(n - 1)
        value *= ac + k
    end
    return value
end

function _integer_snap(z, tol)
    zc = ComplexF64(z)
    abs(imag(zc)) <= tol || return nothing
    rounded = round(Int, real(zc))
    abs(real(zc) - rounded) <= tol || return nothing
    return rounded
end

function _near_nonpositive_integer(z; tol=1.0e-14)
    snapped = _integer_snap(z, tol)
    snapped === nothing && return nothing
    snapped <= 0 || return nothing
    return -snapped
end

function _invgamma_complex(z)
    zc = ComplexF64(z)
    _near_nonpositive_integer(zc) !== nothing && return ComplexF64(0)
    try
        gz = gamma(zc)
        if _finite_complex(gz) && !iszero(gz)
            return inv(gz)
        end
    catch
    end
    try
        value = exp(-loggamma(zc))
        _finite_complex(value) && return ComplexF64(value)
    catch
    end
    return ComplexF64(NaN + NaN * im)
end

function _onef1_series(a, b, z, backend::HyperUBackend)
    term = ComplexF64(1)
    total = ComplexF64(1)
    maxterm = abs(term)
    for n in 0:(backend.series_max - 1)
        den = (b + n) * (n + 1)
        _near_zero(den) && return _hfail(:onef1_series, "series pole in denominator")
        term *= (a + n) * z / den
        total += term
        maxterm = max(maxterm, abs(term))
        rel = _hrelerr(term, total)
        if rel <= backend.reltol
            return _heval(total, :onef1_series,
                max(rel, eps(Float64) * maxterm / max(abs(total), floatmin(Float64))), n + 1)
        end
    end
    return HyperEval(total, :onef1_series, _hrelerr(term, total),
        backend.series_max, :WARN_MAXITER, "series did not meet reltol")
end

function _onef1_direct(a, b, z)
    try
        value = ComplexF64(HGF_1F1(a, b, z))
        _finite_complex(value) || return _hfail(:onef1_direct, "library returned nonfinite")
        return _heval(value, :onef1_direct, eps(Float64), 0)
    catch err
        return _hfail(:onef1_direct, sprint(showerror, err))
    end
end

function _onef1(a, b, z, backend::HyperUBackend)
    direct = _onef1_direct(a, b, z)
    direct.status == :OK && return direct
    abs(z) <= 40 && return _onef1_series(a, b, z, backend)
    return direct
end

# 1F1(a; b; z) Taylor series, after Kummer's transformation when Re z < 0 or |b - a| < |a|.
function _onef1_taylor(a_temp, b_temp, z_temp, backend::HyperUBackend)
    a = ComplexF64(a_temp)
    b = ComplexF64(b_temp)
    z = ComplexF64(z_temp)
    log_prefactor = ComplexF64(0)
    if real(z) < 0 || abs(b - a) < abs(a)
        a = b - a
        z = -z
        log_prefactor = ComplexF64(z_temp)
    end

    term = ComplexF64(1)
    total = ComplexF64(0)
    maxterm = abs(term)
    s = 0
    min_terms = max(4, ceil(Int, max(0.0, -real(a))), ceil(Int, max(0.0, -real(b))))
    while s < min_terms
        total += term
        s += 1
        den = (b + (s - 1)) * s
        _near_zero(den) && return _hfail(:onef1_taylor, "series pole in denominator")
        term *= z * (a + (s - 1)) / den
        maxterm = max(maxterm, abs(term))
    end

    rel = _hrelerr(term, total)
    while rel > eps(Float64) / 100 && s < min(backend.series_max, 500)
        total += term
        s += 1
        den = (b + (s - 1)) * s
        _near_zero(den) && return _hfail(:onef1_taylor, "series pole in denominator")
        term *= z * (a + (s - 1)) / den
        maxterm = max(maxterm, abs(term))
        rel = _hrelerr(term, total)
    end
    total += term

    value = exp(log_prefactor) * total
    err = max(_hrelerr(term, total), eps(Float64) * maxterm / max(abs(total), floatmin(Float64)))
    status = err <= backend.reltol ? :OK : :WARN_MAXITER
    return HyperEval(value, :onef1_taylor, err, s, status, "")
end

function _onef1_for_u(a, b, z, backend::HyperUBackend)
    taylor = _onef1_taylor(a, b, z, backend)
    taylor.status in (:OK, :WARN_MAXITER) && return taylor
    return _onef1(a, b, z, backend)
end

_gamma_ratio(a, b) = exp(loggamma(ComplexF64(a)) - loggamma(ComplexF64(b)))

# log U(a, b, z) from its connection formula in two 1F1 functions.
function _u_connection(a, b, z, backend::HyperUBackend)
    m1 = _onef1_for_u(a, b, z, backend)
    m2 = _onef1_for_u(a - b + 1, 2 - b, z, backend)
    if !(m1.status in (:OK, :WARN_MAXITER)) || !(m2.status in (:OK, :WARN_MAXITER))
        return _hfail(:log_u_connection, "1F1 subroute failed")
    end
    try
        term1 = _gamma_ratio(1 - b, a - b + 1) * m1.value
        term2 = _gamma_ratio(b - 1, a) * z^(1 - b) * m2.value
        value = term1 + term2
        if !_finite_complex(value) || iszero(value)
            return _hfail(:log_u_connection, "connection formula returned zero or nonfinite")
        end
        cancellation = (abs(term1) + abs(term2)) / max(abs(value), floatmin(Float64))
        status = cancellation > 1.0e10 ? :WARN_CONDITIONED : :OK
        relerr = max(m1.estimated_relerr, m2.estimated_relerr, eps(Float64) * cancellation)
        return HyperEval(log(value), :log_u_connection, relerr,
            max(m1.iterations, m2.iterations), status, "cancellation=$(cancellation)")
    catch err
        return _hfail(:log_u_connection, sprint(showerror, err))
    end
end

function _u_direct(a, b, z)
    try
        value = ComplexF64(HypergeometricFunctions.U(a, b, z))
        if !_finite_complex(value) || iszero(value)
            return _hfail(:log_u_direct, "library returned zero or nonfinite")
        end
        return HyperEval(log(value), :log_u_direct, 1.0e-8, 0, :WARN_CONDITIONED,
            "direct library fallback")
    catch err
        return _hfail(:log_u_direct, sprint(showerror, err))
    end
end

# log U near a = b = 1, where U(1, 1, z) = exp(z) E1(z).
function _u_one_one_limit(a, b, z, backend::HyperUBackend)
    ac = ComplexF64(a)
    bc = ComplexF64(b)
    zc = ComplexF64(z)
    da = abs(ac - 1)
    db = abs(bc - 1)
    max(da, db) <= backend.u11_limit_tol ||
        return _hfail(:log_u_one_one_limit, "outside U(1,1,z) tolerance")
    try
        value = exp(zc) * expint(1, zc)
        if !_finite_complex(value) || iszero(value)
            return _hfail(:log_u_one_one_limit, "exp(z)E1(z) returned zero or nonfinite")
        end
        exact = max(da, db) <= 1.0e-14
        status = exact ? :OK : :WARN_LIMIT_APPROX
        relerr = exact ? eps(Float64) :
            max(da, db) * max(1.0, abs(log(zc)), abs(zc)) / max(abs(value), floatmin(Float64))
        return HyperEval(log(value), :log_u_one_one_limit, relerr, 0, status, "")
    catch err
        return _hfail(:log_u_one_one_limit, sprint(showerror, err))
    end
end

# log U(a, n + 1, z) for a nonpositive integer a, where U is a polynomial in z.
function _u_polynomial_a(a, n::Int, z, exact_b::Bool)
    m = _near_nonpositive_integer(a)
    m === nothing && return nothing
    zc = ComplexF64(z)
    total = ComplexF64(0)
    for s in 0:m
        term = binomial(m, s) * _pochhammer_product(n + s + 1, m - s) * (-zc)^s
        total += term
    end
    value = (isodd(m) ? -1.0 : 1.0) * total
    if !_finite_complex(value) || iszero(value)
        return _hfail(:log_u_integer_b_polynomial, "polynomial route returned zero or nonfinite")
    end
    status = exact_b ? :OK : :WARN_LIMIT_APPROX
    return HyperEval(log(value), :log_u_integer_b_polynomial, eps(Float64), m + 1, status, "")
end

# log U(a, n + 1, z) for integer n >= 0 from the logarithmic formula DLMF 13.2.9.
function _u_positive_integer_b(a, n::Int, z, backend::HyperUBackend, exact_b::Bool)
    n < 0 && return _hfail(:log_u_integer_b_limit, "positive-integer-b route requires n >= 0")
    ac = ComplexF64(a)
    zc = ComplexF64(z)

    polynomial = _u_polynomial_a(ac, n, zc, exact_b)
    polynomial !== nothing && return polynomial

    try
        finite_sum = ComplexF64(0)
        for k in 1:n
            finite_sum += _factorial_float(k - 1) *
                _pochhammer_product(1 - ac + k, n - k) /
                _factorial_float(n - k) * zc^(-k)
        end
        finite_part = _invgamma_complex(ac) * finite_sum
        prefactor = (isodd(n + 1) ? -1.0 : 1.0) *
            _invgamma_complex(ac - n) / _factorial_float(n)
        if !_finite_complex(prefactor) || !_finite_complex(finite_part)
            return _hfail(:log_u_integer_b_limit, "integer-b prefactor or finite sum is nonfinite")
        end

        coeff = ComplexF64(1)
        zpow = ComplexF64(1)
        series_sum = ComplexF64(0)
        maxterm = 0.0
        last_rel = Inf
        iterations = 0
        converged = false
        logz = log(zc)
        for k in 0:(backend.series_max - 1)
            bracket = logz + digamma(ac + k) - digamma(1 + k) - digamma(n + k + 1)
            term = coeff * zpow * bracket
            series_sum += term
            maxterm = max(maxterm, abs(term))
            scale = max(abs(series_sum), abs(finite_part), floatmin(Float64))
            last_rel = abs(term) / scale
            iterations = k + 1
            if k >= 2 && last_rel <= backend.reltol
                converged = true
                break
            end
            coeff *= (ac + k) / ((n + 1 + k) * (k + 1))
            zpow *= zc
        end

        infinite_part = prefactor * series_sum
        value = infinite_part + finite_part
        if !_finite_complex(value) || iszero(value)
            return _hfail(:log_u_integer_b_limit, "integer-b formula returned zero or nonfinite")
        end
        cancellation = (abs(infinite_part) + abs(finite_part)) / max(abs(value), floatmin(Float64))
        relerr = max(last_rel,
            eps(Float64) * cancellation * max(1.0, maxterm / max(abs(series_sum), floatmin(Float64))))
        status =
            !converged ? :WARN_MAXITER :
            cancellation > 1.0e10 ? :WARN_CONDITIONED :
            exact_b ? :OK : :WARN_LIMIT_APPROX
        return HyperEval(log(value), :log_u_integer_b_limit, relerr, iterations, status,
            "b_n=$(n), cancellation=$(cancellation)")
    catch err
        return _hfail(:log_u_integer_b_limit, sprint(showerror, err))
    end
end

# log U for b within integer_b_tol of an integer; b <= 0 goes through
# U(a, b, z) = z^(1 - b) U(a - b + 1, 2 - b, z).
function _u_integer_b(a, b, z, backend::HyperUBackend)
    bc = ComplexF64(b)
    snapped_b = _integer_snap(bc, backend.integer_b_tol)
    snapped_b === nothing && return _hfail(:log_u_integer_b_limit, "b outside integer tolerance")
    exact_b = abs(bc - snapped_b) <= 1.0e-14
    if snapped_b >= 1
        return _u_positive_integer_b(a, snapped_b - 1, z, backend, exact_b)
    end

    shift_n = -snapped_b
    shifted = _u_positive_integer_b(ComplexF64(a) + shift_n + 1, shift_n + 1, z, backend, exact_b)
    _usable(shifted.status) || return shifted
    route = shifted.route == :log_u_integer_b_polynomial ?
        :log_u_integer_b_shift_polynomial : :log_u_integer_b_shift_limit
    return HyperEval((shift_n + 1) * log(ComplexF64(z)) + shifted.value, route,
        shifted.estimated_relerr, shifted.iterations, shifted.status, shifted.message)
end

# Large-|z| asymptotic series of U, stopped at its smallest term if it does not converge.
function _u_asymptotic(a, b, z, backend::HyperUBackend)
    term = ComplexF64(1)
    total = ComplexF64(1)
    best_total = total
    best_rel = Inf
    previous_abs = abs(term)
    for n in 0:(backend.asymptotic_max - 1)
        term *= -((a + n) * (a - b + 1 + n)) / ((n + 1) * z)
        total += term
        rel = _hrelerr(term, total)
        if rel < best_rel
            best_rel = rel
            best_total = total
        end
        if rel <= backend.reltol
            return _heval(log(z^(-a) * total), :log_u_asymptotic, rel, n + 1)
        end
        current_abs = abs(term)
        if n > 1 && current_abs > previous_abs
            value = z^(-a) * best_total
            return HyperEval(log(value), :log_u_asymptotic_leastterm, best_rel,
                n + 1, :WARN_ASYMPTOTIC_LEASTTERM, "least-term stop")
        end
        previous_abs = current_abs
    end
    value = z^(-a) * best_total
    return HyperEval(log(value), :log_u_asymptotic_maxiter, best_rel,
        backend.asymptotic_max, :WARN_MAXITER, "asymptotic maxiter")
end

# Ordering of U candidates: usable status, integer-b routes at small |z|, estimated
# error, status, and a penalty for an unconverged asymptotic series at |z| < 10.
function _u_rank(candidate::HyperEval, z, backend::HyperUBackend)
    usable_rank = _usable(candidate.status) ? 0 : 1
    integer_priority =
        (candidate.route in (:log_u_integer_b_limit, :log_u_integer_b_shift_limit,
                             :log_u_integer_b_polynomial, :log_u_integer_b_shift_polynomial) &&
         candidate.status in (:OK, :WARN_LIMIT_APPROX) &&
         abs(z) <= backend.integer_b_priority_max_abs_z) ? -1 : 0
    relerr_rank = isfinite(candidate.estimated_relerr) ? candidate.estimated_relerr : Inf
    status_rank =
        candidate.status == :OK ? 0 :
        candidate.status == :WARN_LIMIT_APPROX ? 1 :
        candidate.status == :WARN_CONDITIONED ? 2 : 3
    asymptotic_penalty =
        (candidate.route in (:log_u_asymptotic_leastterm, :log_u_asymptotic_maxiter) && abs(z) < 10) ? 1 : 0
    return (usable_rank, integer_priority, relerr_rank, status_rank, asymptotic_penalty)
end

function _preferred_integer_u(candidate::HyperEval, z, backend::HyperUBackend)
    return candidate.route in (
        :log_u_integer_b_limit,
        :log_u_integer_b_shift_limit,
        :log_u_integer_b_polynomial,
        :log_u_integer_b_shift_polynomial,
    ) && candidate.status in (:OK, :WARN_LIMIT_APPROX) &&
        abs(z) <= backend.integer_b_priority_max_abs_z
end

function _better_u(best, candidate::HyperEval, z, backend::HyperUBackend)
    _usable(candidate.status) || return best
    best === nothing && return candidate
    return _u_rank(candidate, z, backend) < _u_rank(best, z, backend) ? candidate : best
end

# log U(a, b, z) by the best usable route: special cases, integer-b formulas, the
# U(1, 1, z) limit, the asymptotic series, the connection formula, then the library.
function _eval_log_u(backend::HyperUBackend, a, b, z; integer_route::Bool=true)
    ac = ComplexF64(a)
    bc = ComplexF64(b)
    zc = ComplexF64(z)
    _near_zero(ac) && return _heval(0, :log_u_special_a_zero, 0, 0)
    _near_zero(ac - (bc - 1)) &&
        return _heval(-ac * log(zc), :log_u_special_a_eq_b_minus_1, eps(Float64), 0)

    intb = integer_route ? _u_integer_b(ac, bc, zc, backend) :
        _hfail(:log_u_integer_b_limit, "integer route disabled")
    integer_route && _preferred_integer_u(intb, zc, backend) && return intb
    best = integer_route ? _better_u(nothing, intb, zc, backend) : nothing
    oneone = _u_one_one_limit(ac, bc, zc, backend)
    best = _better_u(best, oneone, zc, backend)
    if abs(zc) > backend.asymptotic_min_abs_z
        asym = _u_asymptotic(ac, bc, zc, backend)
        best = _better_u(best, asym, zc, backend)
    end
    conn = _u_connection(ac, bc, zc, backend)
    best = _better_u(best, conn, zc, backend)
    best === nothing || return best
    direct = _u_direct(ac, bc, zc)
    _usable(direct.status) && return direct
    return _hfail(:log_u_no_route, conn.message)
end

# log U at a = a0 + delta/2, b = integer_b + delta: the near-integer jet when it succeeds,
# otherwise the more accurate of it and the general evaluator. _mst_log_du gives
# log dU/dz = log(-a U(a + 1, b + 1, z)).
function _mst_log_u(a0, integer_b::Int, delta, z)
    near = _near_integer_u(a0, integer_b, delta, z, HYPER_U)
    near.status == :OK && return near
    actual_a = ComplexF64(a0) + ComplexF64(delta) / 2
    actual_b = ComplexF64(integer_b) + ComplexF64(delta)
    outer = _eval_log_u(HYPER_U, actual_a, actual_b, z; integer_route=false)
    !_usable(near.status) && return outer
    !_usable(outer.status) && return near
    return near.estimated_relerr <= outer.estimated_relerr ? near : outer
end

function _mst_log_du(a0, integer_b::Int, delta, z)
    actual_a = ComplexF64(a0) + ComplexF64(delta) / 2
    _near_zero(actual_a) && return HyperEval(
        ComplexF64(-Inf),
        :log_du_a_zero,
        0.0,
        0,
        :OK,
        "dU/dz is exactly zero for a=0",
    )
    shifted = _mst_log_u(ComplexF64(a0) + 1, integer_b + 1, delta, z)
    if _usable(shifted.status)
        return HyperEval(
            shifted.value + log(-actual_a),
            Symbol("log_du_", shifted.route),
            shifted.estimated_relerr,
            shifted.iterations,
            shifted.status,
            shifted.message,
        )
    end
    return _hfail(:log_du_no_route, shifted.message)
end

# A U-function state (value, first and second z-derivatives) normalized to unit size,
# with its log scale.
struct MSTHState
    value::ComplexF64
    deriv::ComplexF64
    second::ComplexF64
    log_scale::Float64
end

# The U-function states h_n, |n| <= MST_N_MAX, of one infinity series at a fixed argument,
# filled lazily by three-term recurrence with direct evaluation as anchor and fallback,
# plus the status of the Miller and negative-n constructions.
mutable struct MSTUSequence
    a0::ComplexF64
    integer_b::Int
    delta::ComplexF64
    argument::ComplexF64
    states::Vector{MSTHState}
    ready::BitVector
    conditions::Vector{Float64}
    errors::Vector{Float64}
    miller_ready::Bool
    miller_nmax::Int
    miller_cutoff::Int
    miller_error::Float64
    miller_status::Symbol
    negative_status::Symbol
    negative_error::Float64
end

@inline _seq_index(n::Integer) = Int(n) + MST_N_MAX + 1

@inline function _term_cond(total, terms...)
    scale = 0.0
    @inbounds for term in terms
        scale = max(scale, abs(term))
    end
    return scale / max(abs(total), floatmin(Float64))
end

@inline function _scaled_product(log_factor, value)
    iszero(value) && return ComplexF64(0)
    return exp(ComplexF64(log_factor) + log(ComplexF64(value)))
end

function _normalize_h(value, deriv, second, log_scale)
    norm = max(abs(value), abs(deriv), abs(second))
    if !isfinite(norm) || iszero(norm)
        nan = ComplexF64(NaN + NaN * im)
        return MSTHState(nan, nan, nan, NaN)
    end
    return MSTHState(
        ComplexF64(value / norm),
        ComplexF64(deriv / norm),
        ComplexF64(second / norm),
        Float64(log_scale + log(norm)),
    )
end

# h_n = c^n U(a0 + n, integer_b + 2n + delta, c) and its c-derivatives from the
# near-integer jet; returns (state, jet result) or nothing.
function _near_h(sequence::MSTUSequence, n::Integer)
    a0 = sequence.a0 + n
    integer_b = sequence.integer_b + 2n
    c = sequence.argument
    triplet = _near_integer_triplet(a0, integer_b, sequence.delta, c, HYPER_U)
    triplet === nothing && return nothing
    log_power = n * log(c)
    phase = cis(imag(log_power))
    value = phase * triplet.value
    deriv = phase * (n * triplet.value / c + triplet.deriv)
    second = phase * (
        n * (n - 1) * triplet.value / c^2 +
        2n * triplet.deriv / c + triplet.second
    )
    state = _normalize_h(
        value,
        deriv,
        second,
        triplet.log_scale + real(log_power),
    )
    if !(_finite_complex(state.value) && _finite_complex(state.deriv) &&
            _finite_complex(state.second) && isfinite(state.log_scale))
        return nothing
    end
    return state, triplet
end

# h_n from log U and log dU of the general evaluator, with its error estimate.
function _legacy_h_eval(sequence::MSTUSequence, n::Integer)
    a0 = sequence.a0 + n
    integer_b = sequence.integer_b + 2n
    delta = sequence.delta
    c = sequence.argument
    log_u = _mst_log_u(a0, integer_b, delta, c)
    log_du = _mst_log_du(a0, integer_b, delta, c)
    if !_usable(log_u.status) || !_usable(log_du.status)
        nan = ComplexF64(NaN + NaN * im)
        return MSTHState(nan, nan, nan, NaN), Inf
    end
    actual_a = a0 + delta / 2
    actual_b_minus_c = (ComplexF64(integer_b) - c) + delta
    du_over_u = _near_zero(actual_a) ? ComplexF64(0) : exp(log_du.value - log_u.value)
    d2_over_u = (actual_a - actual_b_minus_c * du_over_u) / c
    deriv_ratio = n / c + du_over_u
    second_ratio = n * (n - 1) / c^2 + 2n * du_over_u / c + d2_over_u
    log_h = n * log(c) + log_u.value
    scale = real(log_h) + log(max(1.0, abs(deriv_ratio), abs(second_ratio)))
    h = exp(log_h - scale)
    state = _normalize_h(h, h * deriv_ratio, h * second_ratio, scale)
    eta_u = max(log_u.estimated_relerr, eps(Float64))
    eta_du = _near_zero(actual_a) ? 0.0 :
        max(log_du.estimated_relerr, eps(Float64))
    deriv_error = (
        abs(n / c) * eta_u + abs(du_over_u) * eta_du
    ) / max(abs(deriv_ratio), floatmin(Float64))
    c0 = n * (n - 1) / c^2 + actual_a / c
    c1 = 2n / c - actual_b_minus_c / c
    second_error = (
        abs(c0) * eta_u + abs(c1 * du_over_u) * eta_du
    ) / max(abs(second_ratio), floatmin(Float64))
    error = max(eta_u, deriv_error, second_error, eps(Float64))
    return state, Float64(error)
end

# Direct h_n: the near-integer jet when it is accurate, else whichever of it and the
# general evaluator has the smaller error.
function _direct_h_eval(sequence::MSTUSequence, n::Integer)
    near = _near_h(sequence, n)
    if near !== nothing
        near_state, near_result = near
        near_error = max(Float64(near_result.estimated_relerr), eps(Float64))
        near_result.status == :OK && return near_state, near_error
        legacy_state, legacy_error = _legacy_h_eval(sequence, n)
        if !isfinite(legacy_state.log_scale) || near_error <= legacy_error
            return near_state, near_error
        end
        return legacy_state, legacy_error
    end
    return _legacy_h_eval(sequence, n)
end

function MSTUSequence(a0, integer_b::Int, delta, argument)
    a0c = ComplexF64(a0)
    dc = ComplexF64(delta)
    c = ComplexF64(argument)
    iszero(c) && error("MST U sequence requires a nonzero argument")
    states = Vector{MSTHState}(undef, 2MST_N_MAX + 1)
    ready = falses(length(states))
    conditions = fill(Inf, length(states))
    errors = fill(Inf, length(states))
    sequence = MSTUSequence(
        a0c,
        integer_b,
        dc,
        c,
        states,
        ready,
        conditions,
        errors,
        false,
        0,
        0,
        Inf,
        :UNTRIED,
        :UNTRIED,
        Inf,
    )
    for n in (0, 1)
        index = _seq_index(n)
        states[index], errors[index] = _direct_h_eval(sequence, n)
        isfinite(states[index].log_scale) ||
            error("MST U sequence could not construct finite anchors")
        ready[index] = true
        conditions[index] = 1.0
    end
    return sequence
end

@inline _h_norm(state::MSTHState) =
    max(abs(state.value), abs(state.deriv), abs(state.second))

# Coefficients of the three-term recurrence of h_n and the c-derivatives of the middle one.
@inline function _h_coeffs(sequence::MSTUSequence, n::Int)
    a0 = sequence.a0
    integer_b = sequence.integer_b
    delta = sequence.delta
    c = sequence.argument
    b2n = ComplexF64(integer_b + 2n) + delta
    q0 = 8 + ComplexF64(integer_b + 2n)^2 +
        2 * (a0 + n) * c - ComplexF64(integer_b + 2n) * (6 + c)
    q = q0 + delta * (2 * (integer_b + 2n) - 6) + delta^2
    qprime = 2 * a0 - integer_b
    qconstant = q - qprime * c
    A = (ComplexF64(integer_b + n - 2) - a0 + delta / 2) * (b2n - 2)
    B = (b2n - 3) * q / c
    Bprime = -(b2n - 3) * qconstant / c^2
    Bsecond = 2 * (b2n - 3) * qconstant / c^3
    D = (a0 + n - 1 + delta / 2) * (b2n - 4)
    return A, B, Bprime, Bsecond, D
end

# One backward step of the h_n recurrence with its condition number.
function _miller_step(sequence::MSTUSequence, n::Int, hn::MSTHState, hn1::MSTHState)
    A, B, Bprime, Bsecond, D = _h_coeffs(sequence, n)
    abs_A = abs(A)
    if !isfinite(abs_A) || iszero(abs_A)
        nan = ComplexF64(NaN + NaN * im)
        return MSTHState(nan, nan, nan, NaN), Inf
    end
    weight_n = max(abs(D) * _h_norm(hn), floatmin(Float64))
    weight_n1 = max(
        max(abs(B), abs(Bprime), abs(Bsecond)) * _h_norm(hn1),
        floatmin(Float64),
    )
    aligned_scale = max(
        hn.log_scale + log(weight_n),
        hn1.log_scale + log(weight_n1),
    )
    if !isfinite(aligned_scale)
        nan = ComplexF64(NaN + NaN * im)
        return MSTHState(nan, nan, nan, NaN), Inf
    end
    scale_n = exp(hn.log_scale - aligned_scale)
    scale_n1 = exp(hn1.log_scale - aligned_scale)
    value_terms = (
        D * hn.value * scale_n,
        -B * hn1.value * scale_n1,
    )
    deriv_terms = (
        D * hn.deriv * scale_n,
        -Bprime * hn1.value * scale_n1,
        -B * hn1.deriv * scale_n1,
    )
    second_terms = (
        D * hn.second * scale_n,
        -Bsecond * hn1.value * scale_n1,
        -2 * Bprime * hn1.deriv * scale_n1,
        -B * hn1.second * scale_n1,
    )
    value = value_terms[1] + value_terms[2]
    deriv = deriv_terms[1] + deriv_terms[2] + deriv_terms[3]
    second = sum(second_terms)
    condition = max(
        _term_cond(value, value_terms...),
        _term_cond(deriv, deriv_terms...),
        _term_cond(second, second_terms...),
    )
    phase_A = A / abs_A
    state = _normalize_h(
        value / phase_A,
        deriv / phase_A,
        second / phase_A,
        aligned_scale - log(abs_A),
    )
    return state, condition
end

function _h_relerr(candidate::MSTHState, reference::MSTHState, field::Symbol)
    candidate_value = getfield(candidate, field)
    reference_value = getfield(reference, field)
    if !(_finite_complex(candidate_value) && _finite_complex(reference_value) &&
            isfinite(candidate.log_scale) && isfinite(reference.log_scale))
        return Inf
    end
    scale = max(candidate.log_scale, reference.log_scale)
    scaled_candidate = candidate_value * exp(candidate.log_scale - scale)
    scaled_reference = reference_value * exp(reference.log_scale - scale)
    return abs(scaled_candidate - scaled_reference) /
        max(abs(scaled_candidate), abs(scaled_reference), floatmin(Float64))
end

function _h_adjust(state::MSTHState, phase, log_adjust, ell, q)
    return _normalize_h(
        phase * state.value,
        phase * (state.deriv + ell * state.value),
        phase * (state.second + 2 * ell * state.deriv + q * state.value),
        state.log_scale + log_adjust,
    )
end

# Miller's backward recurrence from `cutoff`, normalized to the direct h_0 and checked
# against the direct h_1; returns (states, conditions, error) or nothing.
function _miller_trial(sequence::MSTUSequence, cutoff::Int, requested_n::Int)
    cutoff > requested_n || error("Miller cutoff must exceed requested index")
    trial = Vector{MSTHState}(undef, cutoff + 2)
    conditions = fill(Inf, cutoff + 2)
    trial[cutoff + 2] = MSTHState(0, 0, 0, 0.0)
    trial[cutoff + 1] = MSTHState(1, 0, 0, 0.0)
    for n in (cutoff + 1):-1:2
        state, condition = _miller_step(sequence, n, trial[n + 1], trial[n])
        if !(_finite_complex(state.value) && _finite_complex(state.deriv) &&
                _finite_complex(state.second) && isfinite(state.log_scale))
            return nothing
        end
        trial[n - 1] = state
        conditions[n - 1] = condition
    end

    exact0 = sequence.states[_seq_index(0)]
    h0 = trial[1]
    if iszero(h0.value) || iszero(exact0.value)
        return nothing
    end
    ratio = exact0.value / h0.value
    ratio_abs = abs(ratio)
    (!isfinite(ratio_abs) || iszero(ratio_abs)) && return nothing
    phase = ratio / ratio_abs
    log_adjust = exact0.log_scale - h0.log_scale + log(ratio_abs)
    ell = exact0.deriv / exact0.value - h0.deriv / h0.value
    q = exact0.second / exact0.value - h0.second / h0.value -
        2 * ell * h0.deriv / h0.value
    if !(_finite_complex(ell) && _finite_complex(q) && isfinite(log_adjust))
        return nothing
    end
    normalized = Vector{MSTHState}(undef, requested_n + 1)
    for n in 0:requested_n
        state = trial[n + 1]
        normalized_state = _h_adjust(state, phase, log_adjust, ell, q)
        if !(_finite_complex(normalized_state.value) &&
                _finite_complex(normalized_state.deriv) &&
                _finite_complex(normalized_state.second) &&
                isfinite(normalized_state.log_scale))
            return nothing
        end
        normalized[n + 1] = normalized_state
    end
    exact1 = sequence.states[_seq_index(1)]
    error = max(
        _h_relerr(normalized[2], exact1, :value),
        _h_relerr(normalized[2], exact1, :deriv),
        _h_relerr(normalized[2], exact1, :second),
    )
    return normalized, conditions, error
end

# Fills h_2..h_requested_n by Miller's algorithm, raising the cutoff until the h_1 check
# meets MST_MILLER_TOL.
function _miller_positive!(sequence::MSTUSequence, requested_n::Int)
    sequence.miller_status == :ACCEPTED &&
        sequence.miller_nmax >= requested_n && return sequence
    start = max(requested_n + MST_MILLER_MARGIN, MST_MILLER_MARGIN)
    cutoffs = collect(start:MST_MILLER_STEP:MST_MILLER_MAX)
    if isempty(cutoffs) || last(cutoffs) != MST_MILLER_MAX
        push!(cutoffs, MST_MILLER_MAX)
    end
    best_error = Inf
    best_cutoff = 0
    for cutoff in cutoffs
        result = _miller_trial(sequence, cutoff, requested_n)
        result === nothing && continue
        states, conditions, error = result
        if error < best_error
            best_error = error
            best_cutoff = cutoff
        end
        if error <= MST_MILLER_TOL
            for n in 2:requested_n
                index = _seq_index(n)
                sequence.states[index] = states[n + 1]
                sequence.ready[index] = true
                sequence.conditions[index] = conditions[n + 1]
                sequence.errors[index] = max(
                    error,
                    sequence.errors[_seq_index(0)],
                    sequence.errors[_seq_index(1)],
                    eps(Float64) * conditions[n + 1],
                )
            end
            sequence.miller_ready = true
            sequence.miller_nmax = requested_n
            sequence.miller_cutoff = cutoff
            sequence.miller_error = error
            sequence.miller_status = :ACCEPTED
            return sequence
        end
    end
    sequence.miller_cutoff = best_cutoff
    sequence.miller_error = best_error
    sequence.miller_status = :REJECTED
    return sequence
end

# h_n from h_(n-1) and h_(n-2); when unstable, Miller's algorithm (n > 1) or direct
# evaluation. Returns (state, condition, error).
function _forward_h(sequence::MSTUSequence, n::Int)
    previous2 = _h_state(sequence, n - 2)
    previous1 = _h_state(sequence, n - 1)
    parent_error = max(
        sequence.errors[_seq_index(n - 2)],
        sequence.errors[_seq_index(n - 1)],
    )
    a0 = sequence.a0
    integer_b = sequence.integer_b
    delta = sequence.delta
    c = sequence.argument

    b2n = ComplexF64(integer_b + 2n) + delta
    q0 = 8 + ComplexF64(integer_b + 2n)^2 +
        2 * (a0 + n) * c - ComplexF64(integer_b + 2n) * (6 + c)
    q = q0 + delta * (2 * (integer_b + 2n) - 6) + delta^2
    qprime = 2 * a0 - integer_b
    qconstant = q - qprime * c
    first = (ComplexF64(integer_b + n - 2) - a0 + delta / 2) * (b2n - 2)
    second_factor = b2n - 3
    second = second_factor * q / c
    second_prime = -second_factor * qconstant / c^2
    second_second = 2 * second_factor * qconstant / c^3
    denominator = (a0 + n - 1 + delta / 2) * (b2n - 4)

    previous2_weight = max(abs(first), floatmin(Float64))
    previous1_weight = max(abs(second), abs(second_prime), abs(second_second),
        floatmin(Float64))
    aligned_scale = max(
        previous2.log_scale + log(previous2_weight),
        previous1.log_scale + log(previous1_weight),
    )
    scale2 = exp(previous2.log_scale - aligned_scale)
    scale1 = exp(previous1.log_scale - aligned_scale)
    value_terms = (
        first * previous2.value * scale2,
        second * previous1.value * scale1,
    )
    value_numerator = value_terms[1] + value_terms[2]
    deriv_terms = (
        first * previous2.deriv * scale2,
        second_prime * previous1.value * scale1,
        second * previous1.deriv * scale1,
    )
    deriv_numerator = deriv_terms[1] + deriv_terms[2] + deriv_terms[3]
    second_terms = (
        first * previous2.second * scale2,
        second_second * previous1.value * scale1,
        2 * second_prime * previous1.deriv * scale1,
        second * previous1.second * scale1,
    )
    second_numerator = sum(second_terms)
    condition = max(
        _term_cond(value_numerator, value_terms...),
        _term_cond(deriv_numerator, deriv_terms...),
        _term_cond(second_numerator, second_terms...),
    )
    denominator_abs = abs(denominator)
    if !isfinite(aligned_scale) || !isfinite(denominator_abs) || iszero(denominator_abs)
        direct, direct_error = _direct_h_eval(sequence, n)
        return direct, Inf, direct_error
    end
    denominator_phase = denominator / denominator_abs
    state = _normalize_h(
        value_numerator / denominator_phase,
        deriv_numerator / denominator_phase,
        second_numerator / denominator_phase,
        aligned_scale - log(denominator_abs),
    )
    error = min(Inf, condition * (parent_error + eps(Float64)))
    finite = _finite_complex(state.value) && _finite_complex(state.deriv) &&
        _finite_complex(state.second) && isfinite(state.log_scale)
    if !finite || !isfinite(error) || error > _reanchor_tol(n)
        if n > 1
            _miller_positive!(sequence, n)
            if sequence.miller_status == :ACCEPTED && sequence.miller_nmax >= n
                return sequence.states[_seq_index(n)], condition,
                    sequence.errors[_seq_index(n)]
            end
        end
        direct, direct_error = _direct_h_eval(sequence, n)
        if _finite_complex(direct.value) && _finite_complex(direct.deriv) &&
                _finite_complex(direct.second) && isfinite(direct.log_scale)
            selected_error = _reanchor_error(
                state, error, direct, direct_error)
            return direct, condition, selected_error
        end
    end
    return state, condition, error
end

# h_n from h_(n+1) and h_(n+2); when unstable, direct evaluation.
function _backward_h(sequence::MSTUSequence, n::Int)
    next1 = _h_state(sequence, n + 1)
    next2 = _h_state(sequence, n + 2)
    parent_error = max(
        sequence.errors[_seq_index(n + 1)],
        sequence.errors[_seq_index(n + 2)],
    )
    a0 = sequence.a0
    integer_b = sequence.integer_b
    delta = sequence.delta
    c = sequence.argument

    b2n = ComplexF64(integer_b + 2n) + delta
    r0 = ComplexF64(integer_b)^2 + 4n * (n + 1) +
        ComplexF64(integer_b) * (2 + 4n - c) + 2 * a0 * c
    r = r0 + delta * (2integer_b + 4n + 2) + delta^2
    rprime = 2 * a0 - integer_b
    rconstant = r - rprime * c
    first_factor = b2n + 1
    first = -first_factor * r / c
    first_prime = first_factor * rconstant / c^2
    first_second = -2 * first_factor * rconstant / c^3
    second = (a0 + n + 1 + delta / 2) * b2n
    denominator = (ComplexF64(integer_b + n) - a0 + delta / 2) * (b2n + 2)

    next1_weight = max(abs(first), abs(first_prime), abs(first_second),
        floatmin(Float64))
    next2_weight = max(abs(second), floatmin(Float64))
    aligned_scale = max(
        next1.log_scale + log(next1_weight),
        next2.log_scale + log(next2_weight),
    )
    scale1 = exp(next1.log_scale - aligned_scale)
    scale2 = exp(next2.log_scale - aligned_scale)
    value_terms = (
        first * next1.value * scale1,
        second * next2.value * scale2,
    )
    value_numerator = value_terms[1] + value_terms[2]
    deriv_terms = (
        first_prime * next1.value * scale1,
        first * next1.deriv * scale1,
        second * next2.deriv * scale2,
    )
    deriv_numerator = deriv_terms[1] + deriv_terms[2] + deriv_terms[3]
    second_terms = (
        first_second * next1.value * scale1,
        2 * first_prime * next1.deriv * scale1,
        first * next1.second * scale1,
        second * next2.second * scale2,
    )
    second_numerator = sum(second_terms)
    condition = max(
        _term_cond(value_numerator, value_terms...),
        _term_cond(deriv_numerator, deriv_terms...),
        _term_cond(second_numerator, second_terms...),
    )
    denominator_abs = abs(denominator)
    if !isfinite(aligned_scale) || !isfinite(denominator_abs) || iszero(denominator_abs)
        direct, direct_error = _direct_h_eval(sequence, n)
        return direct, Inf, direct_error
    end
    denominator_phase = denominator / denominator_abs
    state = _normalize_h(
        value_numerator / denominator_phase,
        deriv_numerator / denominator_phase,
        second_numerator / denominator_phase,
        aligned_scale - log(denominator_abs),
    )
    error = min(Inf, condition * (parent_error + eps(Float64)))
    finite = _finite_complex(state.value) && _finite_complex(state.deriv) &&
        _finite_complex(state.second) && isfinite(state.log_scale)
    if !finite || !isfinite(error) || error > _reanchor_tol(n)
        direct, direct_error = _direct_h_eval(sequence, n)
        if _finite_complex(direct.value) && _finite_complex(direct.deriv) &&
                _finite_complex(direct.second) && isfinite(direct.log_scale)
            selected_error = _reanchor_error(
                state, error, direct, direct_error)
            return direct, condition, selected_error
        end
    end
    return state, condition, error
end

@inline function _h_error(candidate::MSTHState, reference::MSTHState)
    return max(
        _h_relerr(candidate, reference, :value),
        _h_relerr(candidate, reference, :deriv),
        _h_relerr(candidate, reference, :second),
    )
end

@inline _reanchor_tol(n::Int) =
    n < 0 ? MST_NEGATIVE_REANCHOR_TOL : MST_STATE_REANCHOR_TOL

# Error of a direct re-anchor: at least the direct error, and at least the smaller of the
# recurrence error and twice the recurrence-direct discrepancy.
function _reanchor_error(state::MSTHState, recurrence_error, direct::MSTHState,
        direct_error)
    discrepancy = _h_error(state, direct)
    guarded_discrepancy = isfinite(discrepancy) ?
        2 * discrepancy + 16eps(Float64) : Inf
    cross_error = isfinite(recurrence_error) ?
        min(recurrence_error, guarded_discrepancy) : guarded_discrepancy
    return max(Float64(direct_error), Float64(cross_error))
end

function _clear_negative!(sequence::MSTUSequence)
    for n in MST_NEGATIVE_SEGMENT_LO:-1
        index = _seq_index(n)
        sequence.ready[index] = false
        sequence.conditions[index] = Inf
        sequence.errors[index] = Inf
    end
    return sequence
end

# Builds the negative-n segment once. An unstable segment is rejected
# (false or a caught exception) and cleared; _backward_h then supplies
# those values.
function _negative_segment!(sequence::MSTUSequence)
    sequence.negative_status == :UNTRIED || return sequence
    sequence.negative_status = :RUNNING
    accepted = try
        _build_negative_segment!(sequence)
    catch
        false
    end
    if !accepted
        _clear_negative!(sequence)
        sequence.negative_error = Inf
        sequence.negative_status = :REJECTED
    end
    return sequence
end

# h_n for MST_NEGATIVE_SEGMENT_LO <= n <= -1 by forward recurrence from two direct
# anchors; accepted when it agrees with the direct h_(-1).
function _build_negative_segment!(sequence::MSTUSequence)
    for n in MST_NEGATIVE_SEGMENT_LO:(MST_NEGATIVE_SEGMENT_LO + 1)
        index = _seq_index(n)
        state, state_error = _direct_h_eval(sequence, n)
        isfinite(state.log_scale) || return false
        sequence.states[index] = state
        sequence.ready[index] = true
        sequence.conditions[index] = 1.0
        sequence.errors[index] = state_error
    end

    for n in (MST_NEGATIVE_SEGMENT_LO + 2):-1
        state, condition, state_error = _forward_h(sequence, n)
        condition <= 2.0 && isfinite(state.log_scale) || return false
        index = _seq_index(n)
        sequence.states[index] = state
        sequence.ready[index] = true
        sequence.conditions[index] = condition
        sequence.errors[index] = state_error
    end
    end_index = _seq_index(-1)
    end_direct, end_error = _direct_h_eval(sequence, -1)
    end_delta = _h_error(sequence.states[end_index], end_direct)
    end_delta <= MST_NEGATIVE_SEGMENT_TOL || return false
    segment_error = 2 * end_delta + 16eps(Float64)
    for n in (MST_NEGATIVE_SEGMENT_LO + 2):-2
        index = _seq_index(n)
        sequence.errors[index] = max(sequence.errors[index], segment_error)
    end
    sequence.states[end_index] = end_direct
    sequence.conditions[end_index] = 1.0
    sequence.errors[end_index] = end_error
    sequence.negative_error = segment_error
    sequence.negative_status = :ACCEPTED
    return true
end

# h_n on demand, cached in the sequence.
function _h_state(sequence::MSTUSequence, n::Integer)
    -MST_N_MAX <= n <= MST_N_MAX || error("MST U sequence index outside supported range")
    index = _seq_index(n)
    sequence.ready[index] && return sequence.states[index]
    MST_NEGATIVE_SEGMENT_LO <= n < 0 &&
        sequence.negative_status == :UNTRIED && _negative_segment!(sequence)
    sequence.ready[index] && return sequence.states[index]
    state, condition, error = n > 1 ? _forward_h(sequence, Int(n)) :
        _backward_h(sequence, Int(n))
    sequence.states[index] = state
    sequence.ready[index] = true
    sequence.conditions[index] = condition
    sequence.errors[index] = error
    return state
end

# U-function sequences of the infinity-out and infinity-in series at MST coordinate x.
function _out_sequence(data::MSTSeriesData, x)
    p = data.params
    epsilon = ComplexF64(p.epsilon)
    delta = 2 * p.nu_offset
    z = epsilon * p.kappa * (1 - ComplexF64(x))
    a0 = ComplexF64(p.l + 1 + p.s) - I * epsilon
    return MSTUSequence(a0, 2 * (p.l + 1), delta, -2I * z)
end

function _in_sequence(data::MSTSeriesData, x)
    p = data.params
    epsilon = ComplexF64(p.epsilon)
    delta = 2 * p.nu_offset
    z = 2I * epsilon * p.kappa * (1 - ComplexF64(x))
    a0 = ComplexF64(p.l + 1 - p.s) + I * epsilon
    return MSTUSequence(a0, 2 * (p.l + 1), delta, z)
end

const MST_TAIL_WINDOW = 3
const MST_TAIL_STREAK = 3
const MST_TAIL_SHARE = 0.25
const MST_SUM_GAMMA = 2eps(Float64)

# Compensated (Neumaier) sums of term triplets with their absolute sums, propagated
# state errors and largest terms.
mutable struct MSTAccum
    value::ComplexF64
    deriv::ComplexF64
    second::ComplexF64
    value_c::ComplexF64
    deriv_c::ComplexF64
    second_c::ComplexF64
    value_abs::Float64
    deriv_abs::Float64
    second_abs::Float64
    value_state::Float64
    deriv_state::Float64
    second_state::Float64
    value_max::Float64
    deriv_max::Float64
    second_max::Float64
    count::Int
end

MSTAccum() = MSTAccum(
    0, 0, 0, 0, 0, 0,
    0.0, 0.0, 0.0,
    0.0, 0.0, 0.0,
    0.0, 0.0, 0.0,
    0,
)

# Convergence monitor of one direction of a series: recent term ratios and the
# geometric tail bound they imply.
mutable struct MSTTail
    direction::Int
    last::NTuple{3,Float64}
    ratios::NTuple{MST_TAIL_WINDOW,NTuple{3,Float64}}
    ratio_count::Int
    bounds::NTuple{3,Float64}
    streak::Int
    index::Int
    status::Symbol
end

MSTTail(direction::Int) = MSTTail(
    direction,
    (NaN, NaN, NaN),
    ntuple(_ -> (NaN, NaN, NaN), MST_TAIL_WINDOW),
    0,
    (Inf, Inf, Inf),
    0,
    0,
    :UNTRIED,
)

# Result of _mst_budget: the sum, its error split into rounding, state and tail parts,
# and the status.
struct MSTBudget
    sum::MSTTriplet
    estimated_relerr::Float64
    nmin::Int
    nmax::Int
    max_term_abs::Float64
    condition::Float64
    roundoff::Float64
    state_error::Float64
    tail_error::Float64
    status::Symbol
    log_scale::Float64
end

@inline _triplet_abs(term::MSTTriplet) =
    (abs(term.value), abs(term.deriv), abs(term.second))

@inline function _neumaier(total, correction, term)
    trial = total + term
    next_c = correction +
        (abs(total) >= abs(term) ? (total - trial) + term : (term - trial) + total)
    return trial, next_c
end

@inline _weighted_error(magnitude, error) =
    iszero(magnitude) ? 0.0 :
    isfinite(error) ? magnitude * max(0.0, error) : Inf

function _accum_add!(acc::MSTAccum, term::MSTTriplet, state_error::Float64)
    acc.value, acc.value_c = _neumaier(acc.value, acc.value_c, term.value)
    acc.deriv, acc.deriv_c = _neumaier(acc.deriv, acc.deriv_c, term.deriv)
    acc.second, acc.second_c = _neumaier(acc.second, acc.second_c, term.second)
    value_abs, deriv_abs, second_abs = _triplet_abs(term)
    acc.value_abs += value_abs
    acc.deriv_abs += deriv_abs
    acc.second_abs += second_abs
    acc.value_state += _weighted_error(value_abs, state_error)
    acc.deriv_state += _weighted_error(deriv_abs, state_error)
    acc.second_state += _weighted_error(second_abs, state_error)
    acc.value_max = max(acc.value_max, value_abs)
    acc.deriv_max = max(acc.deriv_max, deriv_abs)
    acc.second_max = max(acc.second_max, second_abs)
    acc.count += 1
    return acc
end

@inline _accum_sum(acc::MSTAccum) = MSTTriplet(
    acc.value + acc.value_c,
    acc.deriv + acc.deriv_c,
    acc.second + acc.second_c,
)

@inline function _tail_ratio(current, previous)
    iszero(previous) && return iszero(current) ? 0.0 : Inf
    return current / previous
end

function _tail_observe!(tail::MSTTail, term::MSTTriplet, n::Int)
    current = _triplet_abs(term)
    if all(isfinite, tail.last)
        ratio = ntuple(j -> _tail_ratio(current[j], tail.last[j]), 3)
        if tail.ratio_count < MST_TAIL_WINDOW
            tail.ratio_count += 1
            tail.ratios = Base.setindex(tail.ratios, ratio, tail.ratio_count)
        else
            tail.ratios = (tail.ratios[2], tail.ratios[3], ratio)
        end
    end
    tail.last = current
    tail.index = n
    if tail.ratio_count < MST_TAIL_WINDOW
        tail.bounds = (Inf, Inf, Inf)
        return tail
    end
    tail.bounds = ntuple(3) do j
        ratio = max(tail.ratios[1][j], tail.ratios[2][j], tail.ratios[3][j])
        isfinite(ratio) && ratio < 1 ? current[j] * ratio / (1 - ratio) : Inf
    end
    return tail
end

function _tail_small(acc::MSTAccum, tail::MSTTail, reltol, abstol)
    total = _accum_sum(acc)
    scales = (
        max(abs(total.value), acc.value_max),
        max(abs(total.deriv), acc.deriv_max),
        max(abs(total.second), acc.second_max),
    )
    return all(j -> tail.bounds[j] <=
        MST_TAIL_SHARE * (abstol + reltol * scales[j]), 1:3)
end

function _budget_fail(status::Symbol, log_scale::Float64=0.0)
    nan = ComplexF64(NaN + NaN * im)
    return MSTBudget(
        MSTTriplet(nan, nan, nan),
        Inf, 0, 0, Inf, Inf, Inf, Inf, Inf, status, log_scale,
    )
end

function _budget_term(termfun, errorfun, n)
    term = try
        termfun(n)
    catch
        return nothing
    end
    _finite_complex(term.value) && _finite_complex(term.deriv) &&
        _finite_complex(term.second) || return nothing
    state_error = try
        Float64(errorfun(n))
    catch
        Inf
    end
    return term, state_error
end

function _extend_tail!(acc, tail, termfun, errorfun, start, nmax, reltol, abstol)
    n = tail.direction * start
    while abs(n) <= nmax
        evaluated = _budget_term(termfun, errorfun, n)
        if evaluated === nothing
            tail.status = :NONFINITE
            return tail
        end
        term, state_error = evaluated
        _accum_add!(acc, term, state_error)
        _tail_observe!(tail, term, n)
        if _tail_small(acc, tail, reltol, abstol)
            tail.streak += 1
            if tail.streak >= MST_TAIL_STREAK
                tail.status = :OK
                return tail
            end
        else
            tail.streak = 0
        end
        n += tail.direction
    end
    tail.status = :MAX_INDEX
    return tail
end

# Sums termfun(n) outward from n = 0 in both directions until each tail bound stays below
# MST_TAIL_SHARE of the target for MST_TAIL_STREAK terms. The error is rounding plus
# propagated state error (errorfun) plus tail bound, compared with budget_tol.
function _mst_budget(
    termfun,
    errorfun;
    reltol::Float64=MST_REL_ERROR,
    budget_tol::Float64=reltol,
    abstol::Float64=0.0,
    ninit::Int=MST_N_INIT,
    nmax::Int=MST_N_MAX,
    log_scale::Float64=0.0,
)
    acc = MSTAccum()
    positive = MSTTail(1)
    negative = MSTTail(-1)
    center = _budget_term(termfun, errorfun, 0)
    center === nothing && return _budget_fail(:NONFINITE_CENTER, log_scale)
    _accum_add!(acc, center...)

    for n in 1:(ninit - 1)
        pos = _budget_term(termfun, errorfun, n)
        neg = _budget_term(termfun, errorfun, -n)
        (pos === nothing || neg === nothing) &&
            return _budget_fail(:NONFINITE_CORE, log_scale)
        _accum_add!(acc, pos...)
        _tail_observe!(positive, pos[1], n)
        _accum_add!(acc, neg...)
        _tail_observe!(negative, neg[1], -n)
    end

    _extend_tail!(acc, positive, termfun, errorfun, ninit, nmax, reltol, abstol)
    _extend_tail!(acc, negative, termfun, errorfun, ninit, nmax, reltol, abstol)
    total = _accum_sum(acc)
    totals = _triplet_abs(total)
    abs_sums = (acc.value_abs, acc.deriv_abs, acc.second_abs)
    state_abs = (acc.value_state, acc.deriv_state, acc.second_state)
    tail_abs = ntuple(j -> positive.bounds[j] + negative.bounds[j], 3)
    round_abs = ntuple(j -> MST_SUM_GAMMA * abs_sums[j], 3)
    targets = ntuple(j -> abstol + budget_tol * totals[j], 3)
    errors = ntuple(j -> round_abs[j] + state_abs[j] + tail_abs[j], 3)
    denominators = ntuple(j -> max(totals[j], floatmin(Float64)), 3)
    relerrs = ntuple(j -> errors[j] / denominators[j], 3)
    conditions = ntuple(j -> abs_sums[j] / denominators[j], 3)
    round_rel = maximum(ntuple(j -> round_abs[j] / denominators[j], 3))
    state_rel = maximum(ntuple(j -> state_abs[j] / denominators[j], 3))
    tail_rel = maximum(ntuple(j -> tail_abs[j] / denominators[j], 3))
    tails_ok = positive.status == :OK && negative.status == :OK
    budget_ok = all(j -> errors[j] <= targets[j], 1:3)
    status =
        !tails_ok ? :TAIL_NOT_CONVERGED :
        budget_ok ? :OK :
        round_rel >= max(state_rel, tail_rel) ? :ROUNDING_LIMIT :
        state_rel >= tail_rel ? :STATE_LIMIT : :TAIL_LIMIT
    return MSTBudget(
        total,
        maximum(relerrs),
        negative.index,
        positive.index,
        max(acc.value_max, acc.deriv_max, acc.second_max),
        maximum(conditions),
        round_rel,
        state_rel,
        tail_rel,
        status,
        log_scale,
    )
end

# Relative error of the n-th term: coefficient error plus U-state error.
function _term_error(data::MSTSeriesData, sequence::MSTUSequence, n::Integer)
    index = _seq_index(n)
    coefficient_error = _coefficient_error(data, n)
    state_error = sequence.errors[index]
    return min(Inf, coefficient_error + state_error + eps(Float64))
end

@inline function _add_logterm(total::ComplexF64, scale::Float64, logterm)
    value = ComplexF64(logterm)
    magnitude_log = real(value)
    if magnitude_log == -Inf
        return total, scale, 0.0
    end
    isfinite(magnitude_log) && isfinite(imag(value)) ||
        return ComplexF64(NaN, NaN), NaN, NaN
    phase = cis(imag(value))
    if magnitude_log > scale
        return total * exp(scale - magnitude_log) + phase, magnitude_log, 1.0
    end
    scaled_magnitude = exp(magnitude_log - scale)
    return total + scaled_magnitude * phase, scale, scaled_magnitude
end

# Sum of exp(logtermfun(n)) over n = 0, +-1, ... until the new terms fall below reltol,
# in log-norm form.
function _sum_symmetric_logs(logtermfun;
    reltol::Float64=MST_REL_ERROR,
    ninit::Int=MST_N_INIT,
    nmax::Int=MST_N_MAX,
)
    total = ComplexF64(0)
    scale = -Inf
    for n in 0:(ninit - 1)
        total, scale, _ = _add_logterm(total, scale, logtermfun(n))
        if n > 0
            total, scale, _ = _add_logterm(total, scale, logtermfun(-n))
        end
    end

    err_pos = Inf
    err_neg = Inf
    n_pos = ninit
    n_neg = ninit
    while n_pos <= nmax && err_pos > reltol
        total, scale, term_magnitude =
            _add_logterm(total, scale, logtermfun(n_pos))
        err_pos = term_magnitude / max(abs(total), floatmin(Float64))
        n_pos += 1
    end
    while n_neg <= nmax && err_neg > reltol
        total, scale, term_magnitude =
            _add_logterm(total, scale, logtermfun(-n_neg))
        err_neg = term_magnitude / max(abs(total), floatmin(Float64))
        n_neg += 1
    end

    norm = abs(total)
    isfinite(scale) && isfinite(norm) && !iszero(norm) ||
        return LogNormSum(ComplexF64(NaN, NaN), NaN, Inf,
            -(n_neg - 1), n_pos - 1)
    roundoff = eps(Float64) / max(norm, floatmin(Float64))
    return LogNormSum(
        ComplexF64(total / norm),
        Float64(scale + log(norm)),
        max(roundoff, err_pos, err_neg),
        -(n_neg - 1),
        n_pos - 1,
    )
end

# Budgeted sum of the infinity-out (U-function) series; throws MSTCertificateError
# unless the budget is :OK.
function _sum_type2(data::MSTSeriesData, x; reltol::Float64=MST_REL_ERROR,
        sequence=nothing)
    sequence = sequence === nothing ? _out_sequence(data, x) : sequence
    log_scale = _out_term_log_scale(data, x, sequence)
    budget = _mst_budget(
        n -> _out_term_triplet(data, n, x, sequence, log_scale),
        n -> _term_error(data, sequence, n);
        reltol=reltol,
        budget_tol=MST_RESIDUAL_TOL,
        log_scale=log_scale,
    )
    budget.status == :OK || throw(MSTCertificateError(
        budget.status == :TAIL_NOT_CONVERGED ?
            :infinity_out_tail : :infinity_out_sum,
        "MST infinity-out sum rejected: status=$(budget.status), " *
        "error=$(budget.estimated_relerr), condition=$(budget.condition)"))
    return MSTSum(
        budget.sum.value,
        budget.sum.deriv,
        budget.estimated_relerr,
        budget.nmin,
        budget.nmax,
        budget.max_term_abs,
        budget.log_scale,
    )
end

# Prefactor of the infinity-out series and its x-derivative.
function _up_prefactor(p::MSTParams, r::Number)
    z = ComplexF64(mst_z(p, r))
    epsv = ComplexF64(p.epsilon)
    kappa = ComplexF64(p.kappa)
    tau = ComplexF64(p.tau)
    s = ComplexF64(p.s)
    nu = _nu_shift(p)
    value = 2^nu * exp(-pi * epsv) * exp(-I * pi * (nu + 1 + s))
    value *= exp(I * z) * z^(nu + I * (epsv + tau) / 2) *
        (z - epsv * kappa)^(-s - I * (epsv + tau) / 2)
    return value
end

function _dup_prefactor_dx(p::MSTParams, r::Number)
    z = ComplexF64(mst_z(p, r))
    epsv = ComplexF64(p.epsilon)
    kappa = ComplexF64(p.kappa)
    tau = ComplexF64(p.tau)
    s = ComplexF64(p.s)
    nu = _nu_shift(p)
    value = _up_prefactor(p, r)
    value *= epsv * kappa * I *
        (epsv^2 * kappa + 2I * z * (nu - s + I * z) +
         epsv * kappa * (2 * z - 2I * nu + tau))
    return value / (2 * z * (z - epsv * kappa))
end

# P and dP/dx of the infinity-out solution at r, with relative error and log scale.
function _mst_out_p_pair(data::MSTSeriesData, r, sequence=nothing, budget=nothing)
    p = data.params
    x = mst_x(p, r)
    sum = budget === nothing ? _sum_type2(data, x; sequence=sequence) : budget
    value = budget === nothing ? sum.value : sum.sum.value
    derivative = budget === nothing ? sum.derivative : sum.sum.deriv
    estimated_relerr = sum.estimated_relerr
    log_scale = sum.log_scale
    tt = _teukolsky_transformation_module()
    Tfun, _, TpTfun = getfield(tt, :Tx)(p.s, p.epsilon, p.kappa, p.tau)
    Tval = Tfun(x)
    Uval = _up_prefactor(p, r)
    Ux = _dup_prefactor_dx(p, r)
    G = Uval / Tval
    Gx = (Ux - Uval * TpTfun(x)) / Tval
    P = G * value
    Px = Gx * value + G * derivative
    return P, Px, estimated_relerr, log_scale
end

# Budgeted sum of the infinity-in series; throws MSTCertificateError unless :OK.
function _sum_true_in(data::MSTSeriesData, x; reltol::Float64=MST_REL_ERROR,
        sequence=nothing)
    sequence = sequence === nothing ? _in_sequence(data, x) : sequence
    budget = _mst_budget(
        n -> _in_term_triplet(data, n, x, sequence),
        n -> _term_error(data, sequence, n);
        reltol=reltol,
        budget_tol=MST_RESIDUAL_TOL,
    )
    budget.status == :OK || throw(MSTCertificateError(
        budget.status == :TAIL_NOT_CONVERGED ?
            :infinity_in_tail : :infinity_in_sum,
        "MST infinity-in sum rejected: status=$(budget.status), " *
        "error=$(budget.estimated_relerr), condition=$(budget.condition)"))
    return MSTSum(
        budget.sum.value,
        budget.sum.deriv,
        budget.estimated_relerr,
        budget.nmin,
        budget.nmax,
        budget.max_term_abs,
        budget.log_scale,
    )
end

# P and dP/dx of the infinity-in solution at r, with relative error and log scale.
function _mst_in_p_pair(data::MSTSeriesData, r, sequence=nothing, budget=nothing)
    p = data.params
    x = mst_x(p, r)
    series = budget === nothing ? _sum_true_in(data, x; sequence=sequence) : budget
    value = budget === nothing ? series.value : series.sum.value
    derivative = budget === nothing ? series.derivative : series.sum.deriv
    estimated_relerr = series.estimated_relerr
    exponent = -I * p.epsilon + p.s + I * p.tau - 1
    one_minus_x = 1 - ComplexF64(x)
    pref = exp(exponent * log(one_minus_x))
    pref_x = -exponent * pref / one_minus_x
    P = pref * value
    Px = pref_x * value + pref * derivative
    return P, Px, estimated_relerr, series.log_scale
end

# Closure r -> coefficients (A0, A1, B0, B1) with X = A0 P + A1 dP/dx and
# dX/dr* = B0 P + B1 dP/dx.
function _p_converter(params)
    tt = _teukolsky_transformation_module()
    return getfield(tt, :P_to_GSN_coefficients_from_matrix)(
        _teukolsky_from_gsn_matrix(params),
        params.s,
        params.m,
        params.a,
        params.omega,
        params.lambda,
    )
end

# GSN X and dX/dx from (P, dP/dx) and the converter coefficients.
function _p_to_gsn_dx(
    params::DirectGSNParameters,
    P,
    Px,
    direct_x,
    transform;
    radius=nothing,
)
    A0, A1, B0, B1 = transform
    X = A0 * P + A1 * Px
    dXdrstar = B0 * P + B1 * Px
    dxdrstar = radius === nothing ? _direct_dx_drstar(params, direct_x) :
        2params.kappa * direct_x / (radius^2 + params.a^2)
    iszero(dxdrstar) && error("direct GSN MST cannot convert dX/drstar at dx/drstar = 0.")
    return ComplexF64(X), ComplexF64(dXdrstar / dxdrstar)
end

function _p_to_gsn_dx(
    coefficients::DirectCoefficientSet,
    P,
    Px,
    direct_x,
    transform;
    radius=nothing,
)
    return _p_to_gsn_dx(
        coefficients.params, P, Px, direct_x, transform; radius)
end

# Normalization sum of a_n over n, in log-norm form.
function _incoming_raw_factor(data::MSTSeriesData)
    return _sum_symmetric_logs(n -> log_series_coefficient!(data, n))
end

# Compensated sum of exp(logtermfun(n)) over |n| <= nmax in log-norm form, with a
# rounding bound.
function _nia_full_logsum(logtermfun, nmax::Int)
    nmax >= 0 || throw(ArgumentError(
        "MST NIA full-term sum requires nmax >= 0."))
    indices = collect(-nmax:nmax)
    logs = ComplexF64[logtermfun(n) for n in indices]
    scale = maximum(real, logs)
    isfinite(scale) ||
        throw(MSTCertificateError(
            :nia_sum_scale,
            "nonfinite MST NIA full-term sum scale.",
        ))
    total = ComplexF64(0)
    correction = ComplexF64(0)
    absolute_sum = 0.0
    for value in logs
        term = ComplexF64(exp(value - scale))
        y = term - correction
        updated = total + y
        correction = (updated - total) - y
        total = updated
        absolute_sum += abs(term)
    end
    norm = abs(total)
    isfinite(norm) && !iszero(norm) ||
        throw(MSTCertificateError(
            :nia_sum_cancellation,
            "zero or nonfinite MST NIA full-term sum.",
        ))
    count = length(logs)
    gamma = count * eps(Float64) /
        max(1 - count * eps(Float64), eps(Float64))
    return LogNormSum(
        ComplexF64(total / norm),
        Float64(scale + log(norm)),
        Float64(gamma * absolute_sum / norm),
        -nmax,
        nmax,
    )
end

# Double-double l + n + shift + offset and 2(l + n) + shift + 2 offset.
@inline _dd_factor(
    p::MSTParams,
    offset::DDComplex,
    n::Int,
    shift::Int=0,
) =
    dc_add(
        DDComplex(Float64(p.l + n + shift)),
        offset,
    )

@inline _dd_twice(
    p::MSTParams,
    offset::DDComplex,
    n::Int,
    shift::Int=0,
) =
    dc_add(
        DDComplex(Float64(2 * (p.l + n) + shift)),
        dc_add(offset, offset),
    )

@inline _nia_dd_factor(p::MSTParams, n::Int, shift::Int=0) =
    _dd_factor(p, DDComplex(p.nu_offset), n, shift)

@inline _nia_dd_twice_factor(p::MSTParams, n::Int, shift::Int=0) =
    _dd_twice(p, DDComplex(p.nu_offset), n, shift)

# Double-double alpha_n, beta_n and gamma_n of the MST recurrence at nu = l + offset.
function _dd_alpha(n::Int, p::MSTParams, offset::DDComplex)
    epsc = DDComplex(p.epsilon)
    kappa = DDComplex(p.kappa)
    tau = DDComplex(p.tau)
    npnu1 = _dd_factor(p, offset, n, 1)
    spin_factor = _dd_factor(p, offset, n, 1 + p.s)
    ieps = dc_imul(epsc)
    numerator = dc_mul(
        dc_mul(
            dc_mul(dc_imul(dc_mul(epsc, kappa)),
                dc_add(spin_factor, ieps)),
            dc_sub(spin_factor, ieps),
        ),
        dc_add(npnu1, dc_imul(tau)),
    )
    denominator = dc_mul(
        npnu1,
        _dd_twice(p, offset, n, 3),
    )
    return dc_div(numerator, denominator)
end

function _dd_beta(n::Int, p::MSTParams, offset::DDComplex)
    epsc = DDComplex(p.epsilon)
    tau = DDComplex(p.tau)
    product = dc_mul(
        _dd_factor(p, offset, n),
        _dd_factor(p, offset, n, 1),
    )
    eps2 = dc_mul(epsc, epsc)
    eps_kappa_tau = dc_mul(
        dc_mul(epsc, DDComplex(p.kappa)),
        tau,
    )
    degree = p.l + n
    angular = dc_add(
        dc_sub(DDComplex(degree*(degree+1)-p.s*(p.s+1)), DDComplex(p.lambda)),
        dc_mul(offset, dc_add(DDComplex(2degree+1), offset)),
    )
    base = dc_add(angular, dc_add(eps2, eps_kappa_tau))
    correction = dc_div(
        dc_mul(
            eps_kappa_tau,
            dc_add(DDComplex(p.s * p.s), eps2),
        ),
        product,
    )
    return dc_add(base, correction)
end

function _dd_gamma(n::Int, p::MSTParams, offset::DDComplex)
    epsc = DDComplex(p.epsilon)
    kappa = DDComplex(p.kappa)
    tau = DDComplex(p.tau)
    npnu = _dd_factor(p, offset, n)
    spin_factor = _dd_factor(p, offset, n, -p.s)
    ieps = dc_imul(epsc)
    numerator = dc_mul(
        dc_mul(
            dc_mul(dc_neg(dc_imul(dc_mul(epsc, kappa))),
                dc_add(spin_factor, ieps)),
            dc_sub(spin_factor, ieps),
        ),
        dc_sub(npnu, dc_imul(tau)),
    )
    denominator = dc_mul(
        npnu,
        _dd_twice(p, offset, n, -1),
    )
    return dc_div(numerator, denominator)
end

@inline _nia_alpha_dd(n::Int, p::MSTParams) =
    _dd_alpha(n, p, DDComplex(p.nu_offset))

@inline _nia_beta_dd(n::Int, p::MSTParams) =
    _dd_beta(n, p, DDComplex(p.nu_offset))

@inline _nia_gamma_dd(n::Int, p::MSTParams) =
    _dd_gamma(n, p, DDComplex(p.nu_offset))

# Sum over |n| <= nmax of a_n (-1)^n (nu+1+s-i epsilon)_n / (nu+1-s+i epsilon)_n.
function _nia_aminus_sum(data::MSTSeriesData, nmax::Int)
    p = data.params
    return _nia_full_logsum(nmax) do n
        log_series_coefficient!(data, n) +
            I * pi * n +
            log_pochhammer(_nu_shift(p, 1 + p.s) - I * p.epsilon, n) -
            log_pochhammer(_nu_shift(p, 1 - p.s) + I * p.epsilon, n)
    end
end

# MSTSeriesData with a_n for |n| <= MST_N_MAX from double-double backward recurrences
# started at +-cutoff.
function _nia_series_data(params::MSTParams, cutoff::Int)
    cutoff >= MST_N_MAX ||
        throw(ArgumentError(
            "MST NIA recurrence cutoff must be at least $MST_N_MAX."))
    positive = Vector{ComplexF64}(undef, MST_N_MAX)
    negative = Vector{ComplexF64}(undef, MST_N_MAX)

    ratio = DDComplex(0)
    for n in cutoff:-1:1
        denominator = dc_add(
            _nia_beta_dd(n, params),
            dc_mul(_nia_alpha_dd(n, params), ratio),
        )
        denominator_value = dc_value(denominator)
        _finite_complex(denominator_value) && !iszero(denominator_value) ||
            throw(MSTCertificateError(
                :nia_positive_tail,
                "nonfinite MST NIA positive-tail denominator at n=$n.",
            ))
        ratio = dc_div(
            dc_neg(_nia_gamma_dd(n, params)),
            denominator,
        )
        n <= MST_N_MAX && (positive[n] = ComplexF64(dc_value(ratio)))
    end

    ratio = DDComplex(0)
    for n in (-cutoff):-1
        denominator = dc_add(
            _nia_beta_dd(n, params),
            dc_mul(_nia_gamma_dd(n, params), ratio),
        )
        denominator_value = dc_value(denominator)
        _finite_complex(denominator_value) && !iszero(denominator_value) ||
            throw(MSTCertificateError(
                :nia_negative_tail,
                "nonfinite MST NIA negative-tail denominator at n=$n.",
            ))
        ratio = dc_div(
            dc_neg(_nia_alpha_dd(n, params)),
            denominator,
        )
        n >= -MST_N_MAX && (negative[-n] = ComplexF64(dc_value(ratio)))
    end

    data = MSTSeriesData(params)
    for n in 1:MST_N_MAX
        value = positive[n]
        _finite_complex(value) && !iszero(value) ||
            throw(MSTCertificateError(
                :nia_positive_coefficient,
                "nonfinite MST NIA positive coefficient ratio at n=$n.",
            ))
        data.log_coeffs[n] = data.log_coeffs[n - 1] + log(value)
        data.cf_relerrs[n] = eps(Float64)
        data.cf_ok[n] = true
        data.coeff_errors[n] = data.coeff_errors[n - 1] + 2eps(Float64)
    end
    for n in 1:MST_N_MAX
        value = negative[n]
        _finite_complex(value) && !iszero(value) ||
            throw(MSTCertificateError(
                :nia_negative_coefficient,
                "nonfinite MST NIA negative coefficient ratio at n=$(-n).",
            ))
        data.log_coeffs[-n] = data.log_coeffs[-n + 1] + log(value)
        data.cf_relerrs[-n] = eps(Float64)
        data.cf_ok[-n] = true
        data.coeff_errors[-n] = data.coeff_errors[-n + 1] + 2eps(Float64)
    end
    return data
end

@inline function _lognorm_value(value::LogNormSum)
    _finite_complex(value.phase) && isfinite(value.logabs) &&
        !iszero(value.phase) ||
        throw(MSTCertificateError(
            :nia_sum,
            "nonfinite MST NIA normalization sum.",
        ))
    return ComplexF64(value.logabs + I * angle(value.phase))
end

# A series sum in log form: its condition sum|t|/|sum t| and a running bound
# on the relative error of exp(logvalue).
struct MSTAmplitudeSum
    logvalue::ComplexF64
    condition::Float64
    error::Float64
end

# Representation shifts nu + k and acceptance thresholds of the amplitude certificates.
const MST_AMP_SHIFTS = (-4, -3, 2, 3)
const MST_AMP_NMAX = 120
const MST_AMP_CHECK_NMAX = 96
const MST_AMP_SPREAD_MAX = 1.0e-8
const MST_AMP_NEAREST_MAX = 1.0e-9
const MST_AMP_TRUNCATION_MAX = 1.0e-10
const MST_AMP_CONDITION_MAX = 1.0e7
const MST_DD_NU_RESIDUAL_MAX = 1.0e-24
const MST_DD_CONDITION_MAX = 1.0e12

# Double-double ratios a_n / a_(n-1) (positive) and a_(-n) / a_(-n+1) (negative) for
# n = 1..nmax, from backward recurrences started at +-cutoff.
function _dd_ratios(
    p::MSTParams,
    offset::DDComplex;
    cutoff::Int=MST_NIA_TAIL,
    nmax::Int=MST_AMP_NMAX,
)
    cutoff >= nmax || throw(ArgumentError(
        "DD MST recurrence cutoff must be at least its sum order."))
    positive = Vector{DDComplex}(undef, nmax)
    negative = Vector{DDComplex}(undef, nmax)

    ratio = DDComplex(0)
    for n in cutoff:-1:1
        denominator = dc_add(
            _dd_beta(n, p, offset),
            dc_mul(_dd_alpha(n, p, offset), ratio),
        )
        value = dc_value(denominator)
        _finite_complex(value) && !iszero(value) ||
            throw(MSTCertificateError(
                :dd_nu_positive,
                "nonfinite DD MST positive-tail denominator at n=$n.",
            ))
        ratio = dc_div(
            dc_neg(_dd_gamma(n, p, offset)),
            denominator,
        )
        if n <= nmax
            value = dc_value(ratio)
            _finite_complex(value) && !iszero(value) ||
                throw(MSTCertificateError(
                    :dd_nu_positive,
                    "nonfinite DD MST positive-tail ratio at n=$n.",
                ))
            positive[n] = ratio
        end
    end

    ratio = DDComplex(0)
    for n in (-cutoff):-1
        denominator = dc_add(
            _dd_beta(n, p, offset),
            dc_mul(_dd_gamma(n, p, offset), ratio),
        )
        value = dc_value(denominator)
        _finite_complex(value) && !iszero(value) ||
            throw(MSTCertificateError(
                :dd_nu_negative,
                "nonfinite DD MST negative-tail denominator at n=$n.",
            ))
        ratio = dc_div(
            dc_neg(_dd_alpha(n, p, offset)),
            denominator,
        )
        if n >= -nmax
            value = dc_value(ratio)
            _finite_complex(value) && !iszero(value) ||
                throw(MSTCertificateError(
                    :dd_nu_negative,
                    "nonfinite DD MST negative-tail ratio at n=$n.",
                ))
            negative[-n] = ratio
        end
    end
    return positive, negative
end

# Double-double MST equation beta_0 + alpha_0 R_1 + gamma_0 L_(-1) at nu = l + offset.
function _dd_equation(
    p::MSTParams,
    offset::DDComplex;
    cutoff::Int=MST_NIA_TAIL,
)
    positive, negative = _dd_ratios(
        p, offset; cutoff, nmax=1)
    return dc_add(
        _dd_beta(0, p, offset),
        dc_add(
            dc_mul(_dd_alpha(0, p, offset), positive[1]),
            dc_mul(_dd_gamma(0, p, offset), negative[1]),
        ),
    )
end

# Newton refinement of nu - l in double-double; throws unless the scaled residual
# reaches MST_DD_NU_RESIDUAL_MAX.
function _dd_refine(
    p::MSTParams;
    cutoff::Int=MST_NIA_TAIL,
    maxiter::Int=8,
)
    current = DDComplex(p.nu_offset)
    best = current
    best_value = _dd_equation(p, current; cutoff)
    equation_scale = max(
        abs(p.lambda), abs(p.l * (p.l + 1)), 1.0)
    best_residual = abs(dc_value(best_value)) / equation_scale
    step_scale = max(1.0, abs(p.nu_offset))
    h = exp2(-18) * step_scale

    for _ in 1:maxiter
        value = _dd_equation(p, current; cutoff)
        plus = _dd_equation(
            p, dc_add(current, DDComplex(h)); cutoff)
        minus = _dd_equation(
            p, dc_sub(current, DDComplex(h)); cutoff)
        derivative = dc_div(
            dc_sub(plus, minus),
            DDComplex(2h),
        )
        derivative_value = dc_value(derivative)
        _finite_complex(derivative_value) &&
            !iszero(derivative_value) ||
            throw(MSTCertificateError(
                :dd_nu_derivative,
                "nonfinite DD MST nu derivative.",
            ))

        step = dc_div(value, derivative)
        trial = current
        trial_value = value
        accepted = false
        for damping in 0:8
            candidate = dc_sub(
                current,
                dc_mul(step, DDComplex(exp2(-damping))),
            )
            candidate_value = _dd_equation(p, candidate; cutoff)
            if abs(dc_value(candidate_value)) <
                    abs(dc_value(trial_value))
                trial = candidate
                trial_value = candidate_value
                accepted = true
                break
            end
        end
        accepted || break
        current = trial
        residual = abs(dc_value(trial_value)) / equation_scale
        if residual < best_residual
            best = current
            best_value = trial_value
            best_residual = residual
        end
        abs(dc_value(step)) <= 1.0e-25 * step_scale && break
    end

    best_residual <= MST_DD_NU_RESIDUAL_MAX ||
        throw(MSTCertificateError(
            :dd_nu_residual,
            "DD MST nu refinement rejected: residual=$best_residual.",
        ))
    return best, Float64(best_residual)
end

function _dd_sum_result(total::DDComplex, absolute_sum)
    value = dc_value(total)
    _finite_complex(value) && !iszero(value) ||
        throw(MSTCertificateError(
            :dd_amplitude_sum,
            "nonfinite DD MST amplitude sum.",
        ))
    condition = Float64(absolute_sum / abs(value))
    isfinite(condition) || throw(MSTCertificateError(
        :dd_amplitude_sum,
        "nonfinite DD MST amplitude condition.",
    ))
    # Double-double terms: the sum rounds at eps^2 per term (condition
    # times that) before the final rounding to Float64.
    return MSTAmplitudeSum(
        ComplexF64(log(abs(value)) + I * angle(value)),
        condition,
        eps(Float64) * (1 + eps(Float64) * condition),
    )
end

# Double-double sums of base times cumulative products of ratios and weights, over
# |n| <= nmax and |n| <= check_nmax.
function _dd_sum_pair(
    base,
    positive,
    negative,
    positive_weight,
    negative_weight,
    nmax::Int,
    check_nmax::Int,
)
    0 <= check_nmax <= nmax ||
        throw(ArgumentError("invalid DD MST amplitude check order."))
    positive_terms = Vector{DDComplex}(undef, nmax)
    negative_terms = Vector{DDComplex}(undef, nmax)

    term = DDComplex(base)
    for n in 0:(nmax - 1)
        term = dc_mul(
            term,
            dc_mul(positive[n + 1], positive_weight(n)),
        )
        positive_terms[n + 1] = term
    end
    term = DDComplex(base)
    for k in 1:nmax
        n = -k
        term = dc_mul(
            term,
            dc_mul(negative[k], negative_weight(n)),
        )
        negative_terms[k] = term
    end

    function sum_at(limit)
        total = DDComplex(base)
        absolute_sum = abs(base)
        for k in 1:limit
            total = dc_add(total, positive_terms[k])
            total = dc_add(total, negative_terms[k])
            absolute_sum += abs(dc_value(positive_terms[k])) +
                abs(dc_value(negative_terms[k]))
        end
        return _dd_sum_result(total, absolute_sum)
    end
    return sum_at(nmax), sum_at(check_nmax)
end

# Double-double amplitude sums (coefficient sum, A_- sum, K_nu sums, D sum) at full and
# check order.
function _dd_amp_sums(
    data::MSTSeriesData,
    offset::DDComplex;
    nmax::Int=MST_AMP_NMAX,
    check_nmax::Int=MST_AMP_CHECK_NMAX,
)
    p = data.params
    positive, negative = _dd_ratios(
        p, offset; nmax)
    unit_weight(_) = DDComplex(1)
    fsum = _dd_sum_pair(
        1.0 + 0.0im,
        positive,
        negative,
        unit_weight,
        unit_weight,
        nmax,
        check_nmax,
    )

    ieps = dc_imul(DDComplex(p.epsilon))
    aminus_positive(n) = dc_neg(dc_div(
        dc_sub(_dd_factor(p, offset, n, 1 + p.s), ieps),
        dc_add(_dd_factor(p, offset, n, 1 - p.s), ieps),
    ))
    aminus_negative(n) = dc_neg(dc_div(
        dc_add(_dd_factor(p, offset, n, 1 - p.s), ieps),
        dc_sub(_dd_factor(p, offset, n, 1 + p.s), ieps),
    ))
    asum = _dd_sum_pair(
        1.0 + 0.0im,
        positive,
        negative,
        aminus_positive,
        aminus_negative,
        nmax,
        check_nmax,
    )

    epsc = ComplexF64(p.epsilon)
    tau = ComplexF64(p.tau)
    nu = _nu_shift(p)
    up_a = 1 + p.s + I * epsc + nu
    up_b = 1 + 2nu
    up_c = 1 + nu + I * tau
    up_d = 1 - p.s - I * epsc + nu
    up_e = 1 + nu - I * tau
    up_base = exp(
        _amp_lgamma(up_a) +
        _amp_lgamma(up_b) +
        _amp_lgamma(up_c) -
        _amp_lgamma(up_d) -
        _amp_lgamma(up_e)
    )
    up_term = DDComplex(up_base)
    up_total = up_term
    up_absolute = abs(up_base)
    up_check = nothing
    up_check_absolute = NaN
    for n in 0:(nmax - 1)
        up_a_dd = dc_add(
            _dd_factor(p, offset, n, 1 + p.s), ieps)
        up_b_dd = _dd_twice(p, offset, 0, n + 1)
        up_c_dd = dc_add(
            _dd_factor(p, offset, n, 1),
            dc_imul(DDComplex(p.tau)),
        )
        up_d_dd = dc_sub(
            _dd_factor(p, offset, n, 1 - p.s), ieps)
        up_e_dd = dc_sub(
            _dd_factor(p, offset, n, 1),
            dc_imul(DDComplex(p.tau)),
        )
        rational = dc_neg(dc_div(
            dc_mul(dc_mul(up_a_dd, up_b_dd), up_c_dd),
            dc_mul(
                dc_mul(DDComplex(n + 1), up_d_dd),
                up_e_dd,
            ),
        ))
        up_term = dc_mul(
            up_term,
            dc_mul(positive[n + 1], rational),
        )
        up_total = dc_add(up_total, up_term)
        up_absolute += abs(dc_value(up_term))
        if n + 1 == check_nmax
            up_check = up_total
            up_check_absolute = up_absolute
        end
    end

    down_term = DDComplex(1)
    down_total = down_term
    down_absolute = 1.0
    down_check = nothing
    down_check_absolute = NaN
    for k in 1:nmax
        n = -k
        down_a_dd = dc_sub(
            _dd_factor(p, offset, n, 1 + p.s), ieps)
        down_d_dd = dc_add(
            _dd_factor(p, offset, n, 1 - p.s), ieps)
        down_e_dd = _dd_twice(p, offset, 0, n + 2)
        rational = dc_neg(dc_div(
            dc_mul(down_d_dd, down_e_dd),
            dc_mul(down_a_dd, DDComplex(-n)),
        ))
        down_term = dc_mul(
            down_term,
            dc_mul(negative[k], rational),
        )
        down_total = dc_add(down_total, down_term)
        down_absolute += abs(dc_value(down_term))
        if k == check_nmax
            down_check = down_total
            down_check_absolute = down_absolute
        end
    end
    up_check === nothing && throw(MSTCertificateError(
        :dd_amplitude_sum, "missing DD MST upper check sum."))
    down_check === nothing && throw(MSTCertificateError(
        :dd_amplitude_sum, "missing DD MST lower check sum."))
    ksum = (
        (
            _dd_sum_result(up_total, up_absolute),
            _dd_sum_result(down_total, down_absolute),
        ),
        (
            _dd_sum_result(up_check, up_check_absolute),
            _dd_sum_result(down_check, down_check_absolute),
        ),
    )

    itau = dc_imul(DDComplex(p.tau))
    d_positive(n) = dc_div(
        dc_mul(
            dc_add(_dd_factor(p, offset, n, 1 + p.s), ieps),
            dc_add(_dd_factor(p, offset, n, 1), itau),
        ),
        dc_mul(
            dc_sub(_dd_factor(p, offset, n, 1 - p.s), ieps),
            dc_sub(_dd_factor(p, offset, n, 1), itau),
        ),
    )
    d_negative(n) = dc_div(
        dc_mul(
            dc_sub(_dd_factor(p, offset, n, 1 - p.s), ieps),
            dc_sub(_dd_factor(p, offset, n, 1), itau),
        ),
        dc_mul(
            dc_add(_dd_factor(p, offset, n, 1 + p.s), ieps),
            dc_add(_dd_factor(p, offset, n, 1), itau),
        ),
    )
    d_base = exp(
        _amp_lgamma(1 + nu + p.s + I * epsc) +
        _amp_lgamma(1 + nu + I * tau) -
        _amp_lgamma(1 + nu - p.s - I * epsc) -
        _amp_lgamma(1 + nu - I * tau)
    )
    dsum = _dd_sum_pair(
        d_base,
        positive,
        negative,
        d_positive,
        d_negative,
        nmax,
        check_nmax,
    )
    return fsum, asum, ksum, dsum
end

# Double-double version of _amp_pair at nu = l + offset.
function _dd_amp_pair(
    data::MSTSeriesData,
    offset::DDComplex;
    nmax::Int=MST_AMP_NMAX,
    check_nmax::Int=MST_AMP_CHECK_NMAX,
)
    p = data.params
    offset2 = dc_sub(
        DDComplex(Float64(-2 * p.l - 1)),
        offset,
    )
    nu2 = ComplexF64(p.l) + dc_value(offset2)
    data2 = MSTSeriesData(MSTParams(
        p.s, p.l, p.m, p.a, p.omega, p.lambda, nu2))
    fsum, asum, ksum1, dsum1 = _dd_amp_sums(
        data, offset; nmax, check_nmax)
    _, _, ksum2, dsum2 = _dd_amp_sums(
        data2, offset2; nmax, check_nmax)
    full = _amp_build(
        data, data2, fsum[1], asum[1],
        ksum1[1], ksum2[1], dsum1[1], dsum2[1])
    check = _amp_build(
        data, data2, fsum[2], asum[2],
        ksum1[2], ksum2[2], dsum1[2], dsum2[2])
    return full, check
end

function _amp_logsum_values(logs)
    all(z -> isfinite(real(z)) && isfinite(imag(z)), logs) ||
        throw(MSTCertificateError(
            :nia_amplitude_sum,
            "nonfinite analytic MST amplitude summand.",
        ))
    scale = maximum(real, logs)
    total = ComplexF64(0)
    correction = ComplexF64(0)
    absolute_sum = 0.0
    # Each term exp(value - scale) carries the rounding of its exponent,
    # eps (|value| + |scale|), and of exp; compensated summation adds about
    # 2 eps sum|t|.
    weighted = 0.0
    for value in logs
        term = ComplexF64(exp(value - scale))
        y = term - correction
        updated = total + y
        correction = (updated - total) - y
        total = updated
        absolute_sum += abs(term)
        weighted += abs(term) * (1 + abs(value) + abs(scale))
    end
    norm = abs(total)
    isfinite(norm) && !iszero(norm) ||
        throw(MSTCertificateError(
            :nia_amplitude_sum,
            "zero or nonfinite analytic MST amplitude sum.",
        ))
    return MSTAmplitudeSum(
        ComplexF64(scale + log(norm) + I * angle(total)),
        Float64(absolute_sum / norm),
        eps(Float64) * (1 + (weighted + 2 * absolute_sum) / norm),
    )
end

function _amp_logsum(logtermfun, nmin::Int, nmax::Int)
    return _amp_logsum_values(
        ComplexF64[logtermfun(n) for n in nmin:nmax])
end

# Log-sum of a_n times a weight accumulated from exp(ratio_log(n)) outward from n = 0;
# _amp_rec_pair also returns the sum over [check_min, check_max].
function _amp_rec_sum(
    ratio_log,
    data::MSTSeriesData,
    nmin::Int,
    nmax::Int,
    base_extra,
)
    nmin <= 0 <= nmax || throw(ArgumentError(
        "MST amplitude recurrence range must contain zero."))
    logs = Vector{ComplexF64}(undef, nmax - nmin + 1)
    at(n) = n - nmin + 1
    extra = ComplexF64(base_extra)
    logs[at(0)] = log_series_coefficient!(data, 0) + extra
    for n in 0:(nmax - 1)
        extra += ComplexF64(ratio_log(n))
        logs[at(n + 1)] =
            log_series_coefficient!(data, n + 1) + extra
    end
    extra = ComplexF64(base_extra)
    for n in -1:-1:nmin
        extra -= ComplexF64(ratio_log(n))
        logs[at(n)] = log_series_coefficient!(data, n) + extra
    end
    return _amp_logsum_values(logs)
end

# Log-sums of `logs` over the full range and over [check_min, check_max].
function _amp_pair_sum(
    logs,
    nmin::Int,
    check_min::Int,
    check_max::Int,
)
    first_index = check_min - nmin + 1
    last_index = check_max - nmin + 1
    1 <= first_index <= last_index <= length(logs) ||
        throw(ArgumentError("invalid MST amplitude check range."))
    return (
        _amp_logsum_values(logs),
        _amp_logsum_values(@view logs[first_index:last_index]),
    )
end

function _amp_rec_pair(
    ratio_log,
    data::MSTSeriesData,
    nmin::Int,
    nmax::Int,
    check_min::Int,
    check_max::Int,
    base_extra,
)
    nmin <= check_min <= 0 <= check_max <= nmax ||
        throw(ArgumentError(
            "MST amplitude recurrence ranges must contain zero."))
    logs = Vector{ComplexF64}(undef, nmax - nmin + 1)
    at(n) = n - nmin + 1
    extra = ComplexF64(base_extra)
    logs[at(0)] = log_series_coefficient!(data, 0) + extra
    for n in 0:(nmax - 1)
        extra += ComplexF64(ratio_log(n))
        logs[at(n + 1)] =
            log_series_coefficient!(data, n + 1) + extra
    end
    extra = ComplexF64(base_extra)
    for n in -1:-1:nmin
        extra -= ComplexF64(ratio_log(n))
        logs[at(n)] = log_series_coefficient!(data, n) + extra
    end
    return _amp_pair_sum(logs, nmin, check_min, check_max)
end

function _amp_logadd(first::ComplexF64, second::ComplexF64)
    scale = max(real(first), real(second))
    value = exp(first - scale) + exp(second - scale)
    !iszero(value) ||
        throw(MSTCertificateError(
            :nia_amplitude_cancellation,
            "exact cancellation in analytic MST amplitude formula.",
        ))
    return ComplexF64(scale + log(abs(value)) + I * angle(value))
end

_amp_logminus(first::ComplexF64, second::ComplexF64) =
    _amp_logadd(first, second + I * pi)

function _amp_lgamma(z)
    value = ComplexF64(loggamma(ComplexF64(z)))
    isfinite(real(value)) && isfinite(imag(value)) ||
        throw(MSTCertificateError(
            :nia_amplitude_gamma,
            "nonfinite analytic MST loggamma factor.",
        ))
    return value
end

_amp_lfactorial(n::Int) = _amp_lgamma(n + 1)
_amp_lsin(z) = log(ComplexF64(sin(ComplexF64(z))))

# Sum of a_n over |n| <= nmax.
function _amp_coeff_sum(data::MSTSeriesData, nmax::Int)
    return _amp_logsum(-nmax, nmax) do n
        log_series_coefficient!(data, n)
    end
end

function _amp_coeff_pair(
    data::MSTSeriesData,
    nmax::Int,
    check_nmax::Int,
)
    logs = ComplexF64[
        log_series_coefficient!(data, n) for n in -nmax:nmax]
    return _amp_pair_sum(logs, -nmax, -check_nmax, check_nmax)
end

# A_- sum of a_n (-1)^n (nu+1+s-i epsilon)_n / (nu+1-s+i epsilon)_n; the _fast and _pair
# versions build the weights by recurrence.
function _amp_aminus_sum(data::MSTSeriesData, nmax::Int)
    p = data.params
    nu = _nu_shift(p)
    epsc = ComplexF64(p.epsilon)
    return _amp_logsum(-nmax, nmax) do n
        log_series_coefficient!(data, n) +
            I * pi * n +
            log_pochhammer(nu + 1 + p.s - I * epsc, n) -
            log_pochhammer(nu + 1 - p.s + I * epsc, n)
    end
end

function _amp_aminus_sum_fast(data::MSTSeriesData, nmax::Int)
    p = data.params
    nu = _nu_shift(p)
    epsc = ComplexF64(p.epsilon)
    numerator = nu + 1 + p.s - I * epsc
    denominator = nu + 1 - p.s + I * epsc
    return _amp_rec_sum(data, -nmax, nmax, 0.0 + 0.0im) do n
        I * pi +
            log(ComplexF64(numerator + n)) -
            log(ComplexF64(denominator + n))
    end
end

function _amp_aminus_pair(
    data::MSTSeriesData,
    nmax::Int,
    check_nmax::Int,
)
    p = data.params
    nu = _nu_shift(p)
    epsc = ComplexF64(p.epsilon)
    numerator = nu + 1 + p.s - I * epsc
    denominator = nu + 1 - p.s + I * epsc
    return _amp_rec_pair(
        data, -nmax, nmax, -check_nmax, check_nmax, 0.0 + 0.0im,
    ) do n
        I * pi +
            log(ComplexF64(numerator + n)) -
            log(ComplexF64(denominator + n))
    end
end

# Upper (n >= 0) and lower (n <= 0) sums of K_nu.
function _amp_knu_sums(data::MSTSeriesData, nmax::Int)
    p = data.params
    nu = _nu_shift(p)
    epsc = ComplexF64(p.epsilon)
    tau = ComplexF64(p.tau)
    up = _amp_logsum(0, nmax) do n
        I * pi * n +
            _amp_lgamma(1 + n + p.s + I * epsc + nu) +
            _amp_lgamma(1 + n + 2nu) +
            _amp_lgamma(1 + n + nu + I * tau) -
            _amp_lfactorial(n) -
            _amp_lgamma(1 + n - p.s - I * epsc + nu) -
            _amp_lgamma(1 + n + nu - I * tau) +
            log_series_coefficient!(data, n)
    end
    down = _amp_logsum(-nmax, 0) do n
        I * pi * n +
            log_pochhammer(1 + p.s - I * epsc + nu, n) -
            _amp_lfactorial(-n) -
            log_pochhammer(1 - p.s + I * epsc + nu, n) -
            log_pochhammer(2 + 2nu, n) +
            log_series_coefficient!(data, n)
    end
    return up, down
end

function _amp_knu_sums_fast(data::MSTSeriesData, nmax::Int)
    p = data.params
    nu = _nu_shift(p)
    epsc = ComplexF64(p.epsilon)
    tau = ComplexF64(p.tau)

    up_a = 1 + p.s + I * epsc + nu
    up_b = 1 + 2nu
    up_c = 1 + nu + I * tau
    up_d = 1 - p.s - I * epsc + nu
    up_e = 1 + nu - I * tau
    up_base =
        _amp_lgamma(up_a) +
        _amp_lgamma(up_b) +
        _amp_lgamma(up_c) -
        _amp_lgamma(up_d) -
        _amp_lgamma(up_e)
    up = _amp_rec_sum(data, 0, nmax, up_base) do n
        I * pi +
            log(ComplexF64(up_a + n)) +
            log(ComplexF64(up_b + n)) +
            log(ComplexF64(up_c + n)) -
            log(ComplexF64(n + 1)) -
            log(ComplexF64(up_d + n)) -
            log(ComplexF64(up_e + n))
    end

    down_a = 1 + p.s - I * epsc + nu
    down_d = 1 - p.s + I * epsc + nu
    down_e = 2 + 2nu
    down = _amp_rec_sum(data, -nmax, 0, 0.0 + 0.0im) do n
        I * pi +
            log(ComplexF64(down_a + n)) +
            log(ComplexF64(-n)) -
            log(ComplexF64(down_d + n)) -
            log(ComplexF64(down_e + n))
    end
    return up, down
end

function _amp_knu_pair(
    data::MSTSeriesData,
    nmax::Int,
    check_nmax::Int,
)
    p = data.params
    nu = _nu_shift(p)
    epsc = ComplexF64(p.epsilon)
    tau = ComplexF64(p.tau)

    up_a = 1 + p.s + I * epsc + nu
    up_b = 1 + 2nu
    up_c = 1 + nu + I * tau
    up_d = 1 - p.s - I * epsc + nu
    up_e = 1 + nu - I * tau
    up_base =
        _amp_lgamma(up_a) +
        _amp_lgamma(up_b) +
        _amp_lgamma(up_c) -
        _amp_lgamma(up_d) -
        _amp_lgamma(up_e)
    up = _amp_rec_pair(
        data, 0, nmax, 0, check_nmax, up_base,
    ) do n
        I * pi +
            log(ComplexF64(up_a + n)) +
            log(ComplexF64(up_b + n)) +
            log(ComplexF64(up_c + n)) -
            log(ComplexF64(n + 1)) -
            log(ComplexF64(up_d + n)) -
            log(ComplexF64(up_e + n))
    end

    down_a = 1 + p.s - I * epsc + nu
    down_d = 1 - p.s + I * epsc + nu
    down_e = 2 + 2nu
    down = _amp_rec_pair(
        data, -nmax, 0, -check_nmax, 0, 0.0 + 0.0im,
    ) do n
        I * pi +
            log(ComplexF64(down_a + n)) +
            log(ComplexF64(-n)) -
            log(ComplexF64(down_d + n)) -
            log(ComplexF64(down_e + n))
    end
    return (
        (up[1], down[1]),
        (up[2], down[2]),
    )
end

# Sum of a_n Gamma(1+nu+n+s+i epsilon) Gamma(1+nu+n+i tau) /
# (Gamma(1+nu+n-s-i epsilon) Gamma(1+nu+n-i tau)), used in up_inc.
function _amp_dsum(data::MSTSeriesData, nmax::Int)
    p = data.params
    nu = _nu_shift(p)
    epsc = ComplexF64(p.epsilon)
    tau = ComplexF64(p.tau)
    return _amp_logsum(-nmax, nmax) do n
        _amp_lgamma(1 + nu + n + p.s + I * epsc) +
            _amp_lgamma(1 + nu + n + I * tau) -
            _amp_lgamma(1 + nu + n - p.s - I * epsc) -
            _amp_lgamma(1 + nu + n - I * tau) +
            log_series_coefficient!(data, n)
    end
end

function _amp_dsum_fast(data::MSTSeriesData, nmax::Int)
    p = data.params
    nu = _nu_shift(p)
    epsc = ComplexF64(p.epsilon)
    tau = ComplexF64(p.tau)
    numerator_a = 1 + nu + p.s + I * epsc
    numerator_b = 1 + nu + I * tau
    denominator_a = 1 + nu - p.s - I * epsc
    denominator_b = 1 + nu - I * tau
    base =
        _amp_lgamma(numerator_a) +
        _amp_lgamma(numerator_b) -
        _amp_lgamma(denominator_a) -
        _amp_lgamma(denominator_b)
    return _amp_rec_sum(data, -nmax, nmax, base) do n
        log(ComplexF64(numerator_a + n)) +
            log(ComplexF64(numerator_b + n)) -
            log(ComplexF64(denominator_a + n)) -
            log(ComplexF64(denominator_b + n))
    end
end

function _amp_d_pair(
    data::MSTSeriesData,
    nmax::Int,
    check_nmax::Int,
)
    p = data.params
    nu = _nu_shift(p)
    epsc = ComplexF64(p.epsilon)
    tau = ComplexF64(p.tau)
    numerator_a = 1 + nu + p.s + I * epsc
    numerator_b = 1 + nu + I * tau
    denominator_a = 1 + nu - p.s - I * epsc
    denominator_b = 1 + nu - I * tau
    base =
        _amp_lgamma(numerator_a) +
        _amp_lgamma(numerator_b) -
        _amp_lgamma(denominator_a) -
        _amp_lgamma(denominator_b)
    return _amp_rec_pair(
        data, -nmax, nmax, -check_nmax, check_nmax, base,
    ) do n
        log(ComplexF64(numerator_a + n)) +
            log(ComplexF64(numerator_b + n)) -
            log(ComplexF64(denominator_a + n)) -
            log(ComplexF64(denominator_b + n))
    end
end

# log K_nu from its upper and lower sums, with its absolute error.
function _amp_klog(data::MSTSeriesData, sums)
    p = data.params
    nu = _nu_shift(p)
    epsc = ComplexF64(p.epsilon)
    tau = ComplexF64(p.tau)
    epsp = (tau + epsc) / 2
    up, down = sums
    pieces = (
        -nu * log(2.0),
        I * epsc * p.kappa,
        (p.s - nu) * log(epsc * p.kappa),
        _amp_lgamma(1 - p.s - 2I * epsp),
        _amp_lgamma(2 + 2nu),
        -_amp_lgamma(1 - p.s + I * epsc + nu),
        -_amp_lgamma(1 + p.s + I * epsc + nu),
        -_amp_lgamma(1 + nu + I * tau),
        up.logvalue,
        -down.logvalue,
    )
    # Absolute error of the log: the two sums' relative errors and the
    # rounding of each piece and of the additions.
    error = up.error + down.error + eps(Float64) * 2 * sum(abs, pieces)
    return ComplexF64(sum(pieces)), Float64(error)
end

# log(exp(first) + exp(second)) and its absolute error from those of the
# two parts: a cancelling combination amplifies them by
# (|e^first| + |e^second|) / |e^first + e^second|.
function _amp_logadd_error(first, first_error, second, second_error)
    scale = max(real(first), real(second))
    a, b = exp(first - scale), exp(second - scale)
    return (abs(a) * first_error + abs(b) * second_error) / abs(a + b) +
        eps(Float64) * (abs(a) + abs(b)) / abs(a + b)
end

# The MST closed forms for the amplitudes hold for Re epsilon > 0: their
# Coulomb asymptotics are matched on that side of the anti-Stokes lines, and
# for Re epsilon < 0 they miss a Stokes contribution (the reflection
# amplitude of IN), not only a branch of log epsilon. There the amplitudes
# follow from the exact reflection of the Teukolsky equation on real r,
# R_{s,l,m}(omega) = conj R_{s,l,-m}(-conj(omega)) (lambda and nu conjugate,
# and so do the series' sums): the closed forms are evaluated at the
# reflected parameters and conjugated. The split lies on the negative
# imaginary omega axis, the amplitudes' own branch cut, and the sign of the
# zero of Re omega picks the side (the reflection maps -0.0 to +0.0).
_amp_reflected(p::MSTParams) = signbit(real(p.omega))
_amp_reflect(p::MSTParams) = MSTParams(p.s, p.l, -p.m, p.a, -conj(p.omega),
    conj(p.lambda), conj(p.nu), conj(p.nu_offset))
_amp_conj(x::MSTAmplitudeSum) =
    MSTAmplitudeSum(conj(x.logvalue), x.condition, x.error)
_amp_conj(x::Tuple) = map(_amp_conj, x)
_amp_conj(x::NamedTuple) = map(v -> v isa Complex ? conj(v) : v, x)

# Logs of the six asymptotic amplitudes (in/up incidence, transmission, reflection)
# and their absolute errors from the coefficient sum, the A_- sum, the K sums at nu and
# -nu - 1 and the D sums (Sasaki & Tagoshi, Living Rev. Relativ. 6, 6 (2003)).
function _amp_build(
    data::MSTSeriesData,
    data2::MSTSeriesData,
    fsum,
    asum,
    ksum1,
    ksum2,
    dsum1,
    dsum2,
)
    _amp_reflected(data.params) && return _amp_conj(_amp_build(
        MSTSeriesData(_amp_reflect(data.params)),
        MSTSeriesData(_amp_reflect(data2.params)),
        map(_amp_conj, (fsum, asum, ksum1, ksum2, dsum1, dsum2))...))
    p = data.params
    nu = _nu_shift(p)
    nu2 = -nu - 1
    epsc = ComplexF64(p.epsilon)
    tau = ComplexF64(p.tau)
    kappa = ComplexF64(p.kappa)
    epsp = (tau + epsc) / 2
    s = p.s

    k1, k1_error = _amp_klog(data, ksum1)
    k2, k2_error = _amp_klog(data2, ksum2)

    in_trans =
        s * log(4.0) + 2s * log(kappa) +
        I * (epsc + tau) * kappa *
            (0.5 + log(kappa) / (1 + kappa)) +
        fsum.logvalue
    aminus =
        (-s - 1 + I * epsc) * log(2.0) -
        pi * epsc / 2 -
        I * pi * (nu + 1 + s) / 2 +
        asum.logvalue
    up_trans =
        (-1 - 2s) * log(epsc / 2) +
        I * epsc * (log(epsc) - (1 - kappa) / 2) +
        aminus
    in_ref =
        up_trans +
        _amp_logadd(k1, I * pi / 2 + I * pi * nu + k2)

    common_d2 =
        I * kappa * (epsc + tau) *
            (1 + kappa + 2log(kappa)) / (2 * (1 + kappa)) +
        2s * log(2kappa) -
        _amp_lsin(pi * I * (epsc + tau)) +
        fsum.logvalue
    d2 =
        I * pi + common_d2 +
        _amp_lsin(pi * (nu - I * epsc)) +
        _amp_lsin(pi * (nu - I * tau)) -
        _amp_lsin(2pi * nu)
    d22 =
        I * pi + common_d2 +
        _amp_lsin(pi * (nu2 - I * epsc)) +
        _amp_lsin(pi * (nu2 - I * tau)) -
        _amp_lsin(2pi * nu2)
    up_ref =
        -pi * epsc - I * pi * s - _amp_lsin(2pi * nu) +
        _amp_logadd(
            -I * pi * nu +
                _amp_lsin(pi * (nu - s + I * epsc)) -
                k1 + d2,
            -I * pi / 2 +
                _amp_lsin(pi * (nu + s - I * epsc)) -
                k2 + d22,
        )

    aplus =
        (-1 + s - I * epsc) * log(2.0) -
        pi * epsc / 2 +
        I * pi * (1 - s + nu) / 2 +
        _amp_lgamma(1 - s + I * epsc + nu) -
        _amp_lgamma(1 + s - I * epsc + nu) +
        fsum.logvalue
    in_inc =
        -log(epsc / 2) +
        _amp_logminus(
            k1,
            I * pi / 2 - I * pi * nu +
                _amp_lsin(pi * (nu - s + I * epsc)) -
                _amp_lsin(pi * (nu + s - I * epsc)) +
                k2,
        ) -
        I * epsc * (log(epsc) - (1 - kappa) / 2) +
        aplus

    common_d1 =
        -I * kappa * (epsc + tau) *
            (1 + kappa + 2log(kappa)) / (2 * (1 + kappa)) -
        _amp_lsin(pi * I * (epsc + tau)) +
        _amp_lgamma(1 - s - I * (epsc + tau)) -
        _amp_lgamma(1 + s + I * epsc + I * tau)
    d1 =
        common_d1 +
        _amp_lsin(pi * (nu + I * epsc)) +
        _amp_lsin(pi * (nu + I * tau)) -
        _amp_lsin(2pi * nu) +
        dsum1.logvalue
    d12 =
        common_d1 +
        _amp_lsin(pi * (nu2 + I * epsc)) +
        _amp_lsin(pi * (nu2 + I * tau)) -
        _amp_lsin(2pi * nu2) +
        dsum2.logvalue
    up_inc =
        -pi * epsc - I * pi * s - _amp_lsin(2pi * nu) +
        _amp_logadd(
            -I * pi * nu +
                _amp_lsin(pi * (nu - s + I * epsc)) -
                k1 + d1,
            -I * pi / 2 +
                _amp_lsin(pi * (nu + s - I * epsc)) -
                k2 + d12,
        )

    # Absolute errors of the six logs (relative errors of the amplitudes):
    # the sums', the K factors' and the cancellation of each two-term
    # combination; the remaining factors are single elementary functions.
    rounding(x) = eps(Float64) * 2 * (1 + abs(x))
    in_trans_error = fsum.error + rounding(in_trans)
    up_trans_error = asum.error + rounding(up_trans)
    in_ref_error = up_trans_error + _amp_logadd_error(
        k1, k1_error, I * pi / 2 + I * pi * nu + k2, k2_error) +
        rounding(in_ref)
    in_inc_error = fsum.error + _amp_logadd_error(
        k1, k1_error,
        I * pi / 2 - I * pi * nu +
            _amp_lsin(pi * (nu - s + I * epsc)) -
            _amp_lsin(pi * (nu + s - I * epsc)) + k2 + I * pi,
        k2_error) + rounding(in_inc)
    up_inc_error = _amp_logadd_error(
        -I * pi * nu + _amp_lsin(pi * (nu - s + I * epsc)) - k1 + d1,
        k1_error + dsum1.error,
        -I * pi / 2 + _amp_lsin(pi * (nu + s - I * epsc)) - k2 + d12,
        k2_error + dsum2.error) + rounding(up_inc)
    up_ref_error = _amp_logadd_error(
        -I * pi * nu + _amp_lsin(pi * (nu - s + I * epsc)) - k1 + d2,
        k1_error + fsum.error,
        -I * pi / 2 + _amp_lsin(pi * (nu + s - I * epsc)) - k2 + d22,
        k2_error + fsum.error) + rounding(up_ref)
    return (;
        in_inc,
        in_trans,
        in_ref,
        up_inc,
        up_trans,
        up_ref,
        in_inc_error,
        in_trans_error,
        in_ref_error,
        up_inc_error,
        up_trans_error,
        up_ref_error,
        max_condition=maximum((
            fsum.condition,
            asum.condition,
            ksum1[1].condition,
            ksum1[2].condition,
            ksum2[1].condition,
            ksum2[2].condition,
            dsum1.condition,
            dsum2.condition,
        )),
    )
end

# The series at -nu - 1 of a representation, computed on its own.
function _amp_mirror_data(data::MSTSeriesData)
    p = data.params
    return MSTSeriesData(MSTParams(
        p.s, p.l, p.m, p.a, p.omega, p.lambda, -_nu_shift(p) - 1))
end

# Residual of the three-term recurrence at a representation's anchor
# (a_0 = 1), relative to its terms: a representation whose a_n are
# recomputed from an anchor in the tail of |a_k| carries that error into
# every coefficient.
function _amp_anchor_residual(data::MSTSeriesData)
    p = data.params
    up = alpha_mst(0, p) * exp(log_series_coefficient!(data, 1))
    centre = beta_mst(0, p)
    down = gamma_mst(0, p) * exp(log_series_coefficient!(data, -1))
    return Float64(abs(up + centre + down) /
        (abs(up) + abs(centre) + abs(down)))
end

# The coefficients of the representation shifted by d from source (nu' =
# nu + d), taken from source: a'_k = a_{k+d} / a_d, and for its -nu' - 1
# series a''_k = a'_{-k} (a_n^{-nu-1} = a_{-n}^{nu}). Indices beyond
# +-nmax continue by the representation's own recurrence.
function _amp_shared_data(source::MSTSeriesData, params::MSTParams, d::Int,
        nmax::Int)
    data = MSTSeriesData(params)
    mirror = _amp_mirror_data(data)
    anchor = log_series_coefficient!(source, d)
    anchor_error = source.coeff_errors[d]
    for k in -nmax:nmax
        k == 0 && continue
        value = log_series_coefficient!(source, k + d) - anchor
        error = source.coeff_errors[k + d] + anchor_error
        relerr = get(source.cf_relerrs, k + d, eps(Float64))
        ok = get(source.cf_ok, k + d, true)
        for (target, index) in ((data, k), (mirror, -k))
            target.log_coeffs[index] = value
            target.coeff_errors[index] = error
            target.cf_relerrs[index] = relerr
            target.cf_ok[index] = ok
        end
    end
    return data, mirror
end

# Amplitude logs from the sums over |n| <= nmax; `fast` builds the summand weights by
# ratio recurrences.
@noinline function _amp_logs(
    data::MSTSeriesData,
    nmax::Int;
    fast::Bool=true,
    data2::MSTSeriesData=_amp_mirror_data(data),
)
    fsum = _amp_coeff_sum(data, nmax)
    asum = fast ?
        _amp_aminus_sum_fast(data, nmax) :
        _amp_aminus_sum(data, nmax)
    ksum1 = fast ?
        _amp_knu_sums_fast(data, nmax) :
        _amp_knu_sums(data, nmax)
    ksum2 = fast ?
        _amp_knu_sums_fast(data2, nmax) :
        _amp_knu_sums(data2, nmax)
    dsum1 = fast ?
        _amp_dsum_fast(data, nmax) :
        _amp_dsum(data, nmax)
    dsum2 = fast ?
        _amp_dsum_fast(data2, nmax) :
        _amp_dsum(data2, nmax)
    return _amp_build(
        data, data2, fsum, asum, ksum1, ksum2, dsum1, dsum2)
end

# Amplitude logs at order nmax and at the check order check_nmax.
@noinline function _amp_pair(
    data::MSTSeriesData,
    nmax::Int,
    check_nmax::Int,
)
    0 <= check_nmax <= nmax ||
        throw(ArgumentError("invalid MST amplitude check order."))
    p = data.params
    nu2 = -_nu_shift(p) - 1
    data2 = MSTSeriesData(MSTParams(
        p.s, p.l, p.m, p.a, p.omega, p.lambda, nu2))
    fsum = _amp_coeff_pair(data, nmax, check_nmax)
    asum = _amp_aminus_pair(data, nmax, check_nmax)
    ksum1 = _amp_knu_pair(data, nmax, check_nmax)
    ksum2 = _amp_knu_pair(data2, nmax, check_nmax)
    dsum1 = _amp_d_pair(data, nmax, check_nmax)
    dsum2 = _amp_d_pair(data2, nmax, check_nmax)
    full = _amp_build(
        data,
        data2,
        fsum[1],
        asum[1],
        ksum1[1],
        ksum2[1],
        dsum1[1],
        dsum2[1],
    )
    check = _amp_build(
        data,
        data2,
        fsum[2],
        asum[2],
        ksum1[2],
        ksum2[2],
        dsum1[2],
        dsum2[2],
    )
    return full, check
end

# Teukolsky incidence/transmission and reflection/transmission ratios of :IN or :UP.
function _amp_teuk_pair(logs, branch::Symbol)
    if branch == :IN
        return (
            ComplexF64(exp(logs.in_inc - logs.in_trans)),
            ComplexF64(exp(logs.in_ref - logs.in_trans)),
        )
    elseif branch == :UP
        return (
            ComplexF64(exp(logs.up_inc - logs.up_trans)),
            ComplexF64(exp(logs.up_ref - logs.up_trans)),
        )
    end
    throw(ArgumentError("analytic MST amplitude branch must be :IN or :UP."))
end

# Relative errors of the two amplitude ratios of _amp_teuk_pair (the
# conversion to GSN normalisation adds rounding only).
_amp_ratio_errors(logs, branch::Symbol) = branch == :IN ?
    (logs.in_inc_error + logs.in_trans_error,
        logs.in_ref_error + logs.in_trans_error) :
    (logs.up_inc_error + logs.up_trans_error,
        logs.up_ref_error + logs.up_trans_error)

# Converts Teukolsky amplitude ratios to GSN ones with the ConversionFactors amplitudes.
function _amp_gsn_pair(p::MSTParams, branch::Symbol, teuk)
    cf = _conversion_module()
    trans, inc, ref = if branch == :IN
        (
            getfield(cf, :Btrans)(p.s, p.m, p.a, p.omega, p.lambda),
            getfield(cf, :Binc)(p.s, p.m, p.a, p.omega, p.lambda),
            getfield(cf, :Bref)(p.s, p.m, p.a, p.omega, p.lambda),
        )
    else
        (
            getfield(cf, :Ctrans)(p.s, p.m, p.a, p.omega, p.lambda),
            getfield(cf, :Cinc)(p.s, p.m, p.a, p.omega, p.lambda),
            getfield(cf, :Cref)(p.s, p.m, p.a, p.omega, p.lambda),
        )
    end
    return (
        ComplexF64(teuk[1] * trans / inc),
        ComplexF64(teuk[2] * trans / ref),
    )
end

# Max-norm distance between two amplitude pairs, relative to the largest amplitude.
function _amp_distance(first, second)
    scale = max(
        abs(first[1]), abs(first[2]),
        abs(second[1]), abs(second[2]),
        floatmin(Float64),
    )
    return max(
        abs(first[1] - second[1]),
        abs(first[2] - second[2]),
    ) / scale
end

# Relative distance of each amplitude of a pair from the corresponding one of
# another evaluation (the max-norm distance is in units of the larger
# amplitude and says little about a suppressed one).
_amp_component_distance(value, other) =
    Float64.(abs.(value .- other) ./ abs.(value))

# Candidate with the smallest summed distance to the others, with its spread and its
# distance to the nearest other representation.
function _amp_medoid(candidates)
    best = nothing
    best_score = Inf
    best_spread = Inf
    for candidate in candidates
        distances = [
            _amp_distance(candidate.gsn, peer.gsn)
            for peer in candidates
        ]
        score = sum(distances)
        spread = maximum(distances)
        if score < best_score ||
                (score == best_score && spread < best_spread)
            best = candidate
            best_score = score
            best_spread = spread
        end
    end
    peers = [peer for peer in candidates if peer.shift != best.shift]
    nearest_peer = peers[argmin([_amp_distance(best.gsn, peer.gsn)
        for peer in peers])]
    nearest = _amp_distance(best.gsn, nearest_peer.gsn)
    return merge(best, (
        medoid_score=Float64(best_score),
        representation_spread=Float64(best_spread),
        nearest_agreement=Float64(nearest),
        nearest_components=_amp_component_distance(
            best.gsn, nearest_peer.gsn),
    ))
end

# Indices and distance of the closest pair of candidates.
function _amp_pair(candidates)
    pair = nothing
    distance = Inf
    for first_index in 1:(length(candidates) - 1)
        for second_index in (first_index + 1):length(candidates)
            candidate_distance = _amp_distance(
                candidates[first_index].gsn,
                candidates[second_index].gsn,
            )
            if candidate_distance < distance
                pair = (first_index, second_index)
                distance = candidate_distance
            end
        end
    end
    pair === nothing && throw(MSTCertificateError(
        :nia_amplitude_certificate,
        "analytic MST amplitude coherent pair is unavailable.",
    ))
    return pair, Float64(distance)
end

# Distance of a candidate from the same representation summed to MST_AMP_CHECK_NMAX.
function _amp_check(base, shared, candidate, branch)
    data, mirror = shared[candidate.shift]
    logs = _amp_logs(data, MST_AMP_CHECK_NMAX; data2=mirror)
    teuk = _amp_teuk_pair(logs, branch)
    gsn = _amp_gsn_pair(base, branch, teuk)
    return (distance=Float64(_amp_distance(candidate.gsn, gsn)),
        components=_amp_component_distance(candidate.gsn, gsn))
end

# nu continued by Newton on its characteristic equation at the anchor of
# p while the step and the residual decrease: the solver stops at a fixed
# residual, and the coefficients taken from this anchor carry nu's
# remaining error. The Newton step left at the result is that error.
function _amp_refine_nu(p::MSTParams)
    equation = nu -> _mst_eq_value(
        nu - p.l, p.s, p.l, p.m, p.a, p.omega, p.lambda)
    h = exp2(-18) * max(1.0, abs(p.nu))
    current = p.nu
    value = equation(current)
    previous = Inf
    while true
        derivative = (equation(current + h) - equation(current - h)) / (2h)
        step = value / derivative
        abs(step) < previous || return (nu=current, step)
        trial = current - step
        trial_value = equation(trial)
        abs(trial_value) < abs(value) || return (nu=current, step)
        previous = abs(step)
        current = trial
        value = trial_value
    end
end

# Relative change of a representation's amplitudes when nu (and with it
# the shared coefficients) moves by its remaining Newton step.
function _amp_nu_error(base, centered, source_shift, candidate, branch, step)
    moved = MSTParams(base.s, base.l, base.m, base.a, base.omega,
        base.lambda, base.nu - step)
    data, mirror = _amp_shared_data(
        _amp_data(moved, centered - step, source_shift),
        _amp_data(moved, centered - step, candidate.shift).params,
        candidate.shift - source_shift, MST_AMP_NMAX)
    logs = _amp_logs(data, MST_AMP_NMAX; data2=mirror)
    return _amp_component_distance(candidate.gsn,
        _amp_gsn_pair(moved, branch, _amp_teuk_pair(logs, branch)))
end

# MSTSeriesData of the representation nu = centered + shift.
function _amp_data(base::MSTParams, centered, shift)
    return MSTSeriesData(MSTParams(
        base.s,
        base.l,
        base.m,
        base.a,
        base.omega,
        base.lambda,
        centered + shift,
    ))
end

"""
    mst_nia_amplitudes(s, l, m, a, omega, lambda, branch; nu=nothing)

Analytic MST asymptotic amplitudes of `branch` (`:IN` or `:UP`) certified across the
representations nu + k, k in `MST_AMP_SHIFTS`, which share one coefficient set. The
medoid of the candidates is returned when spread, nearest agreement, truncation and
condition pass; otherwise the closest coherent pair is tried, and an
`MSTCertificateError` is thrown if it fails too. The result holds the Teukolsky and GSN
amplitude ratios (`teuk`, `gsn`), the transmission log and the certificate diagnostics.
"""
@noinline function mst_nia_amplitudes(
    s::Integer,
    l::Integer,
    m::Integer,
    a::Real,
    omega::Number,
    lambda,
    branch::Symbol;
    nu=nothing,
)
    omegac = ComplexF64(omega)
    nu_value = nu === nothing ?
        mst_nu(s, l, m, a, omegac, lambda) :
        ComplexF64(nu)
    base = MSTParams(s, l, m, a, omegac, lambda, nu_value)
    centered = nu_value - round(Int, real(nu_value))
    # One coefficient set for all representations: that of the anchor
    # (a representation's or the solver's own nu) that satisfies its
    # recurrence best, with nu refined there, reindexed to each
    # representation's anchor. Recomputing a_n from an anchor in the tail
    # of |a_k| loses digits that the running bounds do not see.
    sources = (MST_AMP_SHIFTS..., round(Int, real(nu_value)))
    source_shift = sources[argmin(map(shift -> _amp_anchor_residual(
        _amp_data(base, centered, shift)), sources))]
    refinement = _amp_refine_nu(_amp_data(base, centered, source_shift).params)
    refined = refinement.nu
    nu_value += refined - (centered + source_shift)
    base = MSTParams(s, l, m, a, omegac, lambda, nu_value)
    centered = refined - source_shift
    source = _amp_data(base, centered, source_shift)
    source_residual = _amp_anchor_residual(source)
    shared = Dict(shift => _amp_shared_data(source,
            _amp_data(base, centered, shift).params,
            shift - source_shift, MST_AMP_NMAX)
        for shift in MST_AMP_SHIFTS)
    candidates = map(MST_AMP_SHIFTS) do shift
        data, mirror = shared[shift]
        logs = _amp_logs(data, MST_AMP_NMAX; data2=mirror)
        teuk = _amp_teuk_pair(logs, branch)
        gsn = _amp_gsn_pair(base, branch, teuk)
        all(_finite_complex, (teuk..., gsn...)) ||
            throw(MSTCertificateError(
                :nia_amplitude_nonfinite,
                "nonfinite analytic MST amplitude candidate.",
            ))
        return (;
            shift,
            nu=nu_value,
            teuk,
            gsn,
            transmission_log=branch == :IN ?
                ComplexF64(logs.in_trans) : ComplexF64(logs.up_trans),
            max_condition=Float64(logs.max_condition),
            gsn_error=_amp_ratio_errors(logs, branch),
            coefficient_residual=source_residual,
            coefficient_source_shift=source_shift,
        )
    end
    medoid = _amp_medoid(candidates)
    medoid_check = _amp_check(base, shared, medoid, branch)
    medoid_truncation = medoid_check.distance
    strict_accepted =
        medoid.representation_spread <= MST_AMP_SPREAD_MAX &&
        medoid.nearest_agreement <= MST_AMP_NEAREST_MAX &&
        medoid_truncation <= MST_AMP_TRUNCATION_MAX &&
        medoid.max_condition <= MST_AMP_CONDITION_MAX
    if strict_accepted
        return merge(medoid, (
            truncation_agreement=medoid_truncation,
            component_agreement=max.(medoid.nearest_components,
                medoid_check.components) .+ _amp_nu_error(base, centered,
                source_shift, medoid, branch, refinement.step),
            nmax=MST_AMP_NMAX,
            check_nmax=MST_AMP_CHECK_NMAX,
            certificate_kind=:global_medoid,
            certificate_accepted=true,
        ))
    end

    pair, pair_agreement = _amp_pair(candidates)
    first = candidates[pair[1]]
    second = candidates[pair[2]]
    first_check = _amp_check(base, shared, first, branch)
    second_check = _amp_check(base, shared, second, branch)
    selected, other, selected_check = (
        first_check.distance,
        first.max_condition,
    ) <= (
        second_check.distance,
        second.max_condition,
    ) ? (first, second, first_check) : (second, first, second_check)
    truncation_agreement = selected_check.distance
    accepted =
        pair_agreement <= MST_AMP_NEAREST_MAX &&
        truncation_agreement <= MST_AMP_TRUNCATION_MAX &&
        selected.max_condition <= MST_AMP_CONDITION_MAX
    accepted ||
        throw(MSTCertificateError(
            :nia_amplitude_certificate,
            "analytic MST amplitude certificate rejected: " *
            "spread=$(medoid.representation_spread), " *
            "nearest=$pair_agreement, " *
            "truncation=$truncation_agreement, " *
            "condition=$(selected.max_condition).",
        ))
    return merge(selected, (
        medoid_score=medoid.medoid_score,
        representation_spread=medoid.representation_spread,
        nearest_agreement=pair_agreement,
        truncation_agreement=Float64(truncation_agreement),
        component_agreement=max.(
            _amp_component_distance(selected.gsn, other.gsn),
            selected_check.components) .+ _amp_nu_error(base, centered,
            source_shift, selected, branch, refinement.step),
        nmax=MST_AMP_NMAX,
        check_nmax=MST_AMP_CHECK_NMAX,
        certificate_kind=:coherent_pair,
        certificate_accepted=true,
    ))
end

# Amplitudes at the principal nu with a truncation check at MST_AMP_CHECK_NMAX.
@noinline function _principal_amplitudes(
    data::MSTSeriesData,
    branch::Symbol,
)
    # nu continued to its noise floor at its own anchor, and its remaining
    # Newton step propagated into the error, as in mst_nia_amplitudes.
    p = data.params
    refinement = _amp_refine_nu(p)
    with_nu(nu) = MSTSeriesData(
        MSTParams(p.s, p.l, p.m, p.a, p.omega, p.lambda, nu))
    data = with_nu(refinement.nu)
    base = data.params
    logs, check_logs =
        _amp_pair(data, MST_AMP_NMAX, MST_AMP_CHECK_NMAX)
    teuk = _amp_teuk_pair(logs, branch)
    gsn = _amp_gsn_pair(base, branch, teuk)
    check_teuk = _amp_teuk_pair(check_logs, branch)
    check_gsn = _amp_gsn_pair(base, branch, check_teuk)
    truncation = Float64(_amp_distance(gsn, check_gsn))
    moved = with_nu(refinement.nu - refinement.step)
    moved_gsn = _amp_gsn_pair(moved.params, branch, _amp_teuk_pair(
        _amp_logs(moved, MST_AMP_NMAX), branch))
    component_agreement = _amp_component_distance(gsn, check_gsn) .+
        _amp_component_distance(gsn, moved_gsn)
    condition = Float64(max(
        logs.max_condition,
        check_logs.max_condition,
    ))
    all(_finite_complex, (teuk..., gsn..., check_teuk..., check_gsn...)) ||
        throw(MSTCertificateError(
            :nia_amplitude_nonfinite,
            "nonfinite principal-nu analytic MST amplitude.",
        ))
    truncation <= MST_AMP_TRUNCATION_MAX ||
        throw(MSTCertificateError(
            :nia_amplitude_certificate,
            "principal-nu analytic MST truncation rejected: " *
            "agreement=$(truncation).",
        ))
    condition <= MST_AMP_CONDITION_MAX ||
        throw(MSTCertificateError(
            :nia_amplitude_certificate,
            "principal-nu analytic MST condition rejected: " *
            "condition=$(condition).",
        ))
    return (
        shift=round(Int, real(base.nu)),
        teuk,
        gsn,
        transmission_log=branch == :IN ?
            ComplexF64(logs.in_trans) : ComplexF64(logs.up_trans),
        check_transmission_log=branch == :IN ?
            ComplexF64(check_logs.in_trans) :
            ComplexF64(check_logs.up_trans),
        medoid_score=0.0,
        representation_spread=0.0,
        nearest_agreement=truncation,
        truncation_agreement=truncation,
        component_agreement,
        max_condition=condition,
        gsn_error=_amp_ratio_errors(logs, branch),
        nmax=MST_AMP_NMAX,
        check_nmax=MST_AMP_CHECK_NMAX,
        certificate_kind=:principal_single_nu,
        certificate_accepted=true,
    )
end

# Amplitudes as in _principal_amplitudes, with nu refined and the sums evaluated in
# double-double.
@noinline function _refined_amplitudes(
    data::MSTSeriesData,
    branch::Symbol,
)
    base = data.params
    offset, nu_residual = _dd_refine(base)
    nu = ComplexF64(base.l) + dc_value(offset)
    refined = MSTSeriesData(MSTParams(
        base.s,
        base.l,
        base.m,
        base.a,
        base.omega,
        base.lambda,
        nu,
        dc_value(offset),
    ))
    logs, check_logs = _dd_amp_pair(refined, offset)
    teuk = _amp_teuk_pair(logs, branch)
    gsn = _amp_gsn_pair(refined.params, branch, teuk)
    check_teuk = _amp_teuk_pair(check_logs, branch)
    check_gsn = _amp_gsn_pair(refined.params, branch, check_teuk)
    truncation = Float64(_amp_distance(gsn, check_gsn))
    component_agreement = _amp_component_distance(gsn, check_gsn)
    condition = Float64(max(
        logs.max_condition,
        check_logs.max_condition,
    ))
    all(_finite_complex, (
        teuk...,
        gsn...,
        check_teuk...,
        check_gsn...,
    )) || throw(MSTCertificateError(
        :dd_amplitude_nonfinite,
        "nonfinite DD-refined MST amplitude.",
    ))
    truncation <= MST_AMP_TRUNCATION_MAX ||
        throw(MSTCertificateError(
            :dd_amplitude_truncation,
            "DD-refined MST amplitude truncation rejected: " *
            "agreement=$truncation.",
        ))
    condition <= MST_DD_CONDITION_MAX ||
        throw(MSTCertificateError(
            :dd_amplitude_condition,
            "DD-refined MST amplitude condition rejected: " *
            "condition=$condition.",
        ))
    return (
        shift=round(Int, real(nu)),
        nu,
        nu_offset_dd=offset,
        nu_offset_low=offset.low,
        nu_residual,
        teuk,
        gsn,
        transmission_log=branch == :IN ?
            ComplexF64(logs.in_trans) : ComplexF64(logs.up_trans),
        check_transmission_log=branch == :IN ?
            ComplexF64(check_logs.in_trans) :
            ComplexF64(check_logs.up_trans),
        medoid_score=0.0,
        representation_spread=0.0,
        nearest_agreement=truncation,
        truncation_agreement=truncation,
        component_agreement,
        max_condition=condition,
        gsn_error=_amp_ratio_errors(logs, branch),
        nmax=MST_AMP_NMAX,
        check_nmax=MST_AMP_CHECK_NMAX,
        certificate_kind=:principal_dd_refined,
        certificate_accepted=true,
    )
end

"""
    mst_principal_amplitudes(coefficients_or_params, branch)

Analytic MST asymptotic amplitudes of `branch` (`:IN` or `:UP`) at the principal nu of
the mode, certified by a truncation and a condition check; throws `MSTCertificateError`
otherwise.
"""
@noinline mst_principal_amplitudes(
    coefficients::DirectCoefficientSet,
    branch::Symbol,
) = _principal_amplitudes(_mst_data(coefficients), branch)

@noinline mst_principal_amplitudes(
    params::DirectGSNParameters,
    branch::Symbol,
) = _principal_amplitudes(_mst_data(params), branch)

# Branch strength and amplitude ratio on the negative imaginary axis from given series
# data, at one truncation order.
function _nia_branch_strength(
    data::MSTSeriesData;
    nmax::Int=MST_NIA_SUM_MAX,
)
    p = data.params
    omega = ComplexF64(p.omega)
    iszero(real(omega)) && imag(omega) < 0 ||
        throw(ArgumentError(
            "MST NIA branch strength requires omega=-i*sigma with sigma>0."))

    fsum = _nia_full_logsum(
        n -> log_series_coefficient!(data, n), nmax)
    aminus = _nia_aminus_sum(data, nmax)
    nu = _nu_shift(p)
    epsc = ComplexF64(p.epsilon)

    log_ratio =
        (2 * p.s - 2I * epsc) * log(2.0) +
        I * pi * (nu + 1) +
        loggamma(1 - p.s + I * epsc + nu) -
        loggamma(1 + p.s - I * epsc + nu) +
        _lognorm_value(fsum) -
        _lognorm_value(aminus)
    monodromy = -expm1(2pi * (epsc - I * nu))
    _finite_complex(monodromy) && !iszero(monodromy) ||
        throw(MSTCertificateError(
            :nia_monodromy,
            "nonfinite or zero MST NIA monodromy factor.",
        ))
    log_strength =
        I * pi / 2 +
        log_ratio +
        2 * p.s * log(omega) -
        2I * epsc * log(epsc) +
        I * epsc * (1 - p.kappa) +
        log(monodromy)
    ratio = ComplexF64(exp(log_ratio))
    strength = ComplexF64(exp(log_strength))
    all(_finite_complex, (ratio, strength)) ||
        throw(MSTCertificateError(
            :nia_branch_strength,
            "nonfinite MST NIA branch strength.",
        ))
    return (
        nu=ComplexF64(nu),
        amplitude_ratio=ratio,
        branch_strength=strength,
        estimated_relerr=max(
            fsum.estimated_relerr,
            aminus.estimated_relerr,
        ),
        coefficient_range=(
            min(fsum.nmin, aminus.nmin),
            max(fsum.nmax, aminus.nmax),
        ),
    )
end

"""
    direct_mst_nia_branch_strength(s, l, m, a, omega, lambda; nu=nothing)

Strength of the branch cut of the MST solution on the negative imaginary frequency axis,
omega = -i sigma with sigma > 0, from series coefficients computed by double-double
backward recurrence. The estimated error combines the sum error, the change between two
truncation orders and the change under one-ulp perturbations of a, omega, lambda and nu;
`certificate_accepted` is true when it is at most `MST_NIA_CERT_TOL`.
"""
function direct_mst_nia_branch_strength(
    s::Integer,
    l::Integer,
    m::Integer,
    a::Real,
    omega::Number,
    lambda;
    nu=nothing,
)
    omegac = ComplexF64(omega)
    nu_value = nu === nothing ?
        mst_nu(s, l, m, a, omegac, lambda) :
        ComplexF64(nu)
    function strength_at(av, omegav, lambdav, nuv, cutoff, nmax)
        params = MSTParams(
            s, l, m, av, omegav, lambdav, nuv,
            ComplexF64(nuv) - ComplexF64(l))
        return _nia_branch_strength(
            _nia_series_data(params, cutoff); nmax)
    end
    fine = strength_at(
        a, omegac, lambda, nu_value, MST_NIA_TAIL, MST_NIA_SUM_MAX)
    coarse = strength_at(
        a, omegac, lambda, nu_value,
        MST_NIA_CHECK_TAIL, MST_NIA_CHECK_SUM_MAX)
    scale = max(
        abs(fine.branch_strength),
        abs(coarse.branch_strength),
        floatmin(Float64),
    )
    adjacent_error =
        abs(fine.branch_strength - coarse.branch_strength) / scale

    af = Float64(a)
    lambdac = ComplexF64(lambda)
    perturbations = (
        (nextfloat(af), omegac, lambdac, nu_value),
        (prevfloat(af), omegac, lambdac, nu_value),
        (af, ComplexF64(0, nextfloat(imag(omegac))), lambdac, nu_value),
        (af, ComplexF64(0, prevfloat(imag(omegac))), lambdac, nu_value),
        (af, omegac,
            ComplexF64(nextfloat(real(lambdac)), imag(lambdac)), nu_value),
        (af, omegac,
            ComplexF64(prevfloat(real(lambdac)), imag(lambdac)), nu_value),
        (af, omegac,
            ComplexF64(real(lambdac), nextfloat(imag(lambdac))), nu_value),
        (af, omegac,
            ComplexF64(real(lambdac), prevfloat(imag(lambdac))), nu_value),
        (af, omegac, lambdac,
            ComplexF64(nextfloat(real(nu_value)), imag(nu_value))),
        (af, omegac, lambdac,
            ComplexF64(prevfloat(real(nu_value)), imag(nu_value))),
        (af, omegac, lambdac,
            ComplexF64(real(nu_value), nextfloat(imag(nu_value)))),
        (af, omegac, lambdac,
            ComplexF64(real(nu_value), prevfloat(imag(nu_value)))),
    )
    input_condition_error = 0.0
    for (av, omegav, lambdav, nuv) in perturbations
        abs(av) < 1.0 || continue
        candidate = try
            strength_at(
                av, omegav, lambdav, nuv,
                MST_NIA_CHECK_TAIL, MST_NIA_CHECK_SUM_MAX)
        catch error
            error isa MSTCertificateError || rethrow()
            input_condition_error = Inf
            break
        end
        candidate_scale = max(
            abs(fine.branch_strength),
            abs(candidate.branch_strength),
            floatmin(Float64),
        )
        input_condition_error = max(
            input_condition_error,
            abs(fine.branch_strength - candidate.branch_strength) /
                candidate_scale,
        )
    end
    estimated_relerr = max(
        fine.estimated_relerr,
        adjacent_error,
        input_condition_error,
    )
    return merge(fine, (
        adjacent_error=Float64(adjacent_error),
        input_condition_error=Float64(input_condition_error),
        estimated_relerr=Float64(estimated_relerr),
        certificate_accepted=estimated_relerr <= MST_NIA_CERT_TOL,
    ))
end

"""
    nia_strength_fast(s, l, m, a, omega, lambda; nu=nothing)

Cheaper `direct_mst_nia_branch_strength`: lower recurrence cutoffs and sum orders, and
an error estimate from two truncation orders only.
"""
function nia_strength_fast(
    s::Integer,
    l::Integer,
    m::Integer,
    a::Real,
    omega::Number,
    lambda;
    nu=nothing,
)
    omegac = ComplexF64(omega)
    nu_value = nu === nothing ?
        mst_nu(s, l, m, a, omegac, lambda) :
        ComplexF64(nu)
    params = MSTParams(
        s, l, m, a, omegac, lambda, nu_value,
        ComplexF64(nu_value) - ComplexF64(l))
    fine = _nia_branch_strength(
        _nia_series_data(params, MST_NIA_FAST_TAIL);
        nmax=MST_NIA_FAST_SUM_MAX,
    )
    coarse = _nia_branch_strength(
        _nia_series_data(params, MST_NIA_FAST_CHECK_TAIL);
        nmax=MST_NIA_FAST_CHECK_SUM_MAX,
    )
    scale = max(
        abs(fine.branch_strength),
        abs(coarse.branch_strength),
        floatmin(Float64),
    )
    adjacent_error =
        abs(fine.branch_strength - coarse.branch_strength) / scale
    estimated_relerr = max(
        fine.estimated_relerr,
        coarse.estimated_relerr,
        adjacent_error,
    )
    return merge(fine, (
        adjacent_error=Float64(adjacent_error),
        estimated_relerr=Float64(estimated_relerr),
        certificate_accepted=estimated_relerr <= MST_NIA_CERT_TOL,
    ))
end

"""
    direct_mst_nia_jump(s, l, m, a, omega, lambda, unit_incidence, unit_reflection; nu=nothing, require_certificate=true)

Jump coefficient across the branch cut on the negative imaginary frequency axis at
omega = -i sigma (sigma > 0), from the branch strength q of
`direct_mst_nia_branch_strength` and the unit-transmission incidence and reflection
amplitudes. Wronskians and the jump are carried in log-polar form. With
`require_certificate`, an uncertified strength throws `MSTCertificateError`.
"""
function direct_mst_nia_jump(
    s::Integer,
    l::Integer,
    m::Integer,
    a::Real,
    omega::Number,
    lambda,
    unit_incidence,
    unit_reflection;
    nu=nothing,
    require_certificate::Bool=true,
)
    omegac = ComplexF64(omega)
    iszero(real(omegac)) && imag(omegac) < 0 ||
        throw(ArgumentError(
            "MST NIA jump requires omega=-i*sigma with sigma>0."))
    incidence = ComplexF64(unit_incidence)
    reflection = ComplexF64(unit_reflection)
    all(_finite_complex, (incidence, reflection)) ||
        throw(ArgumentError(
            "MST NIA jump requires finite unit-transmission amplitudes."))

    strength = direct_mst_nia_branch_strength(
        s, l, m, a, omegac, lambda;
        nu,
    )
    require_certificate && !strength.certificate_accepted &&
        throw(MSTCertificateError(
            :nia_condition,
            "MST NIA jump rejected: estimated relative error=" *
            "$(strength.estimated_relerr), input condition=" *
            "$(strength.input_condition_error).",
        ))
    sigma = -imag(omegac)
    q = strength.branch_strength

    function logpolar(value::ComplexF64)
        magnitude = abs(value)
        isfinite(magnitude) && !iszero(magnitude) ||
            throw(MSTCertificateError(
                :nia_logpolar,
                "zero or nonfinite MST NIA log-polar value.",
            ))
        return (
            phase=ComplexF64(value / magnitude),
            logabs=Float64(log(magnitude)),
        )
    end
    multiply(left, right) = (
        phase=ComplexF64(left.phase * right.phase),
        logabs=Float64(left.logabs + right.logabs),
    )
    function add_logpolar(left, right)
        scale = max(left.logabs, right.logabs)
        total = left.phase * exp(left.logabs - scale) +
            right.phase * exp(right.logabs - scale)
        magnitude = abs(total)
        isfinite(magnitude) && !iszero(magnitude) ||
            throw(MSTCertificateError(
                :nia_wronskian_cancellation,
                "zero or nonfinite MST NIA Wronskian sum.",
            ))
        return (
            phase=ComplexF64(total / magnitude),
            logabs=Float64(scale + log(magnitude)),
        )
    end
    function materialize(value)
        value.logabs > LOG_FLOAT64_MAX && return nothing
        value.logabs < LOG_FLOAT64_MIN && return ComplexF64(0)
        result = ComplexF64(exp(value.logabs) * value.phase)
        _finite_complex(result) ? result : nothing
    end

    positive_scale = (
        phase=ComplexF64(1),
        logabs=Float64(log(2sigma)),
    )
    wronskian_plus_scaled = multiply(
        positive_scale, logpolar(incidence))
    wronskian_minus_scaled = if iszero(reflection)
        wronskian_plus_scaled
    else
        correction_scaled = multiply(
            positive_scale,
            multiply(
                (phase=ComplexF64(-I), logabs=0.0),
                multiply(logpolar(q), logpolar(reflection)),
            ),
        )
        add_logpolar(wronskian_plus_scaled, correction_scaled)
    end
    wronskian_product_scaled = multiply(
        wronskian_plus_scaled, wronskian_minus_scaled)
    numerator_scaled = multiply(
        positive_scale,
        multiply(
            (phase=ComplexF64(-I), logabs=0.0),
            logpolar(q),
        ),
    )
    green_scaled = (
        phase=ComplexF64(
            numerator_scaled.phase /
            wronskian_product_scaled.phase),
        logabs=Float64(
            numerator_scaled.logabs -
            wronskian_product_scaled.logabs),
    )
    green_jump = materialize(green_scaled)
    green_jump === nothing &&
        throw(MSTCertificateError(
            :nia_jump_range,
            "MST NIA jump lies outside the ComplexF64 range.",
        ))
    wronskian_plus = materialize(wronskian_plus_scaled)
    wronskian_minus = materialize(wronskian_minus_scaled)
    wronskian_product = materialize(wronskian_product_scaled)
    return merge(strength, (
        sigma=Float64(sigma),
        up_jump_coefficient=ComplexF64(I * q),
        wronskian_plus,
        wronskian_minus,
        wronskian_product,
        wronskian_plus_scaled,
        wronskian_minus_scaled,
        wronskian_product_scaled,
        green_jump_scaled=green_scaled,
        green_jump_coefficient=ComplexF64(green_jump),
    ))
end

# Log-norm normalization of the outgoing branch from the A_- sum and its prefactor.
function _outgoing_raw_factor(data::MSTSeriesData)
    p = data.params
    fsum = _sum_symmetric_logs() do n
        nc = ComplexF64(n)
        logc = log_series_coefficient!(data, n)
        logfac = log_pochhammer(_nu_shift(p, 1 + p.s) - I * p.epsilon, n) -
            log_pochhammer(_nu_shift(p, 1 - p.s) + I * p.epsilon, n)
        phase = I * pi * (_npnu(p, n) + 1 + p.s - I * p.epsilon)
        logc + logfac + phase
    end
    epsv = ComplexF64(p.epsilon)
    kappa = ComplexF64(p.kappa)
    nu = _nu_shift(p)
    s = ComplexF64(p.s)
    log_prefactor = nu * log(2.0) - pi * epsv - I * pi * (nu + 1 + s) +
        I * epsv * kappa + (nu - s) * log(epsv * kappa) -
        (nu + 1 + s - I * epsv) * log(2I * epsv * kappa)
    return LogNormSum(
        ComplexF64(fsum.phase * cis(imag(log_prefactor))),
        Float64(fsum.logabs + real(log_prefactor)),
        fsum.estimated_relerr,
        fsum.nmin,
        fsum.nmax,
    )
end

# The same normalization from fixed-order A_- sums (full and check order); throws when
# truncation or conditioning exceed `tolerance`.
function _outgoing_raw_fast(data::MSTSeriesData, tolerance::Float64)
    full, check = _amp_aminus_pair(
        data, MST_AMP_NMAX, MST_AMP_CHECK_NMAX)
    truncation = abs(exp(check.logvalue - full.logvalue) - 1)
    condition = max(full.condition, check.condition)
    estimated_relerr = max(truncation, eps(Float64) * condition)
    estimated_relerr <= tolerance || throw(MSTCertificateError(
        :outgoing_factor,
        "MST outgoing normalization rejected: truncation=$(truncation), " *
        "condition=$(condition).",
    ))
    p = data.params
    epsv = ComplexF64(p.epsilon)
    kappa = ComplexF64(p.kappa)
    nu = _nu_shift(p)
    s = ComplexF64(p.s)
    sum_phase = 1im * pi * (nu + 1 + s - 1im * epsv)
    log_prefactor = nu * log(2.0) - pi * epsv -
        1im * pi * (nu + 1 + s) +
        1im * epsv * kappa + (nu - s) * log(epsv * kappa) -
        (nu + 1 + s - 1im * epsv) * log(2im * epsv * kappa)
    value = ComplexF64(full.logvalue + sum_phase + log_prefactor)
    return LogNormSum(
        ComplexF64(cis(imag(value))),
        Float64(real(value)),
        Float64(estimated_relerr),
        -MST_AMP_NMAX,
        MST_AMP_NMAX,
    )
end

# Log factor taking the raw MST series of `branch` to unit normalization: the raw sum
# times f3 or f4 over the conversion amplitude Binc or Bref.
function _raw_unit_logfactor(data::MSTSeriesData, params, branch::Symbol, raw=nothing)
    cf = _conversion_module()
    f3, f4 = _legacy_factors(params)
    raw, legacy, conversion = if branch == :in
        raw = raw === nothing ? _incoming_raw_factor(data) : raw
        raw, f3, ComplexF64(getfield(cf, :Binc)(
            params.s, params.m, params.a, params.omega, params.lambda))
    elseif branch == :out
        raw = raw === nothing ? _outgoing_raw_factor(data) : raw
        raw, f4, ComplexF64(getfield(cf, :Bref)(
            params.s, params.m, params.a, params.omega, params.lambda))
    else
        throw(ArgumentError("MST infinity branch must be :in or :out."))
    end
    _finite_complex(raw.phase) && isfinite(raw.logabs) &&
        _finite_complex(legacy) && !iszero(legacy) &&
        _finite_complex(conversion) && !iszero(conversion) ||
        error("direct GSN MST raw-to-unit log factor is nonfinite.")
    logfactor = raw.logabs + log(legacy) - log(conversion) +
        I * angle(raw.phase)
    return (
        phase=ComplexF64(cis(imag(logfactor))),
        logabs=Float64(real(logfactor)),
        estimated_relerr=raw.estimated_relerr,
    )
end

# n-th term of the hypergeometric (Type-1) physical-in series and its first two
# x-derivatives, weighted by a_n / norm.
function _type1_triplet(data::MSTSeriesData, n::Integer, x, norm::LogNormSum)
    p = data.params
    npnu = _npnu(p, n)
    a = npnu + 1 - 1im * p.tau
    b = -npnu - 1im * p.tau
    c = 1 - p.s - 1im * (p.epsilon + p.tau)
    xc = ComplexF64(x)
    value = ComplexF64(HGF_2F1(a, b, c, xc))
    deriv = ComplexF64(a * b / c * HGF_2F1(a + 1, b + 1, c + 1, xc))
    second = ComplexF64(
        a * (a + 1) * b * (b + 1) / (c * (c + 1)) *
        HGF_2F1(a + 2, b + 2, c + 2, xc),
    )
    log_weight = log_series_coefficient!(data, n) - norm.logabs -
        1im * angle(norm.phase)
    return MSTTriplet(
        _scaled_product(log_weight, value),
        _scaled_product(log_weight, deriv),
        _scaled_product(log_weight, second),
    )
end

function _coefficient_error(data::MSTSeriesData, n::Integer)
    return lock(data.coefficient_lock) do
        get(data.coeff_errors, Int(n), Inf)
    end
end

@inline _type1_error(data::MSTSeriesData, n::Integer) =
    min(Inf, _coefficient_error(data, n) + 8eps(Float64))

const MST_DD_STATE_CHECK_NMAX = 160
const MST_DD_REG_NMAX = 64
const MST_DD_REG_CHECK_NMAX = 62
const MST_DD_REG_SERIES_MAX = 4096
const MST_DD_REG_TOL = 8eps(Float64)
const MST_DD_REG_POLE_RADIUS = 0.125

# Double-double n-th Type-1 term: coefficient times
# 2F1(n+nu+1-i tau, -n-nu-i tau; 1-s-i(epsilon+tau); x) and its first two derivatives.
function _dd_type1_term(
    p::MSTParams,
    n::Int,
    x,
    coefficient::DDComplex,
)
    npnu = _npnu(p, n)
    a = npnu + 1 - I * p.tau
    b = -npnu - I * p.tau
    c = 1 - p.s - I * (p.epsilon + p.tau)
    xc = ComplexF64(x)
    value = ComplexF64(HGF_2F1(a, b, c, xc))
    deriv = ComplexF64(
        a * b / c * HGF_2F1(a + 1, b + 1, c + 1, xc))
    second = ComplexF64(
        a * (a + 1) * b * (b + 1) / (c * (c + 1)) *
        HGF_2F1(a + 2, b + 2, c + 2, xc),
    )
    all(_finite_complex, (value, deriv, second)) ||
        throw(MSTCertificateError(
            :dd_physical_in_term,
            "nonfinite DD physical-in basis term at n=$n.",
        ))
    return (
        dc_mul(coefficient, DDComplex(value)),
        dc_mul(coefficient, DDComplex(deriv)),
        dc_mul(coefficient, DDComplex(second)),
    )
end

# Double-double physical-in sums normalized by sum a_n, at order nmax and check_nmax.
function _dd_type1_pair(
    data::MSTSeriesData,
    offset::DDComplex,
    x;
    nmax::Int=MST_RESIDUAL_NMAX,
    check_nmax::Int=MST_DD_STATE_CHECK_NMAX,
)
    0 <= check_nmax <= nmax ||
        throw(ArgumentError("invalid DD physical-in check order."))
    p = data.params
    positive_ratios, negative_ratios = _dd_ratios(
        p, offset; nmax)
    positive_coefficients = Vector{DDComplex}(undef, nmax)
    negative_coefficients = Vector{DDComplex}(undef, nmax)
    coefficient = DDComplex(1)
    for n in 1:nmax
        coefficient = dc_mul(coefficient, positive_ratios[n])
        positive_coefficients[n] = coefficient
    end
    coefficient = DDComplex(1)
    for n in 1:nmax
        coefficient = dc_mul(coefficient, negative_ratios[n])
        negative_coefficients[n] = coefficient
    end

    zero_term = _dd_type1_term(p, 0, x, DDComplex(1))
    positive_terms = Vector{NTuple{3,DDComplex}}(undef, nmax)
    negative_terms = Vector{NTuple{3,DDComplex}}(undef, nmax)
    for n in 1:nmax
        positive_terms[n] = _dd_type1_term(
            p, n, x, positive_coefficients[n])
        negative_terms[n] = _dd_type1_term(
            p, -n, x, negative_coefficients[n])
    end

    function sum_at(limit)
        norm = DDComplex(1)
        sums = [zero_term[1], zero_term[2], zero_term[3]]
        absolute = [
            abs(dc_value(zero_term[1])),
            abs(dc_value(zero_term[2])),
            abs(dc_value(zero_term[3])),
        ]
        norm_absolute = 1.0
        for n in 1:limit
            norm = dc_add(norm, positive_coefficients[n])
            norm = dc_add(norm, negative_coefficients[n])
            norm_absolute +=
                abs(dc_value(positive_coefficients[n])) +
                abs(dc_value(negative_coefficients[n]))
            for component in 1:3
                sums[component] = dc_add(
                    sums[component], positive_terms[n][component])
                sums[component] = dc_add(
                    sums[component], negative_terms[n][component])
                absolute[component] +=
                    abs(dc_value(positive_terms[n][component])) +
                    abs(dc_value(negative_terms[n][component]))
            end
        end
        norm_value = dc_value(norm)
        _finite_complex(norm_value) && !iszero(norm_value) ||
            throw(MSTCertificateError(
                :dd_physical_in_norm,
                "nonfinite DD physical-in normalization.",
            ))
        values = ntuple(component -> begin
            value = dc_value(sums[component])
            _finite_complex(value) && !iszero(value) ||
                throw(MSTCertificateError(
                    :dd_physical_in_sum,
                    "nonfinite DD physical-in state sum.",
                ))
            ComplexF64(dc_value(dc_div(sums[component], norm)))
        end, 3)
        condition = maximum((
            norm_absolute / abs(norm_value),
            ntuple(component ->
                absolute[component] /
                abs(dc_value(sums[component])), 3)...,
        ))
        return (
            state=MSTTriplet(values...),
            condition=Float64(condition),
        )
    end
    return sum_at(nmax), sum_at(check_nmax)
end

@inline _dd_shift(value::DDComplex, shift) =
    dc_add(value, DDComplex(shift))

# Hypergeometric parameter c = 1 - s - i(epsilon + tau) of the Type-1 series.
@inline function _dd_hyper_c(p::MSTParams)
    return dc_sub(
        DDComplex(1 - p.s),
        dc_imul(dc_add(DDComplex(p.epsilon), DDComplex(p.tau))),
    )
end

# True when c is within MST_DD_REG_POLE_RADIUS of a nonpositive integer, where the
# regularized Type-1 series is used.
@inline function _dd_regime(p::MSTParams)
    c = dc_value(_dd_hyper_c(p))
    nearest = round(Int, -real(c))
    return nearest >= 0 &&
        real(c) <= 0.5 &&
        abs(c + nearest) <= MST_DD_REG_POLE_RADIUS
end

# Double-double power series of the regularized 2F1(a, b; c; z) / Gamma(c). Near a
# nonpositive integer c the series starts past the pole and recurs backward to k = 0.
function _dd_reg_series(
    a::DDComplex,
    b::DDComplex,
    c::DDComplex,
    z::DDComplex,
)
    c_value = dc_value(c)
    nearest = max(0, round(Int, -real(c_value)))
    near_pole = real(c_value) <= 0.5 &&
        abs(c_value + nearest) <= 0.5
    start = near_pole ? nearest + 1 : 0

    term = DDComplex(exp(-loggamma(ComplexF64(c_value + start))))
    for k in 0:(start - 1)
        numerator = dc_mul(
            dc_mul(_dd_shift(a, k), _dd_shift(b, k)),
            z,
        )
        term = dc_mul(
            term,
            dc_div(numerator, DDComplex(k + 1)),
        )
    end
    _finite_complex(dc_value(term)) ||
        throw(MSTCertificateError(
            :dd_regularized_type1,
            "nonfinite regularized Type-1 start term.",
        ))

    lower = Vector{DDComplex}(undef, start + 1)
    lower[start + 1] = term
    for k in (start - 1):-1:0
        denominator = dc_mul(
            dc_mul(_dd_shift(a, k), _dd_shift(b, k)),
            z,
        )
        iszero(dc_value(denominator)) &&
            throw(MSTCertificateError(
                :dd_regularized_type1,
                "singular regularized Type-1 backward recurrence.",
            ))
        numerator = dc_mul(
            term,
            dc_mul(_dd_shift(c, k), DDComplex(k + 1)),
        )
        term = dc_div(numerator, denominator)
        _finite_complex(dc_value(term)) ||
            throw(MSTCertificateError(
                :dd_regularized_type1,
                "nonfinite regularized Type-1 backward term.",
            ))
        lower[k + 1] = term
    end

    total = DDComplex(0)
    absolute_sum = 0.0
    for value in lower
        total = dc_add(total, value)
        absolute_sum += abs(dc_value(value))
    end

    term = lower[end]
    converged = false
    tail_bound = Inf
    iterations = start
    for k in start:(MST_DD_REG_SERIES_MAX - 1)
        numerator = dc_mul(
            dc_mul(_dd_shift(a, k), _dd_shift(b, k)),
            z,
        )
        denominator = dc_mul(
            _dd_shift(c, k),
            DDComplex(k + 1),
        )
        iszero(dc_value(denominator)) &&
            throw(MSTCertificateError(
                :dd_regularized_type1,
                "singular regularized Type-1 forward recurrence.",
            ))
        ratio = dc_div(numerator, denominator)
        term = dc_mul(term, ratio)
        term_value = dc_value(term)
        _finite_complex(term_value) ||
            throw(MSTCertificateError(
                :dd_regularized_type1,
                "nonfinite regularized Type-1 forward term.",
            ))
        total = dc_add(total, term)
        absolute_sum += abs(term_value)
        iterations = k + 1
        ratio_abs = abs(dc_value(ratio))
        total_abs = max(abs(dc_value(total)), floatmin(Float64))
        if iterations >= start + 24 && ratio_abs < 0.9
            tail_bound = abs(term_value) /
                max(1 - ratio_abs, eps(Float64))
            if tail_bound <= MST_DD_REG_TOL * total_abs
                converged = true
                break
            end
        end
    end
    converged || throw(MSTCertificateError(
        :dd_regularized_type1,
        "regularized Type-1 series did not converge.",
    ))
    total_value = dc_value(total)
    scale = max(abs(total_value), floatmin(Float64))
    return (
        value=total,
        condition=Float64(absolute_sum / scale),
        relative_tail=Float64(tail_bound / scale),
        iterations,
    )
end

# Regularized 2F1(a, b; c; z) / Gamma(c) through the Pfaff transformation
# (1 - z)^(-a) 2F1(a, c - b; c; z / (z - 1)), for |z / (z - 1)| <= 1/2.
function _dd_reg_2f1(
    a::DDComplex,
    b::DDComplex,
    c::DDComplex,
    z::DDComplex,
)
    z_value = dc_value(z)
    transformed_value = ComplexF64(z_value / (z_value - 1))
    abs(transformed_value) <= 0.5 + 8eps(Float64) ||
        throw(MSTCertificateError(
            :dd_regularized_type1,
            "regularized Type-1 Pfaff argument is outside its domain.",
        ))
    series = _dd_reg_series(
        a,
        dc_sub(c, b),
        c,
        DDComplex(transformed_value),
    )
    logarithm = ComplexF64(log1p(-z_value))
    exponent_high = -a.high * logarithm
    exponent_low = -a.low * logarithm
    prefactor_high = ComplexF64(exp(exponent_high))
    prefactor = DDComplex(
        prefactor_high,
        prefactor_high * exponent_low,
    )
    value = dc_mul(prefactor, series.value)
    _finite_complex(dc_value(value)) ||
        throw(MSTCertificateError(
            :dd_regularized_type1,
            "nonfinite regularized Type-1 value.",
        ))
    return merge(series, (; value))
end

# Regularized 2F1 of the n-th Type-1 term and its first two x-derivatives.
function _dd_reg_triplet(
    p::MSTParams,
    offset::DDComplex,
    n::Int,
    x,
)
    npnu = _dd_shift(
        dc_add(DDComplex(p.l), offset),
        n,
    )
    tau = DDComplex(p.tau)
    a = dc_sub(_dd_shift(npnu, 1), dc_imul(tau))
    b = dc_sub(dc_neg(npnu), dc_imul(tau))
    c = _dd_hyper_c(p)
    z = DDComplex(x)
    value = _dd_reg_2f1(a, b, c, z)
    derivative = _dd_reg_2f1(
        _dd_shift(a, 1),
        _dd_shift(b, 1),
        _dd_shift(c, 1),
        z,
    )
    second = _dd_reg_2f1(
        _dd_shift(a, 2),
        _dd_shift(b, 2),
        _dd_shift(c, 2),
        z,
    )
    derivative_factor = dc_mul(a, b)
    second_factor = dc_mul(
        dc_mul(a, _dd_shift(a, 1)),
        dc_mul(b, _dd_shift(b, 1)),
    )
    return (
        values=(
            value.value,
            dc_mul(derivative_factor, derivative.value),
            dc_mul(second_factor, second.value),
        ),
        condition=max(
            value.condition,
            derivative.condition,
            second.condition,
        ),
        relative_tail=max(
            value.relative_tail,
            derivative.relative_tail,
            second.relative_tail,
        ),
        iterations=max(
            value.iterations,
            derivative.iterations,
            second.iterations,
        ),
    )
end

# Regularized Type-1 physical-in sums times Gamma(c), at order nmax and check_nmax.
function _dd_reg_pair(
    data::MSTSeriesData,
    offset::DDComplex,
    x;
    nmax::Int=MST_DD_REG_NMAX,
    check_nmax::Int=MST_DD_REG_CHECK_NMAX,
)
    0 <= check_nmax <= nmax ||
        throw(ArgumentError("invalid regularized Type-1 check order."))
    p = data.params
    positive_ratios, negative_ratios = _dd_ratios(
        p, offset; nmax)
    positive_coefficients = Vector{DDComplex}(undef, nmax)
    negative_coefficients = Vector{DDComplex}(undef, nmax)
    coefficient = DDComplex(1)
    for n in 1:nmax
        coefficient = dc_mul(coefficient, positive_ratios[n])
        positive_coefficients[n] = coefficient
    end
    coefficient = DDComplex(1)
    for n in 1:nmax
        coefficient = dc_mul(coefficient, negative_ratios[n])
        negative_coefficients[n] = coefficient
    end

    zero_basis = _dd_reg_triplet(p, offset, 0, x)
    zero_term = zero_basis.values
    positive_terms = Vector{NTuple{3,DDComplex}}(undef, nmax)
    negative_terms = Vector{NTuple{3,DDComplex}}(undef, nmax)
    basis_condition = zero_basis.condition
    relative_tail = zero_basis.relative_tail
    for n in 1:nmax
        positive_basis = _dd_reg_triplet(p, offset, n, x)
        negative_basis = _dd_reg_triplet(p, offset, -n, x)
        positive_terms[n] = ntuple(component ->
            dc_mul(
                positive_coefficients[n],
                positive_basis.values[component],
            ), 3)
        negative_terms[n] = ntuple(component ->
            dc_mul(
                negative_coefficients[n],
                negative_basis.values[component],
            ), 3)
        basis_condition = max(
            basis_condition,
            positive_basis.condition,
            negative_basis.condition,
        )
        relative_tail = max(
            relative_tail,
            positive_basis.relative_tail,
            negative_basis.relative_tail,
        )
    end
    gamma_c = DDComplex(gamma(dc_value(_dd_hyper_c(p))))

    function sum_at(limit)
        norm = DDComplex(1)
        sums = [zero_term[1], zero_term[2], zero_term[3]]
        absolute = [
            abs(dc_value(zero_term[1])),
            abs(dc_value(zero_term[2])),
            abs(dc_value(zero_term[3])),
        ]
        norm_absolute = 1.0
        for n in 1:limit
            norm = dc_add(norm, positive_coefficients[n])
            norm = dc_add(norm, negative_coefficients[n])
            norm_absolute +=
                abs(dc_value(positive_coefficients[n])) +
                abs(dc_value(negative_coefficients[n]))
            for component in 1:3
                sums[component] = dc_add(
                    sums[component],
                    positive_terms[n][component],
                )
                sums[component] = dc_add(
                    sums[component],
                    negative_terms[n][component],
                )
                absolute[component] +=
                    abs(dc_value(positive_terms[n][component])) +
                    abs(dc_value(negative_terms[n][component]))
            end
        end
        norm_value = dc_value(norm)
        _finite_complex(norm_value) && !iszero(norm_value) ||
            throw(MSTCertificateError(
                :dd_regularized_type1,
                "nonfinite regularized Type-1 normalization.",
            ))
        dd_state = ntuple(component ->
            dc_mul(
                gamma_c,
                dc_div(sums[component], norm),
            ), 3)
        values = ntuple(component -> begin
            value = dc_value(dd_state[component])
            _finite_complex(value) ||
                throw(MSTCertificateError(
                    :dd_regularized_type1,
                    "nonfinite regularized Type-1 state sum.",
                ))
            ComplexF64(value)
        end, 3)
        condition = maximum((
            basis_condition,
            norm_absolute / abs(norm_value),
            ntuple(component ->
                absolute[component] /
                max(abs(dc_value(sums[component])),
                    floatmin(Float64)), 3)...,
        ))
        return (
            state=MSTTriplet(values...),
            dd_state=(
                value=dd_state[1],
                deriv=dd_state[2],
                second=dd_state[3],
            ),
            condition=Float64(condition),
            relative_tail=Float64(relative_tail),
        )
    end
    return sum_at(nmax), sum_at(check_nmax)
end

# Double-double converter coefficients (A0, A1, B0, B1), as in _p_converter.
function _dd_transform(params, r)
    kappa = params.kappa
    epsilon = 2 * params.omega
    tau = (epsilon - params.m * params.a) / kappa
    transformation = _teukolsky_transformation_module()
    Tfun, _, TpTfun = getfield(transformation, :Tx)(
        params.s, epsilon, kappa, tau)
    transformed_x = getfield(transformation, :r_to_x)(r, kappa)
    Tval = DDComplex(Tfun(transformed_x))
    Xi = DDComplex(TpTfun(transformed_x))
    half_kappa_inv = DDComplex(inv(2 * kappa))
    matrix = _teukolsky_from_gsn_matrix(params)(r)
    M11, M12 = DDComplex(matrix[1, 1]), DDComplex(matrix[1, 2])
    M21, M22 = DDComplex(matrix[2, 1]), DDComplex(matrix[2, 2])
    detM = dc_sub(
        dc_mul(M11, M22),
        dc_mul(M12, M21),
    )
    M12_half = dc_mul(M12, half_kappa_inv)
    M11_half = dc_mul(M11, half_kappa_inv)
    A0 = dc_div(
        dc_mul(Tval, dc_add(
            M22,
            dc_mul(M12_half, Xi),
        )),
        detM,
    )
    A1 = dc_div(dc_mul(Tval, M12_half), detM)
    B0 = dc_div(
        dc_mul(Tval, dc_sub(
            dc_neg(M21),
            dc_mul(M11_half, Xi),
        )),
        detM,
    )
    B1 = dc_neg(dc_div(
        dc_mul(Tval, M11_half),
        detM,
    ))
    return A0, A1, B0, B1
end

# GSN (X, dX/dx) times `scale` from a double-double P state, with the condition of the
# linear combination.
function _dd_convert(
    params,
    state,
    direct_x,
    r,
    scale,
)
    A0, A1, B0, B1 = _dd_transform(params, r)
    x_term_0 = dc_mul(A0, state.value)
    x_term_1 = dc_mul(A1, state.deriv)
    d_term_0 = dc_mul(B0, state.value)
    d_term_1 = dc_mul(B1, state.deriv)
    X = dc_add(x_term_0, x_term_1)
    dXdrstar = dc_add(d_term_0, d_term_1)
    dxdrstar = _direct_dx_drstar(params, direct_x)
    iszero(dxdrstar) && error(
        "direct GSN MST cannot convert dX/drstar at dx/drstar = 0.")
    applied_scale = DDComplex(scale)
    X = dc_mul(applied_scale, X)
    dXdx = dc_mul(
        applied_scale,
        dc_div(dXdrstar, DDComplex(dxdrstar)),
    )
    X_value = dc_value(X)
    dXdx_value = dc_value(dXdx)
    all(_finite_complex, (X_value, dXdx_value)) ||
        throw(MSTCertificateError(
            :dd_regularized_type1,
            "nonfinite regularized Type-1 GSN state.",
        ))
    combination_condition = max(
        (
            abs(dc_value(x_term_0)) +
            abs(dc_value(x_term_1))
        ) / max(abs(dc_value(X)), floatmin(Float64)),
        (
            abs(dc_value(d_term_0)) +
            abs(dc_value(d_term_1))
        ) / max(abs(dc_value(dXdrstar)), floatmin(Float64)),
    )
    return (
        X=ComplexF64(X_value),
        dXdx=ComplexF64(dXdx_value),
        condition=Float64(combination_condition),
    )
end

# Regularized Type-1 version of _dd_pin_state.
function _dd_reg_pin(
    coefficients::DirectCoefficientSet,
    data::MSTSeriesData,
    offset::DDComplex,
    x,
    r,
    mst_coordinate;
    scale,
)
    full, check = _dd_reg_pair(data, offset, mst_coordinate)
    residual = _p_residual(
        data.params, mst_coordinate, full.state)
    check_residual = _p_residual(
        data.params, mst_coordinate, check.state)
    max(residual, check_residual) <= MST_RESIDUAL_TOL ||
        throw(MSTCertificateError(
            :dd_physical_in_residual,
            "regularized DD physical-in residual rejected: " *
            "full=$residual, check=$check_residual.",
        ))
    full_gsn = _dd_convert(
        coefficients.params,
        full.dd_state,
        x,
        r,
        scale,
    )
    check_gsn = _dd_convert(
        coefficients.params,
        check.dd_state,
        x,
        r,
        scale,
    )
    truncation = Float64(_amp_distance(
        (full_gsn.X, full_gsn.dXdx),
        (check_gsn.X, check_gsn.dXdx),
    ))
    truncation <= MST_AMP_TRUNCATION_MAX ||
        throw(MSTCertificateError(
            :dd_physical_in_truncation,
            "regularized DD physical-in truncation rejected: " *
            "agreement=$truncation.",
        ))
    condition = max(
        full.condition,
        check.condition,
        full_gsn.condition,
        check_gsn.condition,
    )
    condition <= MST_DD_CONDITION_MAX ||
        throw(MSTCertificateError(
            :dd_physical_in_condition,
            "regularized DD physical-in condition rejected: " *
            "condition=$condition.",
        ))
    roundoff = (2MST_DD_REG_NMAX + 1) *
        eps(Float64)^2 * condition
    return (
        X=full_gsn.X,
        dXdx=full_gsn.dXdx,
        estimated_relerr=Float64(max(
            residual,
            check_residual,
            truncation,
            full.relative_tail,
            check.relative_tail,
            roundoff,
        )),
        residual=Float64(residual),
        truncation,
        condition=Float64(condition),
        representation=:regularized_type1,
    )
end

# Double-double physical-in GSN state at direct_x from the Type-1 series (regularized
# near a pole of Gamma(c)), with residual, truncation and condition checks.
function _dd_pin_state(
    coefficients::DirectCoefficientSet,
    data::MSTSeriesData,
    offset::DDComplex,
    direct_x::Real;
    scale,
)
    x = Float64(direct_x)
    0.0 < x < 1.0 || throw(DomainError(x,
        "DD physical-in evaluation requires x in (0, 1)."))
    r = _direct_x_to_r(coefficients.params, x)
    mst_coordinate = mst_x(data.params, r)
    _dd_regime(data.params) && return _dd_reg_pin(
        coefficients,
        data,
        offset,
        x,
        r,
        mst_coordinate;
        scale,
    )
    full, check = _dd_type1_pair(
        data, offset, mst_coordinate)
    residual = _p_residual(
        data.params, mst_coordinate, full.state)
    check_residual = _p_residual(
        data.params, mst_coordinate, check.state)
    max(residual, check_residual) <= MST_RESIDUAL_TOL ||
        throw(MSTCertificateError(
            :dd_physical_in_residual,
            "DD physical-in residual rejected: " *
            "full=$residual, check=$check_residual.",
        ))
    transform = _p_converter(coefficients.params)(r)
    X, dXdx = _p_to_gsn_dx(
        coefficients,
        full.state.value,
        full.state.deriv,
        x,
        transform,
    )
    check_X, check_dXdx = _p_to_gsn_dx(
        coefficients,
        check.state.value,
        check.state.deriv,
        x,
        transform,
    )
    applied_scale = ComplexF64(scale)
    X *= applied_scale
    dXdx *= applied_scale
    check_X *= applied_scale
    check_dXdx *= applied_scale
    truncation = Float64(_amp_distance(
        (X, dXdx), (check_X, check_dXdx)))
    condition = max(full.condition, check.condition)
    truncation <= MST_AMP_TRUNCATION_MAX ||
        throw(MSTCertificateError(
            :dd_physical_in_truncation,
            "DD physical-in truncation rejected: " *
            "agreement=$truncation.",
        ))
    condition <= MST_DD_CONDITION_MAX ||
        throw(MSTCertificateError(
            :dd_physical_in_condition,
            "DD physical-in condition rejected: condition=$condition.",
        ))
    roundoff = (2MST_RESIDUAL_NMAX + 1) *
        eps(Float64)^2 * condition
    return (
        X=ComplexF64(X),
        dXdx=ComplexF64(dXdx),
        estimated_relerr=Float64(max(
            residual,
            check_residual,
            truncation,
            roundoff,
        )),
        residual=Float64(residual),
        truncation,
        condition=Float64(condition),
        representation=:native_type1,
    )
end

# Physical-in P and dP/dx at r from the budgeted Type-1 series normalized by `norm`,
# with a radial-equation residual check.
function _mst_pin_result(data::MSTSeriesData, r, norm::LogNormSum)
    x = mst_x(data.params, r)
    budget = _mst_budget(
        n -> _type1_triplet(data, n, x, norm),
        n -> _type1_error(data, n);
        reltol=MST_REL_ERROR,
        budget_tol=MST_RESIDUAL_TOL,
        nmax=MST_RESIDUAL_NMAX,
    )
    budget.status == :OK || throw(MSTCertificateError(
        :physical_in_sum,
        "MST physical-in sum rejected: status=$(budget.status), " *
        "error=$(budget.estimated_relerr), condition=$(budget.condition)",
    ))
    residual = _p_residual(data.params, x, budget.sum)
    residual <= MST_RESIDUAL_TOL || throw(MSTCertificateError(
        :physical_in_residual,
        "MST physical-in residual rejected: residual=$(residual), " *
        "error=$(budget.estimated_relerr)",
    ))
    return (
        P=budget.sum.value,
        Px=budget.sum.deriv,
        estimated_relerr=Float64(max(
            budget.estimated_relerr, residual)),
        residual=Float64(residual),
        tail=Float64(budget.estimated_relerr),
        condition=Float64(budget.condition),
        nmin=Int(budget.nmin),
        nmax=Int(budget.nmax),
        budget_status=budget.status,
    )
end

function _mst_pin_p(data::MSTSeriesData, r, norm::LogNormSum)
    result = _mst_pin_result(data, r, norm)
    return result.P, result.Px, result.estimated_relerr, result.residual
end

function _mst_state_at(coefficients::DirectCoefficientSet, branch::Symbol, direct_x::Float64)
    params = coefficients.params
    nu_offset = _mst_nu_offset(
        params.s, params.l, params.m, params.a, params.omega, params.lambda)
    nu_value = ComplexF64(params.l) + nu_offset
    data = MSTSeriesData(MSTParams(params.s, params.l, params.m, params.a, params.omega,
        params.lambda, nu_value, nu_offset))
    return _mst_state_at(coefficients, data, branch, direct_x)
end

# MSTSeriesData at the solved nu of the mode.
function _mst_data(params::DirectGSNParameters)
    nu_offset = _mst_nu_offset(
        params.s, params.l, params.m, params.a, params.omega, params.lambda)
    nu_value = ComplexF64(params.l) + nu_offset
    return MSTSeriesData(MSTParams(params.s, params.l, params.m, params.a, params.omega,
        params.lambda, nu_value, nu_offset))
end

_mst_data(coefficients::DirectCoefficientSet) =
    _mst_data(coefficients.params)

function _mst_state_at(
    coefficients::DirectCoefficientSet,
    data::MSTSeriesData,
    branch::Symbol,
    direct_x::Number,
    transform=nothing,
    sequence=nothing,
    budget=nothing,
    factor=nothing;
    radius=nothing,
)
    return _mst_state_at_params(
        coefficients.params,
        data,
        branch,
        direct_x,
        transform,
        sequence,
        budget,
        factor;
        radius,
    )
end

# Unit-normalized GSN state (X, dX/dx, relative error) of :in or :out at direct_x; the
# optional arguments reuse a converter, sequence, budget and normalization factor.
function _mst_state_at_params(
    params::DirectGSNParameters,
    data::MSTSeriesData,
    branch::Symbol,
    direct_x::Number,
    transform=nothing,
    sequence=nothing,
    budget=nothing,
    factor=nothing;
    radius=nothing,
)
    r = radius === nothing ? _direct_x_to_r(params, direct_x) : radius
    P, Px, err, log_scale = branch == :in ?
        _mst_in_p_pair(data, r, sequence, budget) :
        _mst_out_p_pair(data, r, sequence, budget)
    matrix = transform === nothing ? _p_converter(params)(r) : transform
    X, dXdx = _p_to_gsn_dx(params, P, Px, direct_x, matrix; radius)
    factor = factor === nothing ? _raw_unit_logfactor(data, params, branch) : factor
    unit_log_scale = log_scale - factor.logabs
    phase = conj(factor.phase)
    total_error = min(Inf, err + factor.estimated_relerr)
    return _scaled_product(unit_log_scale, phase * X),
        _scaled_product(unit_log_scale, phase * dXdx), Float64(total_error)
end

const MST_RESIDUAL_TOL = 1.0e-13
const MST_RESIDUAL_NMAX = 180
const MST_BOUNDARY_ITER = 80
const MST_BISECTION_ITER = 2
const MST_BISECTION_SCORE_RATIO = 10.0
const MST_Y_SAFETY = 0.8
const MST_SEED_Y_CAP = 1.0e4
const MST_TAU_PREFLIGHT = parse(Float64,
    get(ENV, "DIRECT_GSN_MST_TAU_PREFLIGHT", "1.0e2"))
const MST_TAU_SEED_PRODUCT = 2.0

struct MSTResidual
    residual::Float64
    tail::Float64
    nmin::Int
    nmax::Int
    status::String
end

_finite_triplet(t::MSTTriplet) =
    _finite_complex(t.value) && _finite_complex(t.deriv) && _finite_complex(t.second)

# Coefficients q1, q2 of x(1 - x) P'' + q1 P' + q2 P = 0, the equation checked by
# _p_residual.
function _p_q_coefficients(p::MSTParams, x)
    s = ComplexF64(p.s)
    epsc = ComplexF64(p.epsilon)
    kappa = ComplexF64(p.kappa)
    tau = ComplexF64(p.tau)
    xc = ComplexF64(x)
    q1 = 1 - s - 1im * epsc - 1im * tau -
        2 * (1 - 1im * epsc * kappa - 1im * tau) * xc -
        2im * epsc * kappa * xc^2
    q2 = p.lambda - epsc^2 + 1im * epsc * kappa * (1 - 2s) + s * (s + 1) +
        tau * (1im + tau) +
        2 * epsc * kappa * (-1im + 1im * s + epsc - tau) * xc
    return q1, q2
end

# Relative residual of that equation for a (P, P', P'') triplet.
function _p_residual(p::MSTParams, x, triplet::MSTTriplet)
    xc = ComplexF64(x)
    q1, q2 = _p_q_coefficients(p, xc)
    p2 = xc * (1 - xc)
    terms = (p2 * triplet.second, q1 * triplet.deriv, q2 * triplet.value)
    numer = terms[1] + terms[2] + terms[3]
    denom = abs(terms[1]) + abs(terms[2]) + abs(terms[3])
    return abs(numer) / max(denom, floatmin(Float64))
end

# Common log scale of the infinity-out terms, taken from the n = 0 term.
function _out_term_log_scale(data::MSTSeriesData, x, sequence)
    p = data.params
    epsc = ComplexF64(p.epsilon)
    state = _h_state(sequence, 0)
    log_weight = log_series_coefficient!(data, 0) +
        log_pochhammer(_nu_shift(p, 1 + p.s) - I * epsc, 0) -
        log_pochhammer(_nu_shift(p, 1 - p.s) + I * epsc, 0)
    dc_dx = 2I * epsc * p.kappa
    norm = max(
        abs(state.value),
        abs(dc_dx * state.deriv),
        abs(dc_dx^2 * state.second),
        floatmin(Float64),
    )
    return Float64(real(log_weight) + state.log_scale + log(norm))
end

# n-th term of the infinity-out (infinity-in below) series and its first two
# x-derivatives, built from the U-function state h_n.
function _out_term_triplet(data::MSTSeriesData, n::Integer, x,
        sequence=nothing, common_log_scale::Float64=0.0)
    p = data.params
    epsc = ComplexF64(p.epsilon)
    mst_u = sequence === nothing ? _out_sequence(data, x) : sequence
    state = _h_state(mst_u, n)
    log_weight = log_series_coefficient!(data, n) +
        log_pochhammer(_nu_shift(p, 1 + p.s) - 1im * epsc, n) -
        log_pochhammer(_nu_shift(p, 1 - p.s) + 1im * epsc, n)
    sign = isodd(n) ? -1.0 : 1.0
    dc_dx = 2im * epsc * p.kappa
    scale_log = log_weight + state.log_scale - common_log_scale
    value = _scaled_product(scale_log, sign * state.value)
    deriv = _scaled_product(scale_log, sign * dc_dx * state.deriv)
    second = _scaled_product(scale_log, sign * dc_dx^2 * state.second)
    return MSTTriplet(ComplexF64(value), ComplexF64(deriv), ComplexF64(second))
end

function _in_term_triplet(data::MSTSeriesData, n::Integer, x, sequence=nothing)
    p = data.params
    epsc = ComplexF64(p.epsilon)
    delta = 2 * p.nu_offset
    alpha = epsc * p.kappa
    z = 2im * alpha * (1 - ComplexF64(x))
    zp = -2im * alpha
    mst_u = sequence === nothing ? _in_sequence(data, x) : sequence
    state = _h_state(mst_u, n)
    a0 = ComplexF64(p.l + 1 - p.s) + 1im * epsc
    actual_a = a0 + delta / 2
    log_z = log(z)
    scale_log = log_series_coefficient!(data, n) + a0 * log_z +
        (delta / 2) * log_z + state.log_scale
    value = _scaled_product(scale_log, state.value)
    deriv = _scaled_product(
        scale_log,
        zp * (actual_a * state.value / z + state.deriv),
    )
    second = _scaled_product(
        scale_log,
        zp^2 * (state.second + 2 * actual_a * state.deriv / z +
            actual_a * (actual_a - 1) * state.value / z^2),
    )
    return MSTTriplet(ComplexF64(value), ComplexF64(deriv), ComplexF64(second))
end

function _jh_log_derivatives(p::MSTParams, x)
    epsc = ComplexF64(p.epsilon)
    kappa = ComplexF64(p.kappa)
    tau = ComplexF64(p.tau)
    s = ComplexF64(p.s)
    xc = ComplexF64(x)
    alpha = epsc * kappa
    mu = -s - 1im * (epsc + tau) / 2
    eta = 1im * (epsc - tau) / 2
    L = 1im * alpha + mu / xc - eta / (1 - xc)
    Lp = -mu / xc^2 - eta / (1 - xc)^2
    return L, Lp
end

function _out_log_transform(p::MSTParams, x)
    epsc = ComplexF64(p.epsilon)
    kappa = ComplexF64(p.kappa)
    tau = ComplexF64(p.tau)
    s = ComplexF64(p.s)
    alpha = epsc * kappa
    z = alpha * (1 - ComplexF64(x))
    A = _nu_shift(p) + 1im * (epsc + tau) / 2
    B = -s - 1im * (epsc + tau) / 2
    Lu = -alpha * (1im + A / z + B / (z - alpha))
    Lup = alpha^2 * (-A / z^2 - B / (z - alpha)^2)
    Lt, Ltp = _jh_log_derivatives(p, x)
    return Lu - Lt, Lup - Ltp
end

# P, P' and P'' of the infinity-out (infinity-in below) solution from the series
# triplet and the log-derivative of its prefactor.
function _out_p_triplet(data::MSTSeriesData, r, series::MSTTriplet)
    p = data.params
    x = mst_x(p, r)
    Lg, Lgp = _out_log_transform(p, x)
    P = series.value
    Px = series.deriv + Lg * series.value
    Pxx = series.second + 2 * Lg * series.deriv + (Lg^2 + Lgp) * series.value
    return MSTTriplet(ComplexF64(P), ComplexF64(Px), ComplexF64(Pxx))
end

function _in_p_triplet(data::MSTSeriesData, r, series::MSTTriplet)
    p = data.params
    x = mst_x(p, r)
    exponent = -1im * p.epsilon + p.s + 1im * p.tau - 1
    one_minus_x = 1 - ComplexF64(x)
    L = -exponent / one_minus_x
    Lp = -exponent / one_minus_x^2
    P = series.value
    Px = series.deriv + L * series.value
    Pxx = series.second + 2 * L * series.deriv + (L^2 + Lp) * series.value
    return MSTTriplet(ComplexF64(P), ComplexF64(Px), ComplexF64(Pxx))
end

# Budgeted series sum at r and the radial-equation residual of its transform.
function _budget_residual(data::MSTSeriesData, termfun, transformfun, r, sequence;
        log_scale::Float64=0.0)
    p = data.params
    x = mst_x(p, r)
    budget = _mst_budget(
        n -> termfun(data, n, x),
        n -> _term_error(data, sequence, n);
        reltol=MST_REL_ERROR,
        budget_tol=MST_RESIDUAL_TOL,
        nmax=MST_RESIDUAL_NMAX,
        log_scale=log_scale,
    )
    _finite_triplet(budget.sum) ||
        return MSTResidual(Inf, Inf, budget.nmin, budget.nmax,
            "BUDGET_$(budget.status)"), budget
    scale = max(
        abs(budget.sum.value),
        abs(budget.sum.deriv),
        abs(budget.sum.second),
        floatmin(Float64),
    )
    scaled = MSTTriplet(
        budget.sum.value / scale,
        budget.sum.deriv / scale,
        budget.sum.second / scale,
    )
    P = transformfun(scaled)
    _finite_triplet(P) ||
        return MSTResidual(Inf, budget.estimated_relerr,
            budget.nmin, budget.nmax, "NONFINITE_TRANSFORM"), budget
    residual = _p_residual(p, x, P)
    status =
        budget.status != :OK ? "BUDGET_$(budget.status)" :
        residual <= MST_RESIDUAL_TOL ? "OK" : "NO_RESIDUAL_PASS"
    return MSTResidual(
        residual,
        budget.estimated_relerr,
        budget.nmin,
        budget.nmax,
        status,
    ), budget
end

@inline _residual_ok(result::MSTResidual) =
    result.status == "OK" || result.status == "NOT_REQUIRED"

@inline function _residual_score(result::MSTResidual)
    result.status == "NOT_REQUIRED" && return 0.0
    return max(result.residual, result.tail)
end

# Sum and residual check of the :in or :out series at r; an error gives an Inf residual.
function _branch_check(data::MSTSeriesData, branch::Symbol, r)
    try
        x = mst_x(data.params, r)
        if branch == :out
            sequence = _out_sequence(data, x)
            log_scale = _out_term_log_scale(data, x, sequence)
            termfun = (mst_data, n, point) ->
                _out_term_triplet(mst_data, n, point, sequence, log_scale)
            residual, budget = _budget_residual(data, termfun,
                series -> _out_p_triplet(data, r, series), r, sequence;
                log_scale=log_scale)
            return (residual=residual, sequence=sequence, budget=budget)
        elseif branch == :in
            sequence = _in_sequence(data, x)
            termfun = (mst_data, n, point) ->
                _in_term_triplet(mst_data, n, point, sequence)
            residual, budget = _budget_residual(data, termfun,
                series -> _in_p_triplet(data, r, series), r, sequence)
            return (residual=residual, sequence=sequence, budget=budget)
        end
    catch
    end
    return (
        residual=MSTResidual(Inf, Inf, 0, 0, "ERROR"),
        sequence=nothing,
        budget=nothing,
    )
end

# Residual checks of the requested branches at x = 1 - |omega| y.
function _mst_residual_at_y(
    coefficients::DirectCoefficientSet,
    data::MSTSeriesData,
    y::Float64,
    branches=(:in, :out),
)
    params = coefficients.params
    direct_x = 1.0 - abs(params.omega) * y
    0.0 < direct_x < 1.0 || error("MST seed x outside (0, 1)")
    r = _direct_x_to_r(params, direct_x)
    skipped = MSTResidual(0.0, 0.0, 0, 0, "NOT_REQUIRED")
    skipped_check = (residual=skipped, sequence=nothing, budget=nothing)
    in_check = :in in branches ? _branch_check(data, :in, r) : skipped_check
    out_check = :out in branches ? _branch_check(data, :out, r) : skipped_check
    rin = in_check.residual
    rout = out_check.residual
    score = max(_residual_score(rin), _residual_score(rout))
    finite = isfinite(score) && isfinite(rin.tail) && isfinite(rout.tail)
    valid = finite && _residual_ok(rin) && _residual_ok(rout) &&
        score <= MST_RESIDUAL_TOL
    return (
        valid=valid,
        finite=finite,
        x=direct_x,
        y=y,
        score=score,
        in_residual=rin,
        out_residual=rout,
        in_sequence=in_check.sequence,
        out_sequence=out_check.sequence,
        in_budget=in_check.budget,
        out_budget=out_check.budget,
    )
end

# Factor by which the trial seed y shrinks (independent of the score).
_mst_shrink(score) = sqrt(10.0)

# Seed point for the MST series, in y = (1 - x) / |omega|: tries a capped fraction of the
# matching y, shrinks y until all branch residuals pass, refines the boundary by
# geometric bisection, and returns the seed with its cached sequences and budgets.
function _mst_residual_seed_x(
    coefficients::DirectCoefficientSet,
    data::MSTSeriesData,
    match_x::Float64,
    branches=(:in, :out),
)
    params = coefficients.params
    omega = abs(params.omega)
    omega > 0 || return (
        seed_x=match_x,
        seed_y=0.0,
        safe_y=0.0,
        unsafe_y=NaN,
        evaluations=0,
        in_sequence=nothing,
        out_sequence=nothing,
        in_budget=nothing,
        out_budget=nothing,
        status="ZERO_FREQUENCY",
    )
    y_start = (1.0 - match_x) / omega
    y_start > 0 || return (
        seed_x=match_x,
        seed_y=0.0,
        safe_y=0.0,
        unsafe_y=NaN,
        evaluations=0,
        in_sequence=nothing,
        out_sequence=nothing,
        in_budget=nothing,
        out_budget=nothing,
        status="INVALID_MATCH",
    )

    candidate_y = MST_Y_SAFETY * min(y_start, MST_SEED_Y_CAP)
    tau = abs(data.params.tau)
    if tau > MST_TAU_PREFLIGHT
        candidate_y = min(candidate_y,
            MST_TAU_SEED_PRODUCT / (omega * tau))
    end
    candidate_y = max(100eps(Float64), candidate_y)
    candidate = try
        _mst_residual_at_y(coefficients, data, candidate_y, branches)
    catch
        nothing
    end
    if candidate !== nothing && candidate.valid
        seed_x = 1.0 - omega * candidate_y
        return (
            seed_x=0.0 < seed_x < 1.0 ? seed_x : match_x,
            seed_y=candidate_y,
            safe_y=candidate_y,
            unsafe_y=NaN,
            safe_score=candidate.score,
            unsafe_score=NaN,
            evaluations=1,
            in_sequence=candidate.in_sequence,
            out_sequence=candidate.out_sequence,
            in_budget=candidate.in_budget,
            out_budget=candidate.out_budget,
            status="OK_FAST_SEED",
        )
    end

    y_unsafe = candidate_y
    y_safe = NaN
    unsafe_score = candidate === nothing ? Inf : candidate.score
    safe_score = Inf
    safe = nothing
    evaluations = 1
    y = candidate_y / _mst_shrink(unsafe_score)
    for _ in 2:MST_BOUNDARY_ITER
        evaluations += 1
        ev = try
            _mst_residual_at_y(coefficients, data, y, branches)
        catch
            nothing
        end
        if ev !== nothing && ev.valid
            y_safe = y
            safe_score = ev.score
            safe = ev
            break
        end
        score = ev === nothing ? Inf : ev.score
        y_unsafe = y
        unsafe_score = score
        y /= _mst_shrink(score)
        y <= 100eps(Float64) && break
    end
    safe === nothing && return (
        seed_x=match_x,
        seed_y=candidate_y,
        safe_y=NaN,
        unsafe_y=y_unsafe,
        safe_score=Inf,
        unsafe_score=unsafe_score,
        evaluations=evaluations,
        in_sequence=nothing,
        out_sequence=nothing,
        in_budget=nothing,
        out_budget=nothing,
        status="NO_SAFE_SEED",
    )
    refine_boundary = !isfinite(unsafe_score) ||
        unsafe_score <= MST_BISECTION_SCORE_RATIO * MST_RESIDUAL_TOL
    if isfinite(y_unsafe) && refine_boundary
        for _ in 1:MST_BISECTION_ITER
            mid = sqrt(y_safe * y_unsafe)
            evaluations += 1
            ev = try
                _mst_residual_at_y(coefficients, data, mid, branches)
            catch
                nothing
            end
            if ev !== nothing && ev.valid
                y_safe = mid
                safe_score = ev.score
                safe = ev
            else
                y_unsafe = mid
                unsafe_score = ev === nothing ? Inf : ev.score
            end
        end
    end
    y_seed = max(100eps(Float64), MST_Y_SAFETY * y_safe)
    evaluations += 1
    ev = try
        _mst_residual_at_y(coefficients, data, y_seed, branches)
    catch
        nothing
    end
    chosen_y = ev !== nothing && ev.valid ? y_seed : y_safe
    selected = ev !== nothing && ev.valid ? ev : safe
    seed_x = 1.0 - omega * chosen_y
    valid_seed_x = 0.0 < seed_x < 1.0 ? seed_x : match_x
    return (
        seed_x=valid_seed_x,
        seed_y=chosen_y,
        safe_y=y_safe,
        unsafe_y=y_unsafe,
        safe_score=safe_score,
        unsafe_score=unsafe_score,
        evaluations=evaluations,
        in_sequence=selected === nothing ? nothing : selected.in_sequence,
        out_sequence=selected === nothing ? nothing : selected.out_sequence,
        in_budget=selected === nothing ? nothing : selected.in_budget,
        out_budget=selected === nothing ? nothing : selected.out_budget,
        status=ev !== nothing && ev.valid ? "OK" : "SAFE_BOUNDARY",
    )
end

"""
    direct_mst_plan(coefficients, match_x; branches=(:in, :out))

Plan for MST infinity states of `branches` near `match_x`: the series data, a seed point
where the series pass their residual checks (`selection`), the P-to-GSN converter, the
physical-in normalization and the raw-to-unit factors of each branch.
"""
function direct_mst_plan(
    coefficients::DirectCoefficientSet,
    match_x::Float64;
    branches=(:in, :out),
)
    data = _mst_data(coefficients)
    selection = _mst_residual_seed_x(coefficients, data, match_x, branches)
    pin_norm = :in in branches ? _incoming_raw_factor(data) : nothing
    in_factor = :in in branches ?
        _raw_unit_logfactor(data, coefficients.params, :in, pin_norm) : nothing
    out_factor = :out in branches ?
        _raw_unit_logfactor(data, coefficients.params, :out) : nothing
    return (
        data=data,
        seed_x=selection.seed_x,
        converter=_p_converter(coefficients.params),
        pin_norm=pin_norm,
        in_factor=in_factor,
        out_factor=out_factor,
        in_sequence=selection.in_sequence,
        out_sequence=selection.out_sequence,
        in_budget=selection.in_budget,
        out_budget=selection.out_budget,
        selection=selection,
    )
end

"""
    direct_mst_eval_plan(plan; pin_scale=nothing, out_scale=nothing, logscaled=false, handoff_x=nothing)

The fields of a `direct_mst_plan` result needed for evaluation, with optional scales,
a log-scaling flag, and `handoff_x` replacing the seed point.
"""
direct_mst_eval_plan(
    plan;
    pin_scale=nothing,
    out_scale=nothing,
    logscaled=false,
    handoff_x=nothing,
) = (
    data=plan.data,
    seed_x=handoff_x === nothing ? Float64(plan.seed_x) : Float64(handoff_x),
    converter=plan.converter,
    pin_norm=plan.pin_norm,
    pin_scale=pin_scale,
    in_factor=plan.in_factor,
    out_factor=plan.out_factor,
    out_scale=out_scale,
    logscaled=Bool(logscaled),
)

@inline function _normalized_mst_branch(branch::Symbol)
    return branch in (:in, :ingoing, :IN, :down, :DOWN) ? :in :
        branch in (:out, :outgoing, :UP, :up) ? :out :
        throw(ArgumentError("MST infinity branch must be :in or :out."))
end

# Log-scaled MST state of a normalized branch at x with its relative error; `sequence`
# and `budget` reuse cached work.
function _mst_logscaled_state(
    coefficients::DirectCoefficientSet,
    plan,
    normalized::Symbol,
    x::Number;
    sequence=nothing,
    budget=nothing,
    radius=nothing,
)
    factor = normalized == :in ? plan.in_factor : plan.out_factor
    factor === nothing &&
        error("MST evaluation plan does not contain branch $(normalized).")
    r = radius === nothing ? _direct_x_to_r(coefficients.params, x) : radius
    P, Px, estimated_relerr, series_log_scale = normalized == :in ?
        _mst_in_p_pair(plan.data, r, sequence, budget) :
        _mst_out_p_pair(plan.data, r, sequence, budget)
    X, dXdx = _p_to_gsn_dx(
        coefficients, P, Px, x, plan.converter(r); radius)
    state = direct_logscaled_state(
        conj(factor.phase) * X,
        conj(factor.phase) * dXdx,
        series_log_scale - factor.logabs,
    )
    error_estimate = min(Inf,
        estimated_relerr + factor.estimated_relerr)
    return state, Float64(error_estimate)
end

"""
    direct_mst_principal_plan(coefficients; branches=(:in, :out))

Plan for MST states at the principal nu: series data, P-to-GSN converter and the
raw-to-unit factors of `branches`.
"""
function direct_mst_principal_plan(
    coefficients::DirectCoefficientSet;
    branches=(:in, :out),
)
    normalized = Tuple(_normalized_mst_branch(branch) for branch in branches)
    data = _mst_data(coefficients)
    in_factor = :in in normalized ?
        _raw_unit_logfactor(data, coefficients.params, :in) : nothing
    out_factor = :out in normalized ?
        _raw_unit_logfactor(data, coefficients.params, :out) : nothing
    return (
        data=data,
        converter=_p_converter(coefficients.params),
        in_factor=in_factor,
        out_factor=out_factor,
        branches=normalized,
    )
end

"""
    direct_mst_principal_state(coefficients, plan, branch, direct_x)

MST state of `branch` at a real or complex `direct_x` away from 0 and 1, from a
`direct_mst_principal_plan`. Returns `(X, dXdx, log_scale, estimated_relerr)`.
"""
function direct_mst_principal_state(
    coefficients::DirectCoefficientSet,
    plan,
    branch::Symbol,
    direct_x::Number,
)
    normalized = _normalized_mst_branch(branch)
    normalized in plan.branches ||
        error("principal MST plan does not contain branch $(normalized).")
    x = ComplexF64(direct_x)
    _finite_complex(x) && !iszero(x) && !iszero(1 - x) ||
        throw(DomainError(x,
            "principal MST evaluation requires a finite x away from 0 and 1."))
    state, estimated_relerr = _mst_logscaled_state(
        coefficients, plan, normalized, x)
    materialized = direct_materialize_logscaled_state(state)
    return (
        X=materialized.X,
        dXdx=materialized.dXdx,
        log_scale=state.log_scale,
        estimated_relerr,
    )
end

"""
    MSTPhysicalPlan

Plan returned by `direct_mst_physical_plan`: coefficient set, MST series data, branch
(`:IN` or `:UP`), certified amplitudes, P-to-GSN converter, normalization factors, the
double-double nu offset (or `nothing`) and the amplitude representation.
"""
struct MSTPhysicalPlan{C,D,A,F,P,S,O,R}
    coefficients::C
    data::D
    branch::Symbol
    amplitudes::A
    converter::F
    pin_norm::P
    pin_scale::S
    out_factor::O
    nu_offset_dd::R
    representation::Symbol
end

@inline function _physical_branch(branch::Symbol)
    normalized = _normalized_mst_branch(branch)
    return normalized == :in ? :IN : :UP
end

# MSTSeriesData with a_n for |n| <= MST_N_MAX from the double-double ratios at `offset`.
function _dd_state_data(params::MSTParams, offset::DDComplex)
    positive, negative = _dd_ratios(
        params,
        offset;
        cutoff=MST_NIA_TAIL,
        nmax=MST_N_MAX,
    )
    data = MSTSeriesData(params)
    for n in 1:MST_N_MAX
        value = dc_value(positive[n])
        _finite_complex(value) && !iszero(value) ||
            throw(MSTCertificateError(
                :dd_physical_state,
                "nonfinite positive DD MST state ratio at n=$n.",
            ))
        data.log_coeffs[n] = data.log_coeffs[n - 1] + log(value)
        data.cf_relerrs[n] = eps(Float64)
        data.cf_ok[n] = true
        data.coeff_errors[n] = data.coeff_errors[n - 1] + 2eps(Float64)
    end
    for n in 1:MST_N_MAX
        value = dc_value(negative[n])
        _finite_complex(value) && !iszero(value) ||
            throw(MSTCertificateError(
                :dd_physical_state,
                "nonfinite negative DD MST state ratio at n=$(-n).",
            ))
        data.log_coeffs[-n] = data.log_coeffs[-n + 1] + log(value)
        data.cf_relerrs[-n] = eps(Float64)
        data.cf_ok[-n] = true
        data.coeff_errors[-n] = data.coeff_errors[-n + 1] + 2eps(Float64)
    end
    return data
end

# Certified amplitudes of `branch` from the first representation that passes: principal
# nu, double-double refined nu, then shifted representations (mst_nia_amplitudes).
# Returns (series data, amplitudes, DD offset or nothing, representation).
function _physical_amplitude_plan(data::MSTSeriesData, branch::Symbol)
    errors = String[]
    try
        result = _principal_amplitudes(data, branch)
        return data, result, nothing, :principal
    catch error
        error isa InterruptException && rethrow()
        push!(errors, "principal=$(sprint(showerror, error))")
    end

    try
        result = _refined_amplitudes(data, branch)
        base = data.params
        offset = result.nu_offset_dd
        params = MSTParams(
            base.s,
            base.l,
            base.m,
            base.a,
            base.omega,
            base.lambda,
            result.nu,
            dc_value(offset),
        )
        selected = branch == :UP ?
            _dd_state_data(params, offset) : MSTSeriesData(params)
        return selected, result, offset, :dd_refined
    catch error
        error isa InterruptException && rethrow()
        push!(errors, "dd_refined=$(sprint(showerror, error))")
    end

    try
        base = data.params
        result = mst_nia_amplitudes(
            base.s,
            base.l,
            base.m,
            base.a,
            base.omega,
            base.lambda,
            branch;
            nu=base.nu,
        )
        centered = result.nu - round(Int, real(result.nu))
        selected = _amp_data(base, centered, result.shift)
        return selected, result, nothing, :coherent_shift
    catch error
        error isa InterruptException && rethrow()
        push!(errors, "coherent_shift=$(sprint(showerror, error))")
    end

    throw(MSTCertificateError(
        :physical_plan,
        "no certified physical MST representation: $(join(errors, "; "))",
    ))
end

# Log factor from the Teukolsky transmission amplitude to the GSN one (Btrans or Ctrans).
function _physical_transmission_factor(result, params, branch::Symbol)
    conversion = _conversion_module()
    factor = branch == :IN ?
        getfield(conversion, :Btrans)(
            params.s, params.m, params.a, params.omega, params.lambda) :
        getfield(conversion, :Ctrans)(
            params.s, params.m, params.a, params.omega, params.lambda)
    factor = ComplexF64(factor)
    _finite_complex(factor) && !iszero(factor) ||
        error("physical MST transmission conversion factor is nonfinite.")
    logfactor = ComplexF64(result.transmission_log - log(factor))
    return (
        phase=ComplexF64(cis(imag(logfactor))),
        logabs=Float64(real(logfactor)),
        estimated_relerr=Float64(result.truncation_agreement),
    )
end

"""
    direct_mst_physical_plan(coefficients, branch)

Plan for the physical IN (`branch` `:in`) or UP (`:out`) MST solution, with certified
asymptotic amplitudes and normalization factors; see `direct_mst_physical_state` and
`direct_mst_physical_amplitudes`.
"""
function direct_mst_physical_plan(
    coefficients::DirectCoefficientSet,
    branch::Symbol,
)
    physical = _physical_branch(branch)
    source = _mst_data(coefficients)
    data, amplitudes, offset, representation =
        _physical_amplitude_plan(source, physical)
    pin_norm = physical == :IN && offset === nothing ?
        (hasproperty(amplitudes, :coefficient_norm) ?
            amplitudes.coefficient_norm : _incoming_raw_factor(data)) : nothing
    pin_scale = physical == :IN ?
        _pin_unit_scale(coefficients.params) : nothing
    out_factor = physical == :UP ?
        _physical_transmission_factor(
            amplitudes, coefficients.params, physical) : nothing
    return MSTPhysicalPlan(
        coefficients,
        data,
        physical,
        amplitudes,
        _p_converter(coefficients.params),
        pin_norm,
        pin_scale,
        out_factor,
        offset,
        representation,
    )
end

"""
    direct_mst_physical_state(plan, direct_x)

GSN state `(X, dXdx)` of the plan's physical solution at `direct_x` in (0, 1): IN from
the hypergeometric (Type-1) series, in double-double when the plan holds a refined nu;
UP from the U-function series. Also returns the estimated error, residual, truncation,
condition and series range.
"""
function direct_mst_physical_state(
    plan::MSTPhysicalPlan,
    direct_x::Real,
)
    x = Float64(direct_x)
    0.0 < x < 1.0 || throw(DomainError(x,
        "physical MST evaluation requires x in (0, 1)."))
    coefficients = plan.coefficients
    params = coefficients.params
    r = _direct_x_to_r(params, x)

    if plan.branch == :IN
        if plan.nu_offset_dd !== nothing
            state = _dd_pin_state(
                coefficients,
                plan.data,
                plan.nu_offset_dd,
                x;
                scale=plan.pin_scale,
            )
            order = state.representation == :regularized_type1 ?
                MST_DD_REG_NMAX : MST_RESIDUAL_NMAX
            return merge(state, (
                tail=Float64(state.truncation),
                nmin=-order,
                nmax=order,
                budget_status=:OK,
                amplitude_representation=plan.representation,
            ))
        end
        result = _mst_pin_result(plan.data, r, plan.pin_norm)
        X, dXdx = _p_to_gsn_dx(
            coefficients,
            result.P,
            result.Px,
            x,
            plan.converter(r),
        )
        X = ComplexF64(plan.pin_scale * X)
        dXdx = ComplexF64(plan.pin_scale * dXdx)
        all(_finite_complex, (X, dXdx)) || throw(MSTCertificateError(
            :physical_in_state,
            "physical-IN MST state is nonfinite.",
        ))
        return (;
            X,
            dXdx,
            estimated_relerr=result.estimated_relerr,
            residual=result.residual,
            truncation=result.tail,
            tail=result.tail,
            condition=result.condition,
            representation=:native_type1,
            amplitude_representation=plan.representation,
            nmin=result.nmin,
            nmax=result.nmax,
            budget_status=result.budget_status,
        )
    end

    check = _branch_check(plan.data, :out, r)
    residual = check.residual
    _residual_ok(residual) &&
        residual.residual <= MST_RESIDUAL_TOL ||
        throw(MSTCertificateError(
            :physical_up_state,
            "physical-UP MST state rejected: status=$(residual.status), " *
            "residual=$(residual.residual), tail=$(residual.tail).",
        ))
    X, dXdx, estimated_relerr = _mst_state_at_params(
        params,
        plan.data,
        :out,
        x,
        plan.converter(r),
        check.sequence,
        check.budget,
        plan.out_factor,
    )
    all(_finite_complex, (X, dXdx)) || throw(MSTCertificateError(
        :physical_up_state,
        "physical-UP MST state is nonfinite.",
    ))
    return (
        X=ComplexF64(X),
        dXdx=ComplexF64(dXdx),
        estimated_relerr=Float64(max(
            estimated_relerr,
            residual.residual,
            residual.tail,
        )),
        residual=Float64(residual.residual),
        truncation=Float64(residual.tail),
        tail=Float64(residual.tail),
        condition=Float64(check.budget.condition),
        representation=:type2_u,
        amplitude_representation=plan.representation,
        nmin=Int(residual.nmin),
        nmax=Int(residual.nmax),
        budget_status=check.budget.status,
    )
end

@inline function _amplitude_log(logtransmission, ratio)
    iszero(ratio) && return ComplexF64(-Inf, 0)
    return ComplexF64(logtransmission + log(ComplexF64(ratio)))
end

"""
    direct_mst_physical_amplitudes(plan)

Asymptotic amplitudes of the plan's branch at unit transmission (GSN and Teukolsky
incidence and reflection), the raw Teukolsky logs, nu, and the certificate diagnostics.
"""
function direct_mst_physical_amplitudes(plan::MSTPhysicalPlan)
    result = plan.amplitudes
    transmission_log = ComplexF64(result.transmission_log)
    teuk = (ComplexF64(result.teuk[1]), ComplexF64(result.teuk[2]))
    gsn = (ComplexF64(result.gsn[1]), ComplexF64(result.gsn[2]))
    return (
        transmission=ComplexF64(1),
        incidence=gsn[1],
        reflection=gsn[2],
        teukolsky_transmission=ComplexF64(1),
        teukolsky_incidence=teuk[1],
        teukolsky_reflection=teuk[2],
        raw_teukolsky_logs=(
            transmission=transmission_log,
            incidence=_amplitude_log(transmission_log, teuk[1]),
            reflection=_amplitude_log(transmission_log, teuk[2]),
        ),
        nu=ComplexF64(plan.data.params.nu),
        nu_offset=ComplexF64(plan.data.params.nu_offset),
        representation=plan.representation,
        certificate=result.certificate_kind,
        truncation_agreement=Float64(result.truncation_agreement),
        representation_spread=Float64(result.representation_spread),
        nearest_agreement=Float64(result.nearest_agreement),
        max_condition=Float64(result.max_condition),
        nmax=Int(result.nmax),
        check_nmax=Int(result.check_nmax),
        precision_bits=hasproperty(result, :precision_bits) ?
            Int(result.precision_bits) : 53,
        nu_residual=hasproperty(result, :nu_residual) ?
            Float64(result.nu_residual) : NaN,
    )
end

"""
    direct_mst_monodromy_state(params, direct_x; tolerance=1e-11, nu=nothing)

Infinity-out MST state at `direct_x` in (0, 1) with nu from `mst_nu_complex` (or the
given `nu`), certified at `tolerance`. Returns `(X, dXdx, estimated_relerr, nu, condition)`.
"""
function direct_mst_monodromy_state(
    params::DirectGSNParameters,
    direct_x::Real;
    tolerance::Float64=1.0e-11,
    nu=nothing,
)
    0 < tolerance <= 1.0e-8 ||
        throw(ArgumentError(
            "MST Monodromy state tolerance must be in (0, 1e-8]."))
    x = Float64(direct_x)
    0.0 < x < 1.0 || throw(DomainError(x,
        "MST Monodromy state requires x in (0, 1)."))
    nu_value = nu === nothing ?
        mst_nu_complex(
            params.s,
            params.l,
            params.m,
            params.a,
            params.omega,
            params.lambda,
        ) : ComplexF64(nu)
    data = MSTSeriesData(MSTParams(
        params.s,
        params.l,
        params.m,
        params.a,
        params.omega,
        params.lambda,
        nu_value,
    ))
    r = _direct_x_to_r(params, x)
    mst_coordinate = mst_x(data.params, r)
    sequence = _out_sequence(data, mst_coordinate)
    log_scale = _out_term_log_scale(data, mst_coordinate, sequence)
    budget = _mst_budget(
        n -> _out_term_triplet(
            data, n, mst_coordinate, sequence, log_scale),
        n -> _term_error(data, sequence, n);
        reltol=tolerance,
        budget_tol=tolerance,
        log_scale,
    )
    budget.status == :OK || throw(MSTCertificateError(
        :monodromy_state,
        "MST Monodromy state rejected: status=$(budget.status), " *
        "error=$(budget.estimated_relerr), condition=$(budget.condition).",
    ))
    X, dXdx, estimated_relerr = _mst_state_at_params(
        params,
        data,
        :out,
        x,
        nothing,
        sequence,
        budget,
        nothing,
    )
    return (
        X=ComplexF64(X),
        dXdx=ComplexF64(dXdx),
        estimated_relerr=Float64(estimated_relerr),
        nu=ComplexF64(nu_value),
        condition=Float64(budget.condition),
    )
end

"""
    direct_mst_anchor(coefficients_or_params, branch, direct_x; nu=nothing, tolerance=nothing, fallback_tolerance=nothing)

MST state `(X, dXdx)` of `branch` at `direct_x` in (0, 1). With `tolerance` (at most 1e-8)
the series budget and the radial-equation residual must meet it; `fallback_tolerance`
also accepts a sum limited by rounding, state or tail error up to that bound
(`relaxed_sum = true`). `nu` overrides the solved nu. Throws `MSTCertificateError` when
rejected.
"""
function direct_mst_anchor(
    coefficients::DirectCoefficientSet,
    branch::Symbol,
    direct_x::Real;
    nu=nothing,
    tolerance=nothing,
    fallback_tolerance=nothing,
)
    normalized = _normalized_mst_branch(branch)
    x = Float64(direct_x)
    0.0 < x < 1.0 || throw(DomainError(x,
        "MST anchor requires x in (0, 1)."))
    data = if nu === nothing
        _mst_data(coefficients)
    else
        p = coefficients.params
        MSTSeriesData(MSTParams(
            p.s, p.l, p.m, p.a, p.omega, p.lambda, ComplexF64(nu)))
    end
    if tolerance === nothing
        fallback_tolerance === nothing || throw(ArgumentError(
            "MST anchor fallback tolerance requires a primary tolerance."))
        value, derivative, estimated_relerr =
            _mst_state_at(coefficients, data, normalized, x)
        return (
            X=ComplexF64(value),
            dXdx=ComplexF64(derivative),
            estimated_relerr=Float64(estimated_relerr),
            nu=ComplexF64(data.params.nu),
            nu_offset=ComplexF64(data.params.nu_offset),
        )
    end

    tol = Float64(tolerance)
    0.0 < tol <= 1.0e-8 || throw(ArgumentError(
        "MST anchor tolerance must be in (0, 1e-8]."))
    fallback_tol = if fallback_tolerance === nothing
        nothing
    else
        value = Float64(fallback_tolerance)
        tol <= value <= 1.0e-8 || throw(ArgumentError(
            "MST anchor fallback tolerance must be in [tolerance, 1e-8]."))
        value
    end
    r = _direct_x_to_r(coefficients.params, x)
    point = mst_x(data.params, r)
    if normalized == :out
        sequence = _out_sequence(data, point)
        log_scale = _out_term_log_scale(data, point, sequence)
        termfun = n -> _out_term_triplet(
            data, n, point, sequence, log_scale)
        transformfun = series -> _out_p_triplet(data, r, series)
    else
        sequence = _in_sequence(data, point)
        log_scale = 0.0
        termfun = n -> _in_term_triplet(data, n, point, sequence)
        transformfun = series -> _in_p_triplet(data, r, series)
    end
    budget = _mst_budget(
        termfun,
        n -> _term_error(data, sequence, n);
        reltol=MST_REL_ERROR,
        budget_tol=tol,
        nmax=MST_RESIDUAL_NMAX,
        log_scale,
    )
    relaxed_sum = if budget.status == :OK
        false
    elseif fallback_tol !== nothing &&
            budget.status in (:STATE_LIMIT, :ROUNDING_LIMIT, :TAIL_LIMIT) &&
            isfinite(budget.estimated_relerr) &&
            budget.estimated_relerr <= fallback_tol
        true
    else
        throw(MSTCertificateError(
            :anchor_sum,
            "MST anchor sum rejected: status=$(budget.status), " *
            "error=$(budget.estimated_relerr), condition=$(budget.condition).",
        ))
    end
    scale = max(
        abs(budget.sum.value),
        abs(budget.sum.deriv),
        abs(budget.sum.second),
        floatmin(Float64),
    )
    normalized_sum = MSTTriplet(
        budget.sum.value / scale,
        budget.sum.deriv / scale,
        budget.sum.second / scale,
    )
    residual = _p_residual(data.params, point, transformfun(normalized_sum))
    certificate_tol = relaxed_sum ? fallback_tol : tol
    residual <= certificate_tol || throw(MSTCertificateError(
        :anchor_residual,
        "MST anchor residual rejected: residual=$(residual), " *
        "error=$(budget.estimated_relerr).",
    ))
    factor = relaxed_sum && normalized == :out ?
        _raw_unit_logfactor(
            data,
            coefficients.params,
            :out,
            _outgoing_raw_fast(data, certificate_tol),
        ) : nothing
    value, derivative, estimated_relerr = _mst_state_at(
        coefficients,
        data,
        normalized,
        x,
        nothing,
        sequence,
        budget,
        factor,
    )
    return (
        X=ComplexF64(value),
        dXdx=ComplexF64(derivative),
        estimated_relerr=Float64(max(estimated_relerr, residual)),
        nu=ComplexF64(data.params.nu),
        nu_offset=ComplexF64(data.params.nu_offset),
        condition=Float64(budget.condition),
        residual=Float64(residual),
        relaxed_sum,
        budget_status=budget.status,
    )
end

# Same as above for DirectGSNParameters.
function direct_mst_anchor(
    params::DirectGSNParameters,
    branch::Symbol,
    direct_x::Real;
    nu=nothing,
    tolerance=nothing,
    fallback_tolerance=nothing,
)
    normalized = _normalized_mst_branch(branch)
    x = Float64(direct_x)
    0.0 < x < 1.0 || throw(DomainError(x,
        "MST anchor requires x in (0, 1)."))
    data = if nu === nothing
        _mst_data(params)
    else
        MSTSeriesData(MSTParams(
            params.s,
            params.l,
            params.m,
            params.a,
            params.omega,
            params.lambda,
            ComplexF64(nu),
        ))
    end
    if tolerance === nothing
        fallback_tolerance === nothing || throw(ArgumentError(
            "MST anchor fallback tolerance requires a primary tolerance."))
        value, derivative, estimated_relerr =
            _mst_state_at_params(params, data, normalized, x)
        return (
            X=ComplexF64(value),
            dXdx=ComplexF64(derivative),
            estimated_relerr=Float64(estimated_relerr),
            nu=ComplexF64(data.params.nu),
            nu_offset=ComplexF64(data.params.nu_offset),
        )
    end

    tol = Float64(tolerance)
    0.0 < tol <= 1.0e-8 || throw(ArgumentError(
        "MST anchor tolerance must be in (0, 1e-8]."))
    fallback_tol = if fallback_tolerance === nothing
        nothing
    else
        value = Float64(fallback_tolerance)
        tol <= value <= 1.0e-8 || throw(ArgumentError(
            "MST anchor fallback tolerance must be in [tolerance, 1e-8]."))
        value
    end
    r = _direct_x_to_r(params, x)
    point = mst_x(data.params, r)
    if normalized == :out
        sequence = _out_sequence(data, point)
        log_scale = _out_term_log_scale(data, point, sequence)
        termfun = n -> _out_term_triplet(
            data, n, point, sequence, log_scale)
        transformfun = series -> _out_p_triplet(data, r, series)
    else
        sequence = _in_sequence(data, point)
        log_scale = 0.0
        termfun = n -> _in_term_triplet(data, n, point, sequence)
        transformfun = series -> _in_p_triplet(data, r, series)
    end
    budget = _mst_budget(
        termfun,
        n -> _term_error(data, sequence, n);
        reltol=MST_REL_ERROR,
        budget_tol=tol,
        nmax=MST_RESIDUAL_NMAX,
        log_scale,
    )
    relaxed_sum = if budget.status == :OK
        false
    elseif fallback_tol !== nothing &&
            budget.status in (:STATE_LIMIT, :ROUNDING_LIMIT, :TAIL_LIMIT) &&
            isfinite(budget.estimated_relerr) &&
            budget.estimated_relerr <= fallback_tol
        true
    else
        throw(MSTCertificateError(
            :anchor_sum,
            "MST anchor sum rejected: status=$(budget.status), " *
            "error=$(budget.estimated_relerr), condition=$(budget.condition).",
        ))
    end
    scale = max(
        abs(budget.sum.value),
        abs(budget.sum.deriv),
        abs(budget.sum.second),
        floatmin(Float64),
    )
    normalized_sum = MSTTriplet(
        budget.sum.value / scale,
        budget.sum.deriv / scale,
        budget.sum.second / scale,
    )
    residual = _p_residual(data.params, point, transformfun(normalized_sum))
    certificate_tol = relaxed_sum ? fallback_tol : tol
    residual <= certificate_tol || throw(MSTCertificateError(
        :anchor_residual,
        "MST anchor residual rejected: residual=$(residual), " *
        "error=$(budget.estimated_relerr).",
    ))
    factor = relaxed_sum && normalized == :out ?
        _raw_unit_logfactor(
            data,
            params,
            :out,
            _outgoing_raw_fast(data, certificate_tol),
        ) : nothing
    value, derivative, estimated_relerr = _mst_state_at_params(
        params,
        data,
        normalized,
        x,
        nothing,
        sequence,
        budget,
        factor,
    )
    return (
        X=ComplexF64(value),
        dXdx=ComplexF64(derivative),
        estimated_relerr=Float64(max(estimated_relerr, residual)),
        nu=ComplexF64(data.params.nu),
        nu_offset=ComplexF64(data.params.nu_offset),
        condition=Float64(budget.condition),
        residual=Float64(residual),
        relaxed_sum,
        budget_status=budget.status,
    )
end

"""
    direct_mst_logscaled_seed(coefficients, plan, branch)

Log-scaled MST state of `branch` at the plan's seed point, reusing the plan's cached
sequence and budget; returns `(state, estimated_relerr)`.
"""
function direct_mst_logscaled_seed(
    coefficients::DirectCoefficientSet,
    plan,
    branch::Symbol,
)
    normalized = _normalized_mst_branch(branch)
    sequence = normalized == :in ? plan.in_sequence : plan.out_sequence
    budget = normalized == :in ? plan.in_budget : plan.out_budget
    return _mst_logscaled_state(
        coefficients,
        plan,
        normalized,
        Float64(plan.seed_x);
        sequence=sequence,
        budget=budget,
    )
end

"""
    direct_mst_logscaled_state(coefficients, plan, branch, direct_x; radius=nothing)

MST state of `branch` at `direct_x` in (0, 1) as a `DirectLogScaledState`, with its
estimated relative error.
"""
function direct_mst_logscaled_state(
    coefficients::DirectCoefficientSet,
    plan,
    branch::Symbol,
    direct_x::Real;
    radius=nothing,
)
    normalized = _normalized_mst_branch(branch)
    x = Float64(direct_x)
    0.0 < x < 1.0 || throw(DomainError(x,
        "MST log-scaled evaluation requires x in (0, 1)."))
    return _mst_logscaled_state(coefficients, plan, normalized, x; radius)
end

"""
    direct_mst_scaled_state(coefficients, plan, branch, direct_x; scale, radius=nothing)

`direct_mst_logscaled_state` materialized with `scale`; returns
`(X, dXdx, estimated_relerr)`.
"""
function direct_mst_scaled_state(
    coefficients::DirectCoefficientSet,
    plan,
    branch::Symbol,
    direct_x::Real;
    scale,
    radius=nothing,
)
    state, estimated_relerr = direct_mst_logscaled_state(
        coefficients, plan, branch, direct_x; radius)
    materialized = direct_materialize_logscaled_state(state; scale=scale)
    return (
        X=materialized.X,
        dXdx=materialized.dXdx,
        estimated_relerr=estimated_relerr,
    )
end

"""
    direct_mst_pin_state(coefficients, plan, direct_x; scale=nothing, radius=nothing)

Physical-in MST state at `direct_x` in (0, 1) from the hypergeometric series with the
plan's normalization, times `scale` (default the plan's `pin_scale`, else 1). Returns
`(X, dXdx, estimated_relerr, residual)`.
"""
function direct_mst_pin_state(
    coefficients::DirectCoefficientSet,
    plan,
    direct_x::Real;
    scale=nothing,
    radius=nothing,
)
    x = Float64(direct_x)
    0.0 < x < 1.0 || throw(DomainError(x,
        "MST physical-in evaluation requires x in (0, 1)."))
    plan.pin_norm === nothing &&
        error("MST evaluation plan does not contain the physical-in normalization.")
    r = radius === nothing ? _direct_x_to_r(coefficients.params, x) : Float64(radius)
    P, Px, estimated_relerr, residual = _mst_pin_p(
        plan.data, r, plan.pin_norm)
    X, dXdx = _p_to_gsn_dx(
        coefficients, P, Px, x, plan.converter(r); radius)
    applied_scale = scale === nothing && hasproperty(plan, :pin_scale) ?
        plan.pin_scale : scale
    applied_scale === nothing && (applied_scale = ComplexF64(1))
    applied_scale = ComplexF64(applied_scale)
    return (
        X=applied_scale * X,
        dXdx=applied_scale * dXdx,
        estimated_relerr=estimated_relerr,
        residual=residual,
    )
end

"""
    direct_mst_state(coefficients, plan, branch, direct_x; radius=nothing)

Unit-normalized MST state of `branch` (`:in` or `:out`) at `direct_x` in (0, 1) from a
`direct_mst_plan`; returns `(X, dXdx, estimated_relerr)`.
"""
function direct_mst_state(
    coefficients::DirectCoefficientSet,
    plan,
    branch::Symbol,
    direct_x::Real;
    radius=nothing,
)
    normalized =
        branch in (:in, :ingoing, :IN, :down, :DOWN) ? :in :
        branch in (:out, :outgoing, :UP, :up) ? :out :
        throw(ArgumentError("MST infinity branch must be :in or :out."))
    x = Float64(direct_x)
    0.0 < x < 1.0 || throw(DomainError(x, "MST evaluation requires x in (0, 1)."))
    factor = normalized == :in ? plan.in_factor : plan.out_factor
    factor === nothing && error("MST evaluation plan does not contain branch $(normalized).")
    r = radius === nothing ? _direct_x_to_r(coefficients.params, x) : radius
    transform = plan.converter(r)
    value, derivative, error_estimate = _mst_state_at(
        coefficients,
        plan.data,
        normalized,
        x,
        transform,
        nothing,
        nothing,
        factor;
        radius,
    )
    return (X=value, dXdx=derivative, estimated_relerr=error_estimate)
end

# Low-frequency expansion of nu - l to order epsilon^2
# (Sasaki & Tagoshi, Living Rev. Relativ. 6, 6 (2003)).
function _mst_guess_offset(s, l, epsilon)
    lf = Float64(l)
    sf = Float64(s)
    correction = -2 - sf^2 / (lf * (lf + 1)) +
        ((lf + 1)^2 - sf^2)^2 / ((2 * lf + 1) * (2 * lf + 2) * (2 * lf + 3)) -
        (lf^2 - sf^2)^2 / ((2 * lf - 1) * 2 * lf * (2 * lf + 1))
    return correction * epsilon^2 / (2 * lf + 1)
end

"""
    mst_guess(s, l, epsilon)

Low-frequency estimate of the MST nu, l + O(epsilon^2).
"""
mst_guess(s, l, epsilon) = Float64(l) + _mst_guess_offset(s, l, epsilon)

# MST equation beta_0 + alpha_0 R_1 + gamma_0 L_(-1) at nu = l + nu_offset; nu is a root.
function _mst_eq_value(nu_offset, s, l, m, a, omega, lambda)
    nu_value = ComplexF64(l) + ComplexF64(nu_offset)
    p = MSTParams(s, l, m, a, omega, lambda, nu_value, nu_offset)
    r1 = rn_cf(1, p)
    lm1 = ln_cf(-1, p)
    value = beta_mst(0, p) + alpha_mst(0, p) * r1.value + gamma_mst(0, p) * lm1.value
    (!isfinite(real(value)) || !isfinite(imag(value))) &&
        error("nonfinite MST nu equation at nu-l=$nu_offset")
    return ComplexF64(value)
end

_mst_eq_offset(nu_offset, s, l, m, a, omega, lambda) =
    real(_mst_eq_value(nu_offset, s, l, m, a, omega, lambda))

mst_eq(nu_real, s, l, m, a, omega, lambda) =
    _mst_eq_offset(ComplexF64(nu_real) - ComplexF64(l), s, l, m, a, omega, lambda)

# nu - l for real omega: the guess when negligible; mst_nu_complex when the monodromy
# trace leaves [-1, 1]; otherwise the real root, bracketed around the guess and refined
# by regula falsi with bisection safeguards (secant iteration without a bracket).
function _mst_nu_offset(
    s,
    l,
    m,
    a,
    omega,
    lambda;
    maxiter::Int=60,
    tol::Float64=1.0e-13,
)
    epsilon = 2.0 * Float64(omega)
    iszero(epsilon) && return ComplexF64(0)
    center = _mst_guess_offset(s, l, epsilon)
    if !iszero(center) &&
            abs(center) <= 64eps(Float64) * max(1.0, abs(l))
        return ComplexF64(center)
    end
    # For real omega, nu is real only while cos(2 pi nu) lies in [-1, 1];
    # outside, nu = 1/2 + i nu_i or i nu_i (mod integers) and the root is
    # found in the complex plane from the monodromy.
    trace = _monodromy_trace_data(s, m, a, ComplexF64(omega), lambda)
    # An unresolved excess at +/-1 cannot establish a complex-nu branch.
    abs(real(trace.value)) - 1 <= trace.resolution ||
        return mst_nu_complex(s, l, m, a, omega, lambda) - l

    fcenter = _mst_eq_offset(center, s, l, m, a, omega, lambda)
    iszero(fcenter) && return ComplexF64(center)

    step = max(0.02 * max(abs(center), epsilon^2), eps(Float64) * max(1.0, abs(l)))
    left = center - step
    right = center + step
    fleft = _mst_eq_offset(left, s, l, m, a, omega, lambda)
    fright = _mst_eq_offset(right, s, l, m, a, omega, lambda)

    for _ in 1:24
        (iszero(fleft) || iszero(fright) || signbit(fleft) != signbit(fright)) && break
        step *= 2
        left = center - step
        right = center + step
        fleft = _mst_eq_offset(left, s, l, m, a, omega, lambda)
        fright = _mst_eq_offset(right, s, l, m, a, omega, lambda)
    end

    if iszero(fleft)
        return ComplexF64(left)
    elseif iszero(fright)
        return ComplexF64(right)
    elseif signbit(fleft) != signbit(fright)
        best_x, best_f = abs(fleft) <= abs(fright) ? (left, fleft) : (right, fright)
        previous = NaN
        for _ in 1:min(maxiter, 12)
            denom = fright - fleft
            trial = iszero(denom) ? left + (right - left) / 2 :
                (left * fright - right * fleft) / denom
            if !(left < trial < right) || trial == previous
                trial = left + (right - left) / 2
            end
            (trial == left || trial == right) && return ComplexF64(best_x)
            ftrial = _mst_eq_offset(trial, s, l, m, a, omega, lambda)
            if abs(ftrial) < abs(best_f)
                best_x, best_f = trial, ftrial
            end
            iszero(ftrial) && return ComplexF64(trial)
            if signbit(ftrial) == signbit(fleft)
                left = trial
                fleft = ftrial
            else
                right = trial
                fright = ftrial
            end
            width_tol = max(tol * max(abs(trial), floatmin(Float64)), 8eps(trial))
            abs(right - left) <= width_tol && return ComplexF64(best_x)
            previous = trial
        end
        return ComplexF64(best_x)
    end

    x0 = left
    x1 = right
    f0 = fleft
    f1 = fright
    for _ in 1:maxiter
        denom = f1 - f0
        !iszero(denom) || break
        x2 = x1 - f1 * (x1 - x0) / denom
        isfinite(x2) || break
        f2 = _mst_eq_offset(x2, s, l, m, a, omega, lambda)
        if x2 == x1 || abs(x2 - x1) <= max(tol * max(abs(x2), floatmin(Float64)), 8eps(x2))
            return ComplexF64(x2)
        end
        x0, f0 = x1, f1
        x1, f1 = x2, f2
    end
    return ComplexF64(x1)
end

# nu - l for complex omega by damped Newton iteration from the low-frequency guess.
function _mst_nu_offset(
    s,
    l,
    m,
    a,
    omega::Complex,
    lambda;
    maxiter::Int=60,
    tol::Float64=1.0e-13,
)
    epsilon = 2.0 * ComplexF64(omega)
    iszero(epsilon) && return ComplexF64(0)
    iszero(imag(epsilon)) && return _mst_nu_offset(
        s, l, m, a, real(omega), lambda; maxiter, tol)

    current = ComplexF64(_mst_guess_offset(s, l, epsilon))
    value = _mst_eq_value(current, s, l, m, a, omega, lambda)
    scale = max(abs(ComplexF64(lambda)), abs(l * (l + 1)), 1.0)
    best = current
    best_residual = abs(value) / scale
    best_residual <= tol && return best

    for _ in 1:maxiter
        h = max(1.0e-7, sqrt(eps(Float64)) * max(1.0, abs(current)))
        plus = _mst_eq_value(current + h, s, l, m, a, omega, lambda)
        minus = _mst_eq_value(current - h, s, l, m, a, omega, lambda)
        derivative = (plus - minus) / (2h)
        _finite_complex(derivative) && !iszero(derivative) || break

        step = value / derivative
        accepted = false
        trial = current
        trial_value = value
        for damping in 0:10
            factor = exp2(-damping)
            candidate = current - factor * step
            candidate_value = _mst_eq_value(
                candidate, s, l, m, a, omega, lambda)
            if abs(candidate_value) < abs(value)
                trial = candidate
                trial_value = candidate_value
                accepted = true
                break
            end
        end
        accepted || break
        current = trial
        value = trial_value
        residual = abs(value) / scale
        if residual < best_residual
            best = current
            best_residual = residual
        end
        residual <= tol && return current
        abs(step) <= tol * max(1.0, abs(current)) && break
    end

    best_residual <= 32tol || throw(MSTCertificateError(
        :complex_nu,
        "complex MST nu solve rejected: residual=$(best_residual), " *
        "nu_offset=$(best), omega=$(omega)",
    ))
    return best
end

struct MSTMonodromyParams
    a::ComplexF64
    gamma::ComplexF64
    delta::ComplexF64
    epsilon::ComplexF64
    q::ComplexF64
end

# Confluent Heun parameters of the Teukolsky radial equation for the monodromy of nu.
function _monodromy_params(s, m, a, omega, lambda)
    epsilon = 2omega
    kappa = sqrt(1 - a^2)
    tau = (epsilon - m * a) / kappa
    ach = 2im * epsilon * kappa *
        (1 - s + im * epsilon - im * tau)
    gamma_ch = 1 - s - im * (epsilon + tau)
    delta_ch = 1 + s + im * (epsilon - tau)
    epsilon_ch = 2im * epsilon * kappa
    qch = -(
        -s * (1 + s) + epsilon^2 +
        im * (2s - 1) * epsilon * kappa -
        lambda - tau * (im + tau)
    )
    return MSTMonodromyParams(
        ach, gamma_ch, delta_ch, epsilon_ch, qch)
end

mutable struct MSTMonodromySeries
    coefficients::Vector{ComplexF64}
    previous::ComplexF64
    current::ComplexF64
end

MSTMonodromySeries() = MSTMonodromySeries(ComplexF64[1], 0.0im, 1.0+0.0im)

# Extends in place the kind-1 or kind-2 series coefficients used by _monodromy_stokes
# up to order nmax.
function _monodromy_series!(
    series::MSTMonodromySeries,
    kind::Int,
    params::MSTMonodromyParams,
    nmax::Int,
)
    coefficients = series.coefficients
    ach = params.a
    gamma_ch = params.gamma
    delta_ch = params.delta
    epsilon_ch = params.epsilon
    qch = params.q
    m1 = ach / epsilon_ch - (gamma_ch + delta_ch)
    m2 = -(ach / epsilon_ch)
    previous = series.previous
    current = series.current
    for n in length(coefficients):nmax
        ratio = previous / current
        previous = 1.0 + 0.0im
        cn = ComplexF64(n)
        if kind == 1
            pref2 = (
                (ach - epsilon_ch *
                    (cn - 2 + gamma_ch + delta_ch)) *
                (ach - epsilon_ch * (cn + delta_ch - 1))
            ) / (epsilon_ch * cn)
            pref1 = 1 - delta_ch + 2ach / epsilon_ch +
                epsilon_ch - gamma_ch - cn +
                (
                    -(ach / epsilon_ch)^2 +
                    (-1 + delta_ch) * epsilon_ch +
                    ach * (-1 + delta_ch - epsilon_ch + gamma_ch) /
                        epsilon_ch +
                    qch
                ) / cn
            weight = cn - 1 - m2 + m1
        else
            pref2 = -(
                (ach + epsilon_ch * (cn - 2)) *
                (ach + epsilon_ch * (cn - gamma_ch - 1))
            ) / (epsilon_ch * cn)
            pref1 = (
                (ach + epsilon_ch * (cn - 1)) *
                (ach + epsilon_ch *
                    (cn - delta_ch + epsilon_ch - gamma_ch)) -
                epsilon_ch^2 * qch
            ) / (epsilon_ch^2 * cn)
            weight = cn - 1 + m2 - m1
        end
        current = pref2 * ratio + pref1
        iszero(current) && throw(MSTCertificateError(
            :complex_nu_monodromy,
            "zero complex-monodromy series normalization.",
        ))
        push!(coefficients,current)
        normalization = current
        for j in 0:(n - 1)
            coefficients[j + 1] *= (weight - j) / normalization
        end
        coefficients[n + 1] = 1
    end
    series.previous = previous
    series.current = current
    return coefficients
end

_monodromy_series(kind::Int, params::MSTMonodromyParams, nmax::Int) =
    _monodromy_series!(MSTMonodromySeries(),kind,params,nmax)

# Stokes-multiplier term added to cos(2 pi mu) in _monodromy_trace_data, from the two
# series truncated at nmax.
function _monodromy_stokes(
    params::MSTMonodromyParams,
    nmax::Int,
    first::Vector{ComplexF64}=_monodromy_series(1,params,nmax),
    second::Vector{ComplexF64}=_monodromy_series(2,params,nmax),
)
    ach = params.a
    gamma_ch = params.gamma
    delta_ch = params.delta
    epsilon_ch = params.epsilon
    m1 = ach / epsilon_ch - (gamma_ch + delta_ch)
    m2 = -(ach / epsilon_ch)
    first_gamma = gamma(m1 - m2)
    second_gamma = gamma(m2 - m1)
    split = max(cld(nmax, 2), 7)
    first_sum = sum(@view first[1:(split + 1)])
    second_sum = sum(
        (-1.0)^j * second[j + 1] for j in 0:split)
    product = first_sum * second_sum
    iszero(product) && throw(MSTCertificateError(
        :complex_nu_monodromy,
        "zero complex-monodromy weighted sum.",
    ))
    return (2pi)^2 * (-1.0)^(nmax - 1) /
        (2first_gamma * second_gamma * product)
end

# cos(2 pi nu) from the monodromy of the confluent Heun equation, and the
# centred mu it is measured against.
function _monodromy_trace_data(s, m, a, omega, lambda)
    params = _monodromy_params(s, m, a, omega, lambda)
    epsilon = 2omega
    m1 = params.a / params.epsilon -
        (params.gamma + params.delta)
    m2 = -(params.a / params.epsilon)
    mu = (m2 - m1) / 2
    previous = 0.0 + 0.0im
    stokes = 0.0 + 0.0im
    drift = Inf
    first = MSTMonodromySeries()
    second = MSTMonodromySeries()
    for nmax in 50:50:500
        _monodromy_series!(first,1,params,nmax)
        _monodromy_series!(second,2,params,nmax)
        stokes = _monodromy_stokes(params, nmax, first.coefficients, second.coefficients)
        _finite_complex(stokes) || throw(MSTCertificateError(
            :complex_nu_monodromy,
            "nonfinite complex-monodromy multiplier.",
        ))
        drift = nmax == 50 ? Inf :
            abs(stokes - previous) /
            max(abs(stokes), abs(previous), floatmin(Float64))
        nmax > 50 && drift <= 1.0e-11 && break
        previous = stokes
    end
    drift <= 1.0e-8 || throw(MSTCertificateError(
        :complex_nu_monodromy,
        "complex-monodromy series rejected: drift=$drift.",
    ))
    base = cos(2pi * mu)
    resolution = eps(Float64) * (abs(base) + abs(stokes)) + abs(stokes) * drift
    return (value=base + stokes, resolution=resolution)
end

_monodromy_trace(s, m, a, omega, lambda) =
    _monodromy_trace_data(s, m, a, omega, lambda).value

# Seed for nu from acos of the monodromy trace, reduced mod 1; prefers Im nu opposite in
# sign to Im omega, then the smallest |Re nu|.
function _monodromy_seed(s, m, a, omega, lambda)
    base = acos(_monodromy_trace(s, m, a, omega, lambda)) / (2pi)
    candidates = ComplexF64[]
    for sign in (-1, 1)
        candidate = sign * base
        candidate -= round(real(candidate))
        push!(candidates, candidate)
    end
    opposite_half_plane(candidate) =
        abs(imag(omega)) <= 1.0e-14 ||
        imag(candidate) * imag(omega) <= 0
    sort!(candidates; by=candidate -> (
        opposite_half_plane(candidate) ? 0 : 1,
        abs(real(candidate)),
        -abs(imag(candidate)),
    ))
    return first(candidates)
end

# Finite-difference derivatives of `equation` at `current` for several real and imaginary
# steps; _refine_nu tries each as a Newton direction.
function _nu_derivatives(equation, current, value)
    scale = max(1.0, abs(current))
    result = ComplexF64[]
    for exponent in (18, 20, 22, 24, 26, 28)
        h = exp2(-exponent) * scale
        derivative = (
            equation(current + h) - equation(current - h)
        ) / (2h)
        _finite_complex(derivative) && !iszero(derivative) &&
            push!(result, derivative)
        derivative_i = (
            equation(current + im * h) -
            equation(current - im * h)
        ) / (2im * h)
        _finite_complex(derivative_i) && !iszero(derivative_i) &&
            push!(result, derivative_i)
    end
    h = exp2(-20) * scale * (1 + im)
    derivative = (equation(current + h) - value) / h
    _finite_complex(derivative) && !iszero(derivative) &&
        push!(result, derivative)
    return result
end

# Damped Newton refinement of `seed`; returns (nu, residual scaled by `scale`).
function _refine_nu(
    equation,
    seed,
    scale;
    maxiter::Int=8,
    tol::Float64=3.2e-13,
)
    current = ComplexF64(seed)
    value = equation(current)
    best = current
    best_residual = abs(value) / scale
    best_residual <= tol &&
        return (; nu=best, residual=best_residual)
    for _ in 1:maxiter
        best_trial = current
        best_value = value
        for derivative in _nu_derivatives(equation, current, value)
            step = value / derivative
            _finite_complex(step) || continue
            abs(step) <= 1.0e-4 * max(1.0, abs(current)) || continue
            for damping in 0:6
                candidate = current - exp2(-damping) * step
                candidate_value = equation(candidate)
                if abs(candidate_value) < abs(best_value)
                    best_trial = candidate
                    best_value = candidate_value
                end
            end
        end
        best_trial == current && break
        current = best_trial
        value = best_value
        residual = abs(value) / scale
        if residual < best_residual
            best = current
            best_residual = residual
        end
        residual <= 2.0e-14 && break
    end
    return (; nu=best, residual=best_residual)
end

"""
    mst_nu_complex(s, l, m, a, omega, lambda; tol=3.2e-13)

MST nu seeded from cos(2 pi nu) given by the monodromy of the confluent Heun equation
and refined by damped Newton steps on the MST equation, trying several integer shifts of
the seed. Throws `MSTCertificateError` when no shift reaches `tol`.
"""
function mst_nu_complex(
    s,
    l,
    m,
    a,
    omega,
    lambda;
    tol::Float64=3.2e-13,
)
    omegac = ComplexF64(omega)
    !iszero(omegac) ||
        throw(ArgumentError("MST nu requires a nonzero frequency."))
    seed = _monodromy_seed(s, m, a, omegac, lambda)
    equation = nu -> _mst_eq_value(
        ComplexF64(nu) - ComplexF64(l),
        s, l, m, a, omegac, lambda)
    scale = max(
        abs(ComplexF64(lambda)), abs(l * (l + 1)), 1.0)
    shifts = Int[-(l + 1), -l, -(l + 2)]
    append!(shifts, (0, -4, -3, -2, -1, 1, 2, 3, 4))
    unique!(shifts)
    best = (; nu=ComplexF64(seed), residual=Inf)
    for shift in shifts
        candidate_seed = seed + shift
        candidate = _refine_nu(
            equation, candidate_seed, scale; tol)
        candidate.residual < best.residual && (best = candidate)
        best.residual <= tol && return best.nu
    end
    throw(MSTCertificateError(
        :complex_nu_refinement,
        "complex MST nu refinement rejected: residual=$(best.residual).",
    ))
end

"""
    mst_nu(s, l, m, a, omega, lambda; maxiter=60, tol=1e-13)

Renormalized angular momentum nu of the MST formalism for the mode (s, l, m), spin `a`,
frequency `omega` and eigenvalue `lambda`, as `l + (nu - l)`. Real frequencies use a
bracketed real root (or `mst_nu_complex` when nu is complex); complex frequencies use
damped Newton iteration.
"""
function mst_nu(s, l, m, a, omega, lambda; maxiter::Int=60, tol::Float64=1.0e-13)
    offset = _mst_nu_offset(s, l, m, a, omega, lambda; maxiter=maxiter, tol=tol)
    return ComplexF64(l) + offset
end

# Coefficients of eta(r) in powers of 1/r up to r^-4, and its evaluation at r.
@inline function _eta_coeffs(params)
    coefficient(order) = ComplexF64(eta_coefficient(
        params.s, params.m, params.a, params.omega, params.lambda, order))
    return SVector{5,ComplexF64}(
        coefficient(0),
        coefficient(-1),
        coefficient(-2),
        coefficient(-3),
        coefficient(-4),
    )
end

@inline function _eta_value(coefficients, r)
    u = inv(r)
    return (((coefficients[5] * u + coefficients[4]) * u +
        coefficients[3]) * u + coefficients[2]) * u + coefficients[1]
end

# Closed-form Wronskian denominator of the direct basis at x for :horizon or :infinity,
# from eta(r) (Abel's identity).
function _closed_abel(coefficients::DirectCoefficientSet, endpoint::Symbol, x::Number)
    params = coefficients.params
    kappa = params.kappa
    y = 1.0 - x
    r_plus = 1.0 + kappa
    r_minus = 1.0 - kappa
    r = (r_plus - r_minus * x) / y
    d = kappa + y * y + 2.0 * kappa * kappa * x - kappa * x * x
    eta_coeffs = _eta_coeffs(params)
    eta_r = _eta_value(eta_coeffs, r)

    if endpoint == :horizon
        eta_plus = _eta_value(eta_coeffs, r_plus)
        regular = d / ((1.0 + kappa) * y * y) * eta_r / eta_plus
        eta_h = (2.0 * r_plus * params.omega - params.a * params.m) /
            (2.0 * kappa)
        scale = direct_endpoint_scale(:horizon_in, params) *
            direct_endpoint_scale(:horizon_out, params)
        threshold_log = iszero(eta_h)
        wronskian_factor = threshold_log ?
            one(ComplexF64) : 2.0im * eta_h
        denominator = wronskian_factor * scale * inv(x) * regular
        status = _finite_complex(regular) && _finite_complex(denominator) &&
            !iszero(eta_plus) && !iszero(denominator) ?
            "OK" : "NONFINITE_CLOSED_FORM"
        exponent = status == "OK" ? -log(regular) : ComplexF64(NaN, NaN)
        return (
            exponent=ComplexF64(exponent),
            order=0,
            tail=0.0,
            a0_relerr=0.0,
            status=status,
            method=:closed_eta,
            endpoint=:horizon,
            coordinate=x,
            eta_H=eta_h,
            denominator=ComplexF64(denominator),
            mode=threshold_log ?
                :abel_horizon_log_closed : :abel_horizon_closed,
        )
    elseif endpoint == :infinity
        c0 = eta_coeffs[1]
        regular = d / (2.0 * kappa * kappa * x) * eta_r / c0
        scale = direct_endpoint_scale(:infinity_in, params) *
            direct_endpoint_scale(:infinity_out, params)
        denominator = 4.0im * kappa * params.omega * scale * inv(y * y) * regular
        status = _finite_complex(regular) && _finite_complex(denominator) &&
            !iszero(c0) ? "OK" : "NONFINITE_CLOSED_FORM"
        exponent = status == "OK" ? log(regular) : ComplexF64(NaN, NaN)
        return (
            exponent=ComplexF64(exponent),
            order=0,
            tail=0.0,
            a0_relerr=0.0,
            status=status,
            method=:closed_eta,
            endpoint=:infinity,
            coordinate=y,
            denominator=ComplexF64(denominator),
            mode=:abel_infinity_closed,
        )
    end
    throw(ArgumentError("direct Abel denominator endpoint must be :infinity or :horizon"))
end

"""
    direct_abel_denominator(coefficients, endpoint, match_x)

Closed-form Wronskian denominator of the direct GSN basis at `match_x` (real in (0, 1),
or complex away from 0 and 1) for `endpoint` `:infinity` or `:horizon`. Returns a named
tuple with the log `exponent`, `denominator`, `status` and `mode`.
"""
function direct_abel_denominator(
    coefficients::DirectCoefficientSet,
    endpoint::Symbol,
    match_x::Real,
)
    x = Float64(match_x)
    0.0 < x < 1.0 || throw(ArgumentError("match_x must lie in (0,1)."))
    return _closed_abel(coefficients, endpoint, x)
end

function direct_abel_denominator(
    coefficients::DirectCoefficientSet,
    endpoint::Symbol,
    match_x::Complex,
)
    x = ComplexF64(match_x)
    _finite_complex(x) && !iszero(x) && !iszero(1 - x) ||
        throw(ArgumentError("complex match_x must be finite and avoid 0 and 1."))
    return _closed_abel(coefficients, endpoint, x)
end

end
