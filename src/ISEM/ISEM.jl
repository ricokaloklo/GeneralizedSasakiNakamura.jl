module ISEM

using ..GeneralizedSasakiNakamura: BoundaryCondition, IN, UP, OUT, DOWN
using ..GeneralizedSasakiNakamura: NormalizationConvention, UNIT_GSN_TRANS, UNIT_TEUKOLSKY_TRANS
using ..GeneralizedSasakiNakamura: Mode, GSNRadialFunction, TeukolskyRadialFunction
using ..GeneralizedSasakiNakamura: ConversionFactors, Solutions, r_from_rstar, rstar_from_r
using ..GeneralizedSasakiNakamura: InitialConditions, AsymptoticExpansionCoefficients
using StaticArrays

include("Matching.jl")
using .Matching

export Matching, _Pin, _Pup, _P
export BoundaryCondition, IN, UP, OUT, DOWN
export NormalizationConvention, UNIT_GSN_TRANS, UNIT_TEUKOLSKY_TRANS
export Mode, GSNRadialFunction, TeukolskyRadialFunction, GSN_radial, Teukolsky_radial
export YRadialFunction, Y_radial
export matching_controls

const _DEFAULT_XM = -1.0
const _DEFAULT_RHOM = 1.0
const _DEFAULT_N = 40
const _DEFAULT_TOL = 1e-14
const _GSN_HORIZON_ASYMPTOTIC_ORDER = 12
const _GSN_INFINITY_ASYMPTOTIC_ORDER = 20
const _GSN_HORIZON_DELTA_R_MAX = 0.03
const _GSN_INFINITY_PHASE_MIN = 30.0
const _GSN_ASYMPTOTIC_PATCH_RELATIVE_MISMATCH_MAX = 1e-4

@inline function _is_omega_complex(omega)
    return !iszero(imag(complex(omega)))
end

@inline function _omega_phase_ratio(omega)
    ω = complex(omega)
    denom = max(abs(real(ω)), abs(imag(ω)))
    iszero(denom) && return 0.0
    return abs(imag(ω)) / denom
end

@inline function matching_controls(s, l, m, a, omega; boundary_condition=IN)
    absomega = abs(omega)
    small_real_omega = isreal(omega) && absomega < 0.01
    sfe = (small_real_omega && abs(s) == 2) ? 1 : 0
    lfe = 0
    N = max(Int64(floor(20 - 20 * log10(1 - abs(a)) - 2 * log(l))), 15)
    xm = _DEFAULT_XM
    if boundary_condition == IN && s == 1 && isreal(omega)
        N = max(N, 50)
        xm = -0.6
        TSinInf = 0
    elseif boundary_condition == IN && s == 2 && isreal(omega) && iszero(a) && l == 2
        if m == 2
            N = 55
            xm = -0.6
            sfe = absomega < 0.1 ? 1 : sfe
            TSinInf = 0
        elseif m == 1
            N = 60
            xm = -0.6
            TSinInf = 0
        elseif m == 0
            N = 80
            xm = -0.8
            TSinInf = 0
        elseif m == -1
            N = 25
            xm = -0.4
            TSinInf = 0
        else
            N = 40
            xm = -0.5
            TSinInf = 0
        end
    elseif boundary_condition == IN && s == 2 && isreal(omega) && l == 2 && m != 2
        N = 10
        xm = -0.9
        TSinInf = 0
    elseif boundary_condition == IN && s == 2 && isreal(omega)
        # Direct +2 infinity TSI is inaccurate for small/moderate real frequencies.
        N = max(N, 50)
        xm = -0.6
        TSinInf = 0
    else
        TSinInf = s == 2 && isreal(omega) ? 1 : 0
    end
    TSoutInf = s == -2 && isreal(omega) ? 1 : 0
    TSinHor = 0
    TSoutHor = 0
    return (
        xm = xm,
        rhom = _DEFAULT_RHOM,
        N = N,
        tol = _DEFAULT_TOL,
        sfe = sfe,
        lfe = lfe,
        TSinInf = TSinInf,
        TSoutInf = TSoutInf,
        TSinHor = TSinHor,
        TSoutHor = TSoutHor,
        boundary_condition = boundary_condition,
    )
end

@inline function _resolve_option(user_value, default_value)
    return user_value === nothing ? default_value : user_value
end

@inline function _validate_y_branch(s, boundary_condition)
    if !(s == -2 && boundary_condition == IN || s == 2 && boundary_condition == UP)
        error("Y_radial currently supports only s = -2 with IN and s = +2 with UP.")
    end
end

@inline function _y_phase_argument(r, a, m, omega)
    rp = 1 + sqrt(1 - a^2)
    rm = 1 - sqrt(1 - a^2)
    rs = rstar_from_r(a, r)
    return omega * rs + a * m * log((r - rm) / (r - rp)) / (rp - rm)
end

@inline function _y_phase_derivative(r, a, m, omega)
    return ((r^2 + a^2) * omega - a * m) / (r^2 - 2 * r + a^2)
end

@inline function _y_branch_transmission_amplitude(s, m, a, omega, lambda)
    if s == -2
        return ConversionFactors.Btrans(s, m, a, omega, lambda)
    elseif s == 2
        return ConversionFactors.Ctrans(s, m, a, omega, lambda)
    else
        error("Unsupported Y branch.")
    end
end

struct YRadialFunction
    mode::Mode
    boundary_condition::BoundaryCondition
    transmission_amplitude
    incidence_amplitude
    reflection_amplitude
    P_solution::Union{Function, Missing}
    Teukolsky_solution::Union{Function, Missing}
    X_solution::Union{Function, Missing}
    Y_scalar_solution::Function
    Y_solution::Function
    normalization_convention::NormalizationConvention
end

function Base.show(io::IO, ::MIME"text/plain", y_func::YRadialFunction)
    println(io, "YRadialFunction(")
    print(io, "    mode = "); show(io, "text/plain", y_func.mode); println(io, ",")
    println(io, "    boundary_condition = $(y_func.boundary_condition),")
    println(io, "    transmission_amplitude = $(y_func.transmission_amplitude),")
    println(io, "    incidence_amplitude = $(y_func.incidence_amplitude),")
    println(io, "    reflection_amplitude = $(y_func.reflection_amplitude),")
    println(io, "    normalization_convention = $(y_func.normalization_convention)")
    print(io, ")")
end

function Base.show(io::IO, y_func::YRadialFunction)
    print(io, "YRadialFunction(mode = Mode(s = $(y_func.mode.s), l = $(y_func.mode.l), m = $(y_func.mode.m), a = $(y_func.mode.a), omega = $(y_func.mode.omega), lambda = $(y_func.mode.lambda)), boundary_condition = $(y_func.boundary_condition))")
end

@inline function _isem_parameters(s, l, m, a, omega)
    params = Matching.Parameters.isem_parameters(Matching.Parameters.TeukolskyParameters(s, l, m, a, omega))
    return params.kappa, params.epsilon, params.tau, params.lambda
end

@inline function _rho_parameters(s, l, m, a, omega)
    params = Matching.Parameters.rho_parameters(Matching.Parameters.TeukolskyParameters(s, l, m, a, omega))
    return params.kappa, params.epsilon, params.tau, params.lambda, params.z
end

@inline _kerr_delta(r, a) = r^2 + a^2 - 2 * r
@inline _kerr_rplus(a) = 1 + sqrt(1 - a^2)
@inline _kerr_omega_horizon(a) = a / (2 * _kerr_rplus(a))

@inline function _eval_inf_ansatz(coeffs, omega, r)
    zinv = inv(omega * r)
    zpow = one(zinv)
    ans = zero(eltype(coeffs))
    @inbounds for c in coeffs
        ans += c * zpow
        zpow *= zinv
    end
    return ans
end

@inline function _eval_dinf_ansatz_dr(coeffs, omega, r)
    ans = zero(eltype(coeffs))
    omega_inv = inv(omega)
    omega_pow = omega_inv
    r_pow_inv = inv(r^2)
    r_inv = inv(r)
    @inbounds for i in 1:(length(coeffs) - 1)
        ans -= i * coeffs[i + 1] * omega_pow * r_pow_inv
        omega_pow *= omega_inv
        r_pow_inv *= r_inv
    end
    return ans
end

@inline function _eval_hor_ansatz(coeffs, dr)
    drpow = one(dr)
    ans = zero(eltype(coeffs))
    @inbounds for c in coeffs
        ans += c * drpow
        drpow *= dr
    end
    return ans
end

@inline function _eval_dhor_ansatz_dr(coeffs, dr)
    ans = zero(eltype(coeffs))
    drpow = one(dr)
    @inbounds for i in 1:(length(coeffs) - 1)
        ans += i * coeffs[i + 1] * drpow
        drpow *= dr
    end
    return ans
end

@inline function _lazy_coefficients!(coeffs, builder)
    if coeffs[] === nothing
        coeffs[] = builder()
    end
    return coeffs[]
end

@inline function _gsn_horizon_ansatz(mode::Mode, boundary_condition::BoundaryCondition, transmission_amplitude, incidence_amplitude, reflection_amplitude, rs; order = _GSN_HORIZON_ASYMPTOTIC_ORDER)
    s = mode.s
    m = mode.m
    a = mode.a
    omega = mode.omega
    lambda = mode.lambda
    r = r_from_rstar(a, rs)
    p = omega - m * _kerr_omega_horizon(a)

    gin(r0) = InitialConditions.gansatz(
        ord -> AsymptoticExpansionCoefficients.ingoing_coefficient_at_hor(s, m, a, omega, lambda, ord),
        a,
        r0;
        order = order,
    )
    dgin_dr(r0) = InitialConditions.dgansatz_dr(
        ord -> AsymptoticExpansionCoefficients.ingoing_coefficient_at_hor(s, m, a, omega, lambda, ord),
        a,
        r0;
        order = order,
    )
    gout(r0) = InitialConditions.gansatz(
        ord -> AsymptoticExpansionCoefficients.outgoing_coefficient_at_hor(s, m, a, omega, lambda, ord),
        a,
        r0;
        order = order,
    )
    dgout_dr(r0) = InitialConditions.dgansatz_dr(
        ord -> AsymptoticExpansionCoefficients.outgoing_coefficient_at_hor(s, m, a, omega, lambda, ord),
        a,
        r0;
        order = order,
    )

    Δfac = _kerr_delta(r, a) / (r^2 + a^2)
    if boundary_condition == IN
        phase = exp(-1im * p * rs)
        X = transmission_amplitude * gin(r) * phase
        Xp = transmission_amplitude * (Δfac * dgin_dr(r) - 1im * p * gin(r)) * phase
        return X, Xp
    elseif boundary_condition == UP
        phase_out = exp(1im * p * rs)
        phase_in = exp(-1im * p * rs)
        X = incidence_amplitude * gout(r) * phase_out + reflection_amplitude * gin(r) * phase_in
        Xp = incidence_amplitude * (Δfac * dgout_dr(r) + 1im * p * gout(r)) * phase_out +
             reflection_amplitude * (Δfac * dgin_dr(r) - 1im * p * gin(r)) * phase_in
        return X, Xp
    end
    return nothing
end

@inline function _gsn_infinity_ansatz(mode::Mode, boundary_condition::BoundaryCondition, transmission_amplitude, incidence_amplitude, reflection_amplitude, rs; order = _GSN_INFINITY_ASYMPTOTIC_ORDER)
    s = mode.s
    m = mode.m
    a = mode.a
    omega = mode.omega
    lambda = mode.lambda
    r = r_from_rstar(a, rs)

    fin(r0) = InitialConditions.fansatz(
        ord -> AsymptoticExpansionCoefficients.ingoing_coefficient_at_inf(s, m, a, omega, lambda, ord),
        omega,
        r0;
        order = order,
    )
    dfin_dr(r0) = InitialConditions.dfansatz_dr(
        ord -> AsymptoticExpansionCoefficients.ingoing_coefficient_at_inf(s, m, a, omega, lambda, ord),
        omega,
        r0;
        order = order,
    )
    fout(r0) = InitialConditions.fansatz(
        ord -> AsymptoticExpansionCoefficients.outgoing_coefficient_at_inf(s, m, a, omega, lambda, ord),
        omega,
        r0;
        order = order,
    )
    dfout_dr(r0) = InitialConditions.dfansatz_dr(
        ord -> AsymptoticExpansionCoefficients.outgoing_coefficient_at_inf(s, m, a, omega, lambda, ord),
        omega,
        r0;
        order = order,
    )

    Δfac = _kerr_delta(r, a) / (r^2 + a^2)
    if boundary_condition == UP
        phase = exp(1im * omega * rs)
        X = transmission_amplitude * fout(r) * phase
        Xp = transmission_amplitude * (Δfac * dfout_dr(r) + 1im * omega * fout(r)) * phase
        return X, Xp
    elseif boundary_condition == IN
        phase_out = exp(1im * omega * rs)
        phase_in = exp(-1im * omega * rs)
        X = reflection_amplitude * fout(r) * phase_out + incidence_amplitude * fin(r) * phase_in
        Xp = reflection_amplitude * (Δfac * dfout_dr(r) + 1im * omega * fout(r)) * phase_out +
             incidence_amplitude * (Δfac * dfin_dr(r) - 1im * omega * fin(r)) * phase_in
        return X, Xp
    end
    return nothing
end

@inline function _relative_mismatch(a, b)
    scale = max(abs(a), abs(b), eps(float(real(abs(a) + abs(b)))))
    return abs(a - b) / scale
end

@inline function _continuous_gsn_patch_or_raw(gsn_solution, rs, X, Xp)
    raw = gsn_solution(rs)
    mismatch = max(
        _relative_mismatch(X, raw[1]),
        _relative_mismatch(Xp, raw[2]),
    )
    mismatch > _GSN_ASYMPTOTIC_PATCH_RELATIVE_MISMATCH_MAX && return raw
    return (X, Xp, 0.0)
end

function _with_gsn_asymptotic_patches(gsn_solution, mode::Mode, boundary_condition::BoundaryCondition, transmission_amplitude, incidence_amplitude, reflection_amplitude; horizon_delta_r_max = _GSN_HORIZON_DELTA_R_MAX, infinity_phase_min = _GSN_INFINITY_PHASE_MIN)
    _is_omega_complex(mode.omega) && return gsn_solution
    s = mode.s
    m = mode.m
    a = mode.a
    omega = mode.omega
    lambda = mode.lambda
    rp = _kerr_rplus(a)
    p = omega - m * _kerr_omega_horizon(a)

    hor_order = _GSN_HORIZON_ASYMPTOTIC_ORDER
    inf_order = _GSN_INFINITY_ASYMPTOTIC_ORDER
    real_type = promote_type(typeof(float(real(a))), typeof(float(real(omega))), typeof(float(real(lambda))))
    coeff_type = Complex{real_type}
    hor_in_coeffs = Ref{Union{Nothing, Vector{coeff_type}}}(nothing)
    hor_out_coeffs = Ref{Union{Nothing, Vector{coeff_type}}}(nothing)
    inf_in_coeffs = Ref{Union{Nothing, Vector{coeff_type}}}(nothing)
    inf_out_coeffs = Ref{Union{Nothing, Vector{coeff_type}}}(nothing)
    get_hor_in_coeffs() = _lazy_coefficients!(hor_in_coeffs, () -> coeff_type[AsymptoticExpansionCoefficients.ingoing_coefficient_at_hor(s, m, a, omega, lambda, ord) for ord in 0:hor_order])
    get_hor_out_coeffs() = _lazy_coefficients!(hor_out_coeffs, () -> coeff_type[AsymptoticExpansionCoefficients.outgoing_coefficient_at_hor(s, m, a, omega, lambda, ord) for ord in 0:hor_order])
    get_inf_in_coeffs() = _lazy_coefficients!(inf_in_coeffs, () -> coeff_type[AsymptoticExpansionCoefficients.ingoing_coefficient_at_inf(s, m, a, omega, lambda, ord) for ord in 0:inf_order])
    get_inf_out_coeffs() = _lazy_coefficients!(inf_out_coeffs, () -> coeff_type[AsymptoticExpansionCoefficients.outgoing_coefficient_at_inf(s, m, a, omega, lambda, ord) for ord in 0:inf_order])

    return rs -> begin
        r = r_from_rstar(a, rs)
        Δfac = _kerr_delta(r, a) / (r^2 + a^2)
        if r - rp <= horizon_delta_r_max
            dr = r - rp
            coeffs_in = get_hor_in_coeffs()
            gin = _eval_hor_ansatz(coeffs_in, dr)
            dgin = _eval_dhor_ansatz_dr(coeffs_in, dr)
            if boundary_condition == IN
                phase = exp(-1im * p * rs)
                X = transmission_amplitude * gin * phase
                Xp = transmission_amplitude * (Δfac * dgin - 1im * p * gin) * phase
                return _continuous_gsn_patch_or_raw(gsn_solution, rs, X, Xp)
            elseif boundary_condition == UP
                coeffs_out = get_hor_out_coeffs()
                gout = _eval_hor_ansatz(coeffs_out, dr)
                dgout = _eval_dhor_ansatz_dr(coeffs_out, dr)
                phase_out = exp(1im * p * rs)
                phase_in = exp(-1im * p * rs)
                X = incidence_amplitude * gout * phase_out + reflection_amplitude * gin * phase_in
                Xp = incidence_amplitude * (Δfac * dgout + 1im * p * gout) * phase_out +
                     reflection_amplitude * (Δfac * dgin - 1im * p * gin) * phase_in
                return _continuous_gsn_patch_or_raw(gsn_solution, rs, X, Xp)
            end
        elseif abs(omega * r) >= infinity_phase_min
            coeffs_out = get_inf_out_coeffs()
            fout = _eval_inf_ansatz(coeffs_out, omega, r)
            dfout = _eval_dinf_ansatz_dr(coeffs_out, omega, r)
            if boundary_condition == UP
                phase = exp(1im * omega * rs)
                X = transmission_amplitude * fout * phase
                Xp = transmission_amplitude * (Δfac * dfout + 1im * omega * fout) * phase
                return _continuous_gsn_patch_or_raw(gsn_solution, rs, X, Xp)
            elseif boundary_condition == IN
                coeffs_in = get_inf_in_coeffs()
                fin = _eval_inf_ansatz(coeffs_in, omega, r)
                dfin = _eval_dinf_ansatz_dr(coeffs_in, omega, r)
                phase_out = exp(1im * omega * rs)
                phase_in = exp(-1im * omega * rs)
                X = reflection_amplitude * fout * phase_out + incidence_amplitude * fin * phase_in
                Xp = reflection_amplitude * (Δfac * dfout + 1im * omega * fout) * phase_out +
                     incidence_amplitude * (Δfac * dfin - 1im * omega * fin) * phase_in
                return _continuous_gsn_patch_or_raw(gsn_solution, rs, X, Xp)
            end
        end
        return gsn_solution(rs)
    end
end

@inline function _teukolsky_solution_real(P, s, a, omega, m, lambda, kappa, epsilon, tau)
    T, Tp, _ = Matching.TeukolskyTransformation.Tx(s, epsilon, kappa, tau)
    function Rsoln(r)
        x = Matching.TeukolskyTransformation.r_to_x(r, kappa)
        Pval, Px, Pxx, error = P(x)
        R = T(x) * Pval
        Rp = - (Tp(x) * Pval + T(x) * Px) / (2 * kappa)
        Rpp = Matching.TeukolskyTransformation.d2R(r, R, Rp, s, a, omega, m, lambda)
        return (R, Rp, Rpp, error)
    end
    return Rsoln
end

@inline function _teukolsky_solution_contour(P, s, a, omega, m, lambda, kappa, epsilon, tau, z)
    T, Tp, _ = Matching.TeukolskyTransformation.Tx(s, epsilon, kappa, tau)
    function Rsoln(r)
        x = Matching.TeukolskyTransformation.r_to_x(r, kappa)
        rho = - x * z
        r = Matching.TeukolskyTransformation.x_to_r(x, kappa)
        Pval, Pρ, Pρρ, error = P(rho)
        Px = -z * Pρ
        R = T(x) * Pval
        Rp = - (Tp(x) * Pval + T(x) * Px) / (2 * kappa)
        Rpp = Matching.TeukolskyTransformation.d2R(r, R, Rp, s, a, omega, m, lambda)
        return (R, Rp, Rpp, error)
    end
    return Rsoln
end

@inline function _build_teukolsky_function(s, l, m, a, omega, boundary_condition; xm=nothing, rhom=nothing, N=nothing, tol=nothing, sfe=nothing, lfe=nothing, TSinInf=nothing, TSoutInf=nothing, TSinHor=nothing, TSoutHor=nothing, return_metadata=false)
    if boundary_condition == IN
        P, (incidence_amplitude, reflection_amplitude), metadata = _Pin(s, l, m, a, omega; xm=xm, rhom=rhom, N=N, tol=tol, sfe=sfe, lfe=lfe, TSinInf=TSinInf, TSoutInf=TSoutInf, TSinHor=TSinHor, components=true)
    elseif boundary_condition == UP
        P, (incidence_amplitude, reflection_amplitude), metadata = _Pup(s, l, m, a, omega; xm=xm, rhom=rhom, N=N, tol=tol, sfe=sfe, lfe=lfe, TSoutInf=TSoutInf, TSinHor=TSinHor, TSoutHor=TSoutHor, components=true)
    else
        error("ISEM currently provides the IN and UP basis solutions only.")
    end

    if _is_omega_complex(omega)
        kappa, epsilon, tau, lambda, z = _rho_parameters(s, l, m, a, omega)
        Rsoln = _teukolsky_solution_contour(P, s, a, omega, m, lambda, kappa, epsilon, tau, z)
        mode = Mode(s, l, m, a, omega, lambda)
        teuk_func = TeukolskyRadialFunction(mode, boundary_condition, one(incidence_amplitude), incidence_amplitude, reflection_amplitude, P, missing, Rsoln, UNIT_TEUKOLSKY_TRANS)
        return return_metadata ? (teuk_func, metadata) : teuk_func
    end

    kappa, epsilon, tau, lambda = _isem_parameters(s, l, m, a, omega)
    Rsoln = _teukolsky_solution_real(P, s, a, omega, m, lambda, kappa, epsilon, tau)
    mode = Mode(s, l, m, a, omega, lambda)
    teuk_func = TeukolskyRadialFunction(mode, boundary_condition, one(incidence_amplitude), incidence_amplitude, reflection_amplitude, P, missing, Rsoln, UNIT_TEUKOLSKY_TRANS)
    return return_metadata ? (teuk_func, metadata) : teuk_func
end

@inline function _combine_down(Rin::TeukolskyRadialFunction, Rup::TeukolskyRadialFunction)
    Binc = Rin.incidence_amplitude
    Bref = Rin.reflection_amplitude
    Cinc = Rup.incidence_amplitude
    Cref = Rup.reflection_amplitude

    function teuk_func(coord)
        Rin_vals = Rin.Teukolsky_solution(coord)
        Rup_vals = Rup.Teukolsky_solution(coord)
        R = Binc^-1 * (Rin_vals[1] - Bref * Rup_vals[1])
        dR = Binc^-1 * (Rin_vals[2] - Bref * Rup_vals[2])
        d2R = Binc^-1 * (Rin_vals[3] - Bref * Rup_vals[3])
        err = max(Rin_vals[4], Rup_vals[4])
        return (R, dR, d2R, err)
    end

    return TeukolskyRadialFunction(
        Rin.mode,
        DOWN,
        one(Binc),
        Binc^-1 - Bref * Cref / Binc,
        -Bref * Cinc / Binc,
        missing,
        missing,
        teuk_func,
        UNIT_TEUKOLSKY_TRANS,
    )
end

@inline function _combine_out(Rin::TeukolskyRadialFunction, Rup::TeukolskyRadialFunction)
    Binc = Rin.incidence_amplitude
    Cinc = Rup.incidence_amplitude
    Cref = Rup.reflection_amplitude

    function teuk_func(coord)
        Rin_vals = Rin.Teukolsky_solution(coord)
        Rup_vals = Rup.Teukolsky_solution(coord)
        R = Cinc^-1 * (Rup_vals[1] - Cref * Rin_vals[1])
        dR = Cinc^-1 * (Rup_vals[2] - Cref * Rin_vals[2])
        d2R = Cinc^-1 * (Rup_vals[3] - Cref * Rin_vals[3])
        err = max(Rin_vals[4], Rup_vals[4])
        return (R, dR, d2R, err)
    end

    return TeukolskyRadialFunction(
        Rin.mode,
        OUT,
        one(Cinc),
        Cinc^-1 - Rin.reflection_amplitude * Cref / Cinc,
        -Binc * Cref / Cinc,
        missing,
        missing,
        teuk_func,
        UNIT_TEUKOLSKY_TRANS,
    )
end

@inline function _gsn_amplitudes_from_teukolsky(teuk_func::TeukolskyRadialFunction)
    s = teuk_func.mode.s
    m = teuk_func.mode.m
    a = teuk_func.mode.a
    omega = teuk_func.mode.omega
    lambda = teuk_func.mode.lambda

    if teuk_func.boundary_condition == IN
        trans_factor = ConversionFactors.Btrans(s, m, a, omega, lambda)
        incidence_factor = ConversionFactors.Binc(s, m, a, omega, lambda)
        reflection_factor = ConversionFactors.Bref(s, m, a, omega, lambda)
    elseif teuk_func.boundary_condition == UP
        trans_factor = ConversionFactors.Ctrans(s, m, a, omega, lambda)
        incidence_factor = ConversionFactors.Cinc(s, m, a, omega, lambda)
        reflection_factor = ConversionFactors.Cref(s, m, a, omega, lambda)
    elseif teuk_func.boundary_condition == DOWN
        trans_factor = ConversionFactors.Binc(s, m, a, omega, lambda)
        incidence_factor = ConversionFactors.Btrans(s, m, a, omega, lambda)
        reflection_factor = ConversionFactors.Cinc(s, m, a, omega, lambda)
    elseif teuk_func.boundary_condition == OUT
        trans_factor = ConversionFactors.Cinc(s, m, a, omega, lambda)
        incidence_factor = ConversionFactors.Bref(s, m, a, omega, lambda)
        reflection_factor = ConversionFactors.Binc(s, m, a, omega, lambda)
    else
        error("Does not understand the boundary condition applied to the solution")
    end

    transmission_amplitude = teuk_func.transmission_amplitude / trans_factor
    incidence_amplitude = teuk_func.incidence_amplitude / incidence_factor
    reflection_amplitude = teuk_func.reflection_amplitude / reflection_factor
    return transmission_amplitude, incidence_amplitude, reflection_amplitude
end

@inline _metadata_property(metadata, key::Symbol, default) = metadata !== missing && hasproperty(metadata, key) ? getproperty(metadata, key) : default

@inline function _isem_parameter_summary(; xm, rhom, N, tol, sfe, lfe, TSinInf, TSoutInf, TSinHor, TSoutHor, metadata=missing)
    return (
        xm = _metadata_property(metadata, :xsplit, _metadata_property(metadata, :xm_match, xm)),
        rhom = _metadata_property(metadata, :rhom, rhom),
        N = _metadata_property(metadata, :N, N),
        tol = tol,
        sfe = sfe,
        lfe = lfe,
        TSinInf = TSinInf,
        TSoutInf = TSoutInf,
        TSinHor = TSinHor,
        TSoutHor = TSoutHor,
        xm_requested = _metadata_property(metadata, :xm, xm),
        xm_min = _metadata_property(metadata, :xm_min, missing),
        xm_max = _metadata_property(metadata, :xm_max, missing),
        xm_match = _metadata_property(metadata, :xm_match, missing),
        split_mismatch = _metadata_property(metadata, :split_mismatch, missing),
        N_candidates = _metadata_property(metadata, :N_candidates, missing),
    )
end

@inline function _build_gsn_function_from_teukolsky(teuk_func::TeukolskyRadialFunction; isem_parameters = missing, use_asymptotic_patches = true, gsn_horizon_delta_r_max = _GSN_HORIZON_DELTA_R_MAX, gsn_infinity_phase_min = _GSN_INFINITY_PHASE_MIN)
    mode = teuk_func.mode
    teukolsky_from_gsn_matrix = r -> begin
        if r isa BigFloat
            return Solutions.Teukolsky_radial_function_from_Sasaki_Nakamura_function_conversion_matrix(mode.s, mode.m, big(mode.a), big(mode.omega), big(mode.lambda), r)
        end
        return Solutions.Teukolsky_radial_function_from_Sasaki_Nakamura_function_conversion_matrix(mode.s, mode.m, mode.a, mode.omega, mode.lambda, r)
    end
    if teuk_func.P_solution !== missing && !_is_omega_complex(mode.omega)
        gsn_solution = Matching.TeukolskyTransformation.P_to_GSN_solution_from_matrix(rs -> r_from_rstar(mode.a, rs), teukolsky_from_gsn_matrix, teuk_func.P_solution, mode.s, mode.m, mode.a, mode.omega, mode.lambda)
    else
        teuk_to_gsn_matrix = r -> inv(teukolsky_from_gsn_matrix(r))
        gsn_solution = Matching.TeukolskyTransformation.Teukolsky_to_GSN_solution(rs -> r_from_rstar(mode.a, rs), teuk_to_gsn_matrix, teuk_func.Teukolsky_solution)
    end
    transmission_amplitude, incidence_amplitude, reflection_amplitude = _gsn_amplitudes_from_teukolsky(teuk_func)
    gsn_transmission_amplitude = transmission_amplitude
    raw_gsn_solution = gsn_solution
    gsn_solution_unit_trans = rs -> begin
        vals = raw_gsn_solution(rs)
        return vals ./ gsn_transmission_amplitude
    end
    incidence_amplitude /= gsn_transmission_amplitude
    reflection_amplitude /= gsn_transmission_amplitude
    transmission_amplitude = one(gsn_transmission_amplitude)
    gsn_solution = use_asymptotic_patches ?
        _with_gsn_asymptotic_patches(
            gsn_solution_unit_trans,
            mode,
            teuk_func.boundary_condition,
            transmission_amplitude,
            incidence_amplitude,
            reflection_amplitude,
            horizon_delta_r_max = gsn_horizon_delta_r_max,
            infinity_phase_min = gsn_infinity_phase_min,
        ) :
        gsn_solution_unit_trans

    return GSNRadialFunction(
        mode,
        teuk_func.boundary_condition,
        missing,
        missing,
        missing,
        missing,
        missing,
        transmission_amplitude,
        incidence_amplitude,
        reflection_amplitude,
        isem_parameters,
        missing,
        gsn_solution,
        UNIT_GSN_TRANS,
        "ISEM",
    )
end

function Teukolsky_radial(s::Int, l::Int, m::Int, a, omega, boundary_condition::BoundaryCondition; xm=nothing, rhom=nothing, N=nothing, tol=nothing, sfe=nothing, lfe=nothing, TSinInf=nothing, TSoutInf=nothing, TSinHor=nothing, TSoutHor=nothing)
    if boundary_condition == DOWN || boundary_condition == OUT
        Rin = Teukolsky_radial(s, l, m, a, omega, IN; xm=xm, rhom=rhom, N=N, tol=tol, sfe=sfe, lfe=lfe, TSinInf=TSinInf, TSoutInf=TSoutInf, TSinHor=TSinHor, TSoutHor=TSoutHor)
        Rup = Teukolsky_radial(s, l, m, a, omega, UP; xm=xm, rhom=rhom, N=N, tol=tol, sfe=sfe, lfe=lfe, TSinInf=TSinInf, TSoutInf=TSoutInf, TSinHor=TSinHor, TSoutHor=TSoutHor)
        return boundary_condition == DOWN ? _combine_down(Rin, Rup) : _combine_out(Rin, Rup)
    end

    defaults = matching_controls(s, l, m, a, omega; boundary_condition = boundary_condition)
    xm = _resolve_option(xm, defaults.xm)
    rhom = _resolve_option(rhom, defaults.rhom)
    N = _resolve_option(N, defaults.N)
    tol = _resolve_option(tol, defaults.tol)
    sfe = _resolve_option(sfe, defaults.sfe)
    lfe = _resolve_option(lfe, defaults.lfe)
    TSinInf = _resolve_option(TSinInf, defaults.TSinInf)
    TSoutInf = _resolve_option(TSoutInf, defaults.TSoutInf)
    TSinHor = _resolve_option(TSinHor, defaults.TSinHor)
    TSoutHor = _resolve_option(TSoutHor, defaults.TSoutHor)
    return _build_teukolsky_function(s, l, m, a, omega, boundary_condition; xm=xm, rhom=rhom, N=N, tol=tol, sfe=sfe, lfe=lfe, TSinInf=TSinInf, TSoutInf=TSoutInf, TSinHor=TSinHor, TSoutHor=TSoutHor)
end

function GSN_radial(s::Int, l::Int, m::Int, a, omega, boundary_condition::BoundaryCondition; xm=nothing, rhom=nothing, N=nothing, tol=nothing, sfe=nothing, lfe=nothing, TSinInf=nothing, TSoutInf=nothing, TSinHor=nothing, TSoutHor=nothing, use_gsn_asymptotic_patches=false, gsn_horizon_delta_r_max=_GSN_HORIZON_DELTA_R_MAX, gsn_infinity_phase_min=_GSN_INFINITY_PHASE_MIN)
    defaults = matching_controls(s, l, m, a, omega; boundary_condition = boundary_condition)
    xm = _resolve_option(xm, defaults.xm)
    rhom = _resolve_option(rhom, defaults.rhom)
    N = _resolve_option(N, defaults.N)
    tol = _resolve_option(tol, defaults.tol)
    sfe = _resolve_option(sfe, defaults.sfe)
    lfe = _resolve_option(lfe, defaults.lfe)
    TSinInf = _resolve_option(TSinInf, defaults.TSinInf)
    TSoutInf = _resolve_option(TSoutInf, defaults.TSoutInf)
    TSinHor = _resolve_option(TSinHor, defaults.TSinHor)
    TSoutHor = _resolve_option(TSoutHor, defaults.TSoutHor)
    teuk_func, metadata = _build_teukolsky_function(s, l, m, a, omega, boundary_condition; xm=xm, rhom=rhom, N=N, tol=tol, sfe=sfe, lfe=lfe, TSinInf=TSinInf, TSoutInf=TSoutInf, TSinHor=TSinHor, TSoutHor=TSoutHor, return_metadata=true)
    isem_parameters = _isem_parameter_summary(; xm=xm, rhom=rhom, N=N, tol=tol, sfe=sfe, lfe=lfe, TSinInf=TSinInf, TSoutInf=TSoutInf, TSinHor=TSinHor, TSoutHor=TSoutHor, metadata=metadata)
    return _build_gsn_function_from_teukolsky(teuk_func; isem_parameters=isem_parameters, use_asymptotic_patches=use_gsn_asymptotic_patches, gsn_horizon_delta_r_max=gsn_horizon_delta_r_max, gsn_infinity_phase_min=gsn_infinity_phase_min)
end

function GSN_radial(s::Int, l::Int, m::Int, a, omega; xm=_DEFAULT_XM, rhom=_DEFAULT_RHOM, N=_DEFAULT_N, tol=_DEFAULT_TOL, sfe=nothing, lfe=nothing, TSinInf=nothing, TSoutInf=nothing, TSinHor=nothing, TSoutHor=nothing, use_gsn_asymptotic_patches=false, gsn_horizon_delta_r_max=_GSN_HORIZON_DELTA_R_MAX, gsn_infinity_phase_min=_GSN_INFINITY_PHASE_MIN)
    Xin = GSN_radial(s, l, m, a, omega, IN; xm=xm, rhom=rhom, N=N, tol=tol, sfe=sfe, lfe=lfe, TSinInf=TSinInf, TSoutInf=TSoutInf, TSinHor=TSinHor, TSoutHor=TSoutHor, use_gsn_asymptotic_patches=use_gsn_asymptotic_patches, gsn_horizon_delta_r_max=gsn_horizon_delta_r_max, gsn_infinity_phase_min=gsn_infinity_phase_min)
    Xup = GSN_radial(s, l, m, a, omega, UP; xm=xm, rhom=rhom, N=N, tol=tol, sfe=sfe, lfe=lfe, TSinInf=TSinInf, TSoutInf=TSoutInf, TSinHor=TSinHor, TSoutHor=TSoutHor, use_gsn_asymptotic_patches=use_gsn_asymptotic_patches, gsn_horizon_delta_r_max=gsn_horizon_delta_r_max, gsn_infinity_phase_min=gsn_infinity_phase_min)
    return Xin, Xup
end

function Teukolsky_radial(s::Int, l::Int, m::Int, a, omega; xm=_DEFAULT_XM, rhom=_DEFAULT_RHOM, N=_DEFAULT_N, tol=_DEFAULT_TOL, sfe=nothing, lfe=nothing, TSinInf=nothing, TSoutInf=nothing, TSinHor=nothing, TSoutHor=nothing)
    Rin = Teukolsky_radial(s, l, m, a, omega, IN; xm=xm, rhom=rhom, N=N, tol=tol, sfe=sfe, lfe=lfe, TSinInf=TSinInf, TSoutInf=TSoutInf, TSinHor=TSinHor, TSoutHor=TSoutHor)
    Rup = Teukolsky_radial(s, l, m, a, omega, UP; xm=xm, rhom=rhom, N=N, tol=tol, sfe=sfe, lfe=lfe, TSinInf=TSinInf, TSoutInf=TSoutInf, TSinHor=TSinHor, TSoutHor=TSoutHor)
    return Rin, Rup
end

@inline function _build_y_solution(s, l, m, a, omega, boundary_condition; xm=nothing, rhom=nothing, N=nothing, tol=nothing, sfe=nothing, lfe=nothing, TSinInf=nothing, TSoutInf=nothing, TSinHor=nothing, TSoutHor=nothing)
    _validate_y_branch(s, boundary_condition)
    iszero(omega) && error("Y_radial currently requires nonzero omega.")
    isreal(omega) || error("Y_radial currently supports only real omega.")

    teuk_func = Teukolsky_radial(s, l, m, a, omega, boundary_condition; xm=xm, rhom=rhom, N=N, tol=tol, sfe=sfe, lfe=lfe, TSinInf=TSinInf, TSoutInf=TSoutInf, TSinHor=TSinHor, TSoutHor=TSoutHor)
    lambda = teuk_func.mode.lambda
    kappa = sqrt(1 - a^2)
    branch_transmission_amplitude = _y_branch_transmission_amplitude(s, m, a, omega, lambda)
    if s == -2 && boundary_condition == IN
        func1, func2, func3, func4, func5 = Matching.TeukolskyTransformation.P_to_Y_coeffs_minus2_in(m, a, omega, lambda, branch_transmission_amplitude)
    elseif s == 2 && boundary_condition == UP
        func1, func2, func3, func4, func5 = Matching.TeukolskyTransformation.P_to_Y_coeffs_plus2_up(m, a, omega, lambda, branch_transmission_amplitude)
    else
        error("Unsupported Y branch.")
    end

    x_of_r = r -> Matching.TeukolskyTransformation.r_to_x(r, kappa)
    pullback_P = r -> teuk_func.P_solution(x_of_r(r))
    Xsoln = r -> begin
        P, Pp, _, _ = pullback_P(r)
        return func4(r) * P + func5(r) * Pp
    end

    function Yscalar(r)
        P, _, _, _ = pullback_P(r)
        return func1(r) * P
    end

    function Ysoln(r)
        P, Pp, _, error = pullback_P(r)
        Y = func1(r) * P
        Yp = func2(r) * P + func3(r) * Pp
        X = func4(r) * P + func5(r) * Pp
        return (Y, Yp, X, error)
    end

    mode = Mode(s, l, m, a, omega, lambda)
    return YRadialFunction(mode, boundary_condition, branch_transmission_amplitude, branch_transmission_amplitude * teuk_func.incidence_amplitude, branch_transmission_amplitude * teuk_func.reflection_amplitude, teuk_func.P_solution, teuk_func.Teukolsky_solution, Xsoln, Yscalar, Ysoln, UNIT_TEUKOLSKY_TRANS)
end

function Y_radial(s::Int, l::Int, m::Int, a, omega, boundary_condition::BoundaryCondition; xm=nothing, rhom=nothing, N=nothing, tol=nothing, sfe=nothing, lfe=nothing, TSinInf=nothing, TSoutInf=nothing, TSinHor=nothing, TSoutHor=nothing)
    defaults = matching_controls(s, l, m, a, omega; boundary_condition = boundary_condition)
    xm = _resolve_option(xm, defaults.xm)
    rhom = _resolve_option(rhom, defaults.rhom)
    N = _resolve_option(N, defaults.N)
    tol = _resolve_option(tol, defaults.tol)
    sfe = _resolve_option(sfe, defaults.sfe)
    lfe = _resolve_option(lfe, defaults.lfe)
    TSinInf = _resolve_option(TSinInf, defaults.TSinInf)
    TSoutInf = _resolve_option(TSoutInf, defaults.TSoutInf)
    TSinHor = _resolve_option(TSinHor, defaults.TSinHor)
    TSoutHor = _resolve_option(TSoutHor, defaults.TSoutHor)
    return _build_y_solution(s, l, m, a, omega, boundary_condition; xm=xm, rhom=rhom, N=N, tol=tol, sfe=sfe, lfe=lfe, TSinInf=TSinInf, TSoutInf=TSoutInf, TSinHor=TSinHor, TSoutHor=TSoutHor)
end

(y_func::YRadialFunction)(coord) = y_func.Y_scalar_solution(coord)

end
