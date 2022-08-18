using SpinWeightedSpheroidalHarmonics

# Mathematica being Mathematica
const Pi = pi
const I = 1im

function _sourceterm_without_phasing_circularorbit_nearhorizon_seriesexpansion_zerothorder(s::Int, l::Int, m::Int, a, omega, En, Lz; swsh_piover2=nothing, psptheta_piover2=nothing, p2sptheta2_piover2=nothing)
    if isnothing(swsh_piover2) || isnothing(psptheta_piover2) || isnothing(p2sptheta2_piover2)
        swsh_piover2 = spin_weighted_spheroidal_harmonic(s, l, m, a*omega, pi/2, 0)
        psptheta_piover2 = spin_weighted_spheroidal_harmonic(s, l, m, a*omega, pi/2, 0; theta_derivative=1)
        p2sptheta2_piover2 = spin_weighted_spheroidal_harmonic(s, l, m, a*omega, pi/2, 0; theta_derivative=2)
    end

    # Some useful terms to pre-compute
    sqrt_one_minus_a_sq = sqrt(1 - a^2)
    horizon_term = (a*Lz - 2*En*(1 + sqrt_one_minus_a_sq))

    #=
    We have derived/shown the explicit expression for
    s = +2 ONLY
    =#
    if s == 2
        # s = +2
        return begin
            -((1/(horizon_term*(1 + sqrt_one_minus_a_sq)*((-a)*Lz*(1 + sqrt_one_minus_a_sq) + En*(4 - 2*a^2 + 4*sqrt_one_minus_a_sq))*Pi))*(4*((-a)*En + Lz)^2*sqrt_one_minus_a_sq*(1 + sqrt_one_minus_a_sq)^2*swsh_piover2*(I*a*horizon_term*m - 2*a*Lz*sqrt_one_minus_a_sq + 4*En*sqrt_one_minus_a_sq*(1 + sqrt_one_minus_a_sq) - 2*I*horizon_term*(1 + sqrt_one_minus_a_sq)*omega) - 
            ((-a)*En + Lz)^2*(1 + sqrt_one_minus_a_sq)^2*swsh_piover2*(8*(-1 + a^2)*horizon_term + a^2*horizon_term*m^2 + 6*I*a*horizon_term*m*sqrt_one_minus_a_sq - 4*a*horizon_term*m*(1 + sqrt_one_minus_a_sq)*omega - 12*I*horizon_term*sqrt_one_minus_a_sq*(1 + sqrt_one_minus_a_sq)*omega - 4*(1 + sqrt_one_minus_a_sq)^2*((-a)*Lz + 2*En*(1 + sqrt_one_minus_a_sq))*omega^2) + 
            I*(1 + sqrt_one_minus_a_sq)^2*(I*a*horizon_term*m - 2*a*Lz*sqrt_one_minus_a_sq + 4*En*sqrt_one_minus_a_sq*(1 + sqrt_one_minus_a_sq) - 2*I*horizon_term*(1 + sqrt_one_minus_a_sq)*omega)*(2*((-a)*En + Lz)^2*swsh_piover2*(a*m - 2*(1 + sqrt_one_minus_a_sq)*omega) - (2*a*En - 2*Lz)*(-horizon_term - a*Lz + 2*En*(1 + sqrt_one_minus_a_sq))*(psptheta_piover2 + swsh_piover2*(m - a*omega))) - 
            ((-a)*Lz + 2*En*(1 + sqrt_one_minus_a_sq))*(((-a)*En + Lz)^2*(1 + sqrt_one_minus_a_sq)^2*swsh_piover2*(a*m - 2*(1 + sqrt_one_minus_a_sq)*omega)*((-a)*m + 2*omega + 2*sqrt_one_minus_a_sq*(I + omega)) + (2*a*En - 2*Lz)*(1 + sqrt_one_minus_a_sq)^2*(-horizon_term - a*Lz + 2*En*(1 + sqrt_one_minus_a_sq))*(a*m - 2*(1 + sqrt_one_minus_a_sq)*omega)*(psptheta_piover2 + swsh_piover2*(m - a*omega)) - 
            (horizon_term + a*Lz - 2*En*(1 + sqrt_one_minus_a_sq))^2*((2 - a^2 + 2*sqrt_one_minus_a_sq)*p2sptheta2_piover2 + (1 + sqrt_one_minus_a_sq)*((1 + sqrt_one_minus_a_sq)*(2*m*psptheta_piover2 + (-2 + m^2)*swsh_piover2) - 2*a*(psptheta_piover2 + m*swsh_piover2)*(I + omega + sqrt_one_minus_a_sq*omega) + a^2*swsh_piover2*omega*(2*I + omega + sqrt_one_minus_a_sq*omega))))))
        end
    else
        # Throw an error, this spin weight is not supported
        throw(DomainError(s, "Currently only spin weight s of +2 is supported"))
    end
end

function _sourceterm_without_phasing_circularorbit_nearhorizon_seriesexpansion_firstorder(s::Int, l::Int, m::Int, a, omega, En, Lz; swsh_piover2=nothing, psptheta_piover2=nothing, p2sptheta2_piover2=nothing)
    if isnothing(swsh_piover2) || isnothing(psptheta_piover2) || isnothing(p2sptheta2_piover2)
        swsh_piover2 = spin_weighted_spheroidal_harmonic(s, l, m, a*omega, pi/2, 0)
        psptheta_piover2 = spin_weighted_spheroidal_harmonic(s, l, m, a*omega, pi/2, 0; theta_derivative=1)
        p2sptheta2_piover2 = spin_weighted_spheroidal_harmonic(s, l, m, a*omega, pi/2, 0; theta_derivative=2)
    end

    # Some useful terms to pre-compute
    sqrt_one_minus_a_sq = sqrt(1 - a^2)
    horizon_term = (a*Lz - 2*En*(1 + sqrt_one_minus_a_sq))

    #=
    We have derived/shown the explicit expression for
    s = +2 ONLY
    =#
    if s == 2
        # s = +2
        return begin
            (1/(4*sqrt_one_minus_a_sq*(a*Lz*(1 + sqrt_one_minus_a_sq) + 2*En*(a^2 - 2*(1 + sqrt_one_minus_a_sq)))^2*Pi))*(8*((-a)*En + Lz)^2*(1 + sqrt_one_minus_a_sq)*(a*Lz*(-1 + a^2 + sqrt_one_minus_a_sq) + 2*En*(6*(1 + sqrt_one_minus_a_sq) + a^4*(4 + sqrt_one_minus_a_sq) - 2*a^2*(5 + 4*sqrt_one_minus_a_sq)))*(swsh_piover2) - 
            (1/horizon_term^5)*(4*((-a)*En + Lz)^2*(-2*horizon_term^5*(1 + sqrt_one_minus_a_sq)*(a*Lz*(-1 + sqrt_one_minus_a_sq) + 2*En*(2 - 6*sqrt_one_minus_a_sq - 8*sqrt_one_minus_a_sq^2 + a^2*(-1 + 4*sqrt_one_minus_a_sq + sqrt_one_minus_a_sq^2)))*((-I)*a*m + 2*sqrt_one_minus_a_sq + 2*I*(1 + sqrt_one_minus_a_sq)*omega) - 2*horizon_term^3*sqrt_one_minus_a_sq*((-a)*Lz*(1 + sqrt_one_minus_a_sq) + En*(4 - 2*a^2 + 4*sqrt_one_minus_a_sq))*
            (-4*horizon_term^2*sqrt_one_minus_a_sq - I*m*(2*a^5 + a^3*(-10 - 4*En^2 + Lz^2 - 6*sqrt_one_minus_a_sq) + 8*a*(1 + En^2)*(1 + sqrt_one_minus_a_sq) + 8*En*Lz*(1 + sqrt_one_minus_a_sq) + a^4*En*Lz*(4 + sqrt_one_minus_a_sq) - 4*a^2*En*Lz*(4 + 3*sqrt_one_minus_a_sq)) - sqrt_one_minus_a_sq*(1 + sqrt_one_minus_a_sq)*(4*horizon_term*En*(1 + sqrt_one_minus_a_sq) + 2*sqrt_one_minus_a_sq*(((-a)*En + Lz)^2 + (1 + sqrt_one_minus_a_sq)^2)) + 
            2*(10*a^3*En*Lz - 14*a*En*Lz*(1 + sqrt_one_minus_a_sq) + (4 + 16*En^2 + Lz^2)*(1 + sqrt_one_minus_a_sq) + a^2*(-7 + En^2*(-7 + sqrt_one_minus_a_sq) - 5*sqrt_one_minus_a_sq + 2*Lz^2*sqrt_one_minus_a_sq) + a^4*(3 + sqrt_one_minus_a_sq - En^2*(5 + sqrt_one_minus_a_sq))) + 4*I*(1 + sqrt_one_minus_a_sq)*(a^4*(1 + En^2) - 8*a*En*Lz*(1 + sqrt_one_minus_a_sq) + (4 + 16*En^2 + Lz^2)*(1 + sqrt_one_minus_a_sq) + a^3*En*Lz*(4 + sqrt_one_minus_a_sq) - 
            a^2*(5 + 3*sqrt_one_minus_a_sq + En^2*(13 + 5*sqrt_one_minus_a_sq)))*omega))*(swsh_piover2)) - (1/horizon_term^3)*(2*((-a)*En + Lz)^2*(horizon_term^3*(1 + sqrt_one_minus_a_sq)*((-a)*Lz + 8*En*(1 + sqrt_one_minus_a_sq) - a^2*En*(4 + sqrt_one_minus_a_sq))*(8*(-1 + a^2) + a^2*m^2 + 6*I*a*m*sqrt_one_minus_a_sq - 4*a*m*(1 + sqrt_one_minus_a_sq)*omega - 12*I*sqrt_one_minus_a_sq*(1 + sqrt_one_minus_a_sq)*omega + 4*(1 + sqrt_one_minus_a_sq)^2*omega^2) + 
            (1/(a^2*Lz^2 - 4*a*En*Lz*(1 + sqrt_one_minus_a_sq) + En^2*(8 - 4*a^2 + 8*sqrt_one_minus_a_sq)))*(horizon_term*((-a)*Lz*(1 + sqrt_one_minus_a_sq) + En*(4 - 2*a^2 + 4*sqrt_one_minus_a_sq))*(8*(1 - a^2)*(1 + sqrt_one_minus_a_sq)*(((-a)*En + Lz)^2*sqrt_one_minus_a_sq + 2*horizon_term*En*(1 + sqrt_one_minus_a_sq) + sqrt_one_minus_a_sq*(1 + sqrt_one_minus_a_sq)^2)*(a^2*Lz^2 - 4*a*En*Lz*(1 + sqrt_one_minus_a_sq) + En^2*(8 - 4*a^2 + 8*sqrt_one_minus_a_sq)) - 
            12*horizon_term*(a*Lz*(-1 + a^2 - sqrt_one_minus_a_sq) + 4*En*(1 + sqrt_one_minus_a_sq) - 2*a^2*En*(2 + sqrt_one_minus_a_sq))*(a^2*Lz^2 - 4*a*En*Lz*(1 + sqrt_one_minus_a_sq) + En^2*(8 - 4*a^2 + 8*sqrt_one_minus_a_sq)) - a*m^2*(a^2*Lz^2 - 4*a*En*Lz*(1 + sqrt_one_minus_a_sq) + En^2*(8 - 4*a^2 + 8*sqrt_one_minus_a_sq))*((-a^5)*(-3 + En^2) + 16*En*Lz*(1 + sqrt_one_minus_a_sq) - a*(-12 + Lz^2)*(1 + sqrt_one_minus_a_sq) + a^4*En*Lz*(4 + sqrt_one_minus_a_sq) - 
            8*a^2*En*Lz*(3 + 2*sqrt_one_minus_a_sq) + a^3*(-15 + 2*Lz^2 - 9*sqrt_one_minus_a_sq + 5*En^2*(1 + sqrt_one_minus_a_sq))) - 4*sqrt_one_minus_a_sq^2*(1 + sqrt_one_minus_a_sq)*(a^2*Lz^2 - 4*a*En*Lz*(1 + sqrt_one_minus_a_sq) + En^2*(8 - 4*a^2 + 8*sqrt_one_minus_a_sq))*(4*horizon_term*En*(1 + sqrt_one_minus_a_sq) + 2*sqrt_one_minus_a_sq*(((-a)*En + Lz)^2 + (1 + sqrt_one_minus_a_sq)^2)) - 
            4*I*m*((-sqrt_one_minus_a_sq)*(1 + sqrt_one_minus_a_sq)*(a^2*Lz^2 - 4*a*En*Lz*(1 + sqrt_one_minus_a_sq) + En^2*(8 - 4*a^2 + 8*sqrt_one_minus_a_sq))*(2*horizon_term*(a*En + Lz*sqrt_one_minus_a_sq) - a*(((-a)*En + Lz)^2*sqrt_one_minus_a_sq + 2*horizon_term*En*(1 + sqrt_one_minus_a_sq) + sqrt_one_minus_a_sq*(1 + sqrt_one_minus_a_sq)^2)) + a*(a^4*Lz^4*(1 + sqrt_one_minus_a_sq) + 8*a^3*En*Lz^3*(a^2 - 2*(1 + sqrt_one_minus_a_sq)) - 32*a*En^3*Lz*(a^4 + 8*(1 + sqrt_one_minus_a_sq) - 4*a^2*(2 + sqrt_one_minus_a_sq)) - 
            24*a^2*En^2*Lz^2*(-4*(1 + sqrt_one_minus_a_sq) + a^2*(3 + sqrt_one_minus_a_sq)) + 16*En^4*(16*(1 + sqrt_one_minus_a_sq) + a^4*(5 + sqrt_one_minus_a_sq) - 4*a^2*(5 + 3*sqrt_one_minus_a_sq)))) - 2*I*m*(a*(1 + sqrt_one_minus_a_sq)*(a^4*Lz^4 - 8*a^3*En*Lz^3*(1 + sqrt_one_minus_a_sq) - 24*a^2*En^2*Lz^2*(a^2 - 2*(1 + sqrt_one_minus_a_sq)) + 16*En^4*(a^4 + 8*(1 + sqrt_one_minus_a_sq) - 4*a^2*(2 + sqrt_one_minus_a_sq)) + 
            32*a*En^3*Lz*(-4*(1 + sqrt_one_minus_a_sq) + a^2*(3 + sqrt_one_minus_a_sq))) + a*sqrt_one_minus_a_sq*(a^2*Lz^2 - 4*a*En*Lz*(1 + sqrt_one_minus_a_sq) + En^2*(8 - 4*a^2 + 8*sqrt_one_minus_a_sq))*(2*horizon_term^2 + sqrt_one_minus_a_sq*(1 + sqrt_one_minus_a_sq)*(a^2*En^2 - 2*a*En*Lz + Lz^2 + (1 + sqrt_one_minus_a_sq)^2) + horizon_term*(2*a^2*En - a*Lz + En*(sqrt_one_minus_a_sq - sqrt_one_minus_a_sq^3)))) - 
            4*m*(1 + sqrt_one_minus_a_sq)*(a^2*Lz^2 - 4*a*En*Lz*(1 + sqrt_one_minus_a_sq) + En^2*(8 - 4*a^2 + 8*sqrt_one_minus_a_sq))*(-3*a*(1 + sqrt_one_minus_a_sq)*(((-a)*En + Lz)^2*sqrt_one_minus_a_sq + 2*horizon_term*En*(1 + sqrt_one_minus_a_sq) + sqrt_one_minus_a_sq*(1 + sqrt_one_minus_a_sq)^2) + 2*horizon_term*(-2*a*horizon_term + a^2*Lz + 5*a*En*(1 + sqrt_one_minus_a_sq) + Lz*sqrt_one_minus_a_sq*(1 + sqrt_one_minus_a_sq) - a^3*En*(4 + sqrt_one_minus_a_sq)))*omega + 
            12*I*(1 + sqrt_one_minus_a_sq)*(256*En^4 - 320*a^2*En^4 + 80*a^4*En^4 - 256*a*En^3*Lz + 256*a^3*En^3*Lz - 32*a^5*En^3*Lz + 96*a^2*En^2*Lz^2 - 72*a^4*En^2*Lz^2 - 16*a^3*En*Lz^3 + 8*a^5*En*Lz^3 + a^4*Lz^4 + 256*En^4*sqrt_one_minus_a_sq - 192*a^2*En^4*sqrt_one_minus_a_sq + 16*a^4*En^4*sqrt_one_minus_a_sq - 256*a*En^3*Lz*sqrt_one_minus_a_sq + 128*a^3*En^3*Lz*sqrt_one_minus_a_sq + 
            96*a^2*En^2*Lz^2*sqrt_one_minus_a_sq - 24*a^4*En^2*Lz^2*sqrt_one_minus_a_sq - 16*a^3*En*Lz^3*sqrt_one_minus_a_sq + a^4*Lz^4*sqrt_one_minus_a_sq + sqrt_one_minus_a_sq*(a^2*Lz^2 - 4*a*En*Lz*(1 + sqrt_one_minus_a_sq) + En^2*(8 - 4*a^2 + 8*sqrt_one_minus_a_sq))*(2*horizon_term^2 + sqrt_one_minus_a_sq*(1 + sqrt_one_minus_a_sq)*(a^2*En^2 - 2*a*En*Lz + Lz^2 + (1 + sqrt_one_minus_a_sq)^2) + horizon_term*(2*a^2*En - a*Lz + En*(sqrt_one_minus_a_sq - sqrt_one_minus_a_sq^3))))*omega - 
            12*(1 + sqrt_one_minus_a_sq)^2*(a^2*Lz^2 - 4*a*En*Lz*(1 + sqrt_one_minus_a_sq) + En^2*(8 - 4*a^2 + 8*sqrt_one_minus_a_sq))*(a^4*(1 + En^2) - 8*a*En*Lz*(1 + sqrt_one_minus_a_sq) + (4 + 16*En^2 + Lz^2)*(1 + sqrt_one_minus_a_sq) + a^3*En*Lz*(4 + sqrt_one_minus_a_sq) - a^2*(5 + 3*sqrt_one_minus_a_sq + En^2*(13 + 5*sqrt_one_minus_a_sq)))*omega^2)))*(swsh_piover2)) + 
            2*((-a)*En + Lz)^2*(a*m - 2*(1 + sqrt_one_minus_a_sq)*omega)*(-((1 + sqrt_one_minus_a_sq)*(a*Lz - 8*En*(1 + sqrt_one_minus_a_sq) + a^2*En*(4 + sqrt_one_minus_a_sq))*((-a)*m + 2*omega + 2*sqrt_one_minus_a_sq*(I + omega))) - 2*(a*Lz*(1 + sqrt_one_minus_a_sq) + 2*En*(a^2 - 2*(1 + sqrt_one_minus_a_sq)))*(-I + 2*(-2 + a^2)*omega - sqrt_one_minus_a_sq*(3*I + 4*omega)))*(swsh_piover2) + 
            (2*I*(2*a*En - 2*Lz)*(I*(1 + sqrt_one_minus_a_sq)*(-2*horizon_term^2*(a*Lz - 8*En*(1 + sqrt_one_minus_a_sq) + a^2*En*(4 + sqrt_one_minus_a_sq)) + (((-a)*En + Lz)^2*sqrt_one_minus_a_sq + 4*horizon_term*En*(1 + sqrt_one_minus_a_sq) + sqrt_one_minus_a_sq*(1 + sqrt_one_minus_a_sq)^2)*((-a)*Lz*(1 + sqrt_one_minus_a_sq) + En*(4 - 2*a^2 + 4*sqrt_one_minus_a_sq)))*(a*m - 2*(1 + sqrt_one_minus_a_sq)*omega) - 
            4*horizon_term^2*((-a)*Lz*(1 + sqrt_one_minus_a_sq) + En*(4 - 2*a^2 + 4*sqrt_one_minus_a_sq))*(2*sqrt_one_minus_a_sq - I*(1 + sqrt_one_minus_a_sq)^2*omega))*((psptheta_piover2) + (m - a*omega)*(swsh_piover2)))/horizon_term + 
            (1/horizon_term^2)*(4*I*((-a)*Lz*(1 + sqrt_one_minus_a_sq) + En*(4 - 2*a^2 + 4*sqrt_one_minus_a_sq))*(-4*horizon_term^2*sqrt_one_minus_a_sq - I*m*(2*a^5 + a^3*(-10 - 4*En^2 + Lz^2 - 6*sqrt_one_minus_a_sq) + 8*a*(1 + En^2)*(1 + sqrt_one_minus_a_sq) + 8*En*Lz*(1 + sqrt_one_minus_a_sq) + a^4*En*Lz*(4 + sqrt_one_minus_a_sq) - 4*a^2*En*Lz*(4 + 3*sqrt_one_minus_a_sq)) - 
            sqrt_one_minus_a_sq*(1 + sqrt_one_minus_a_sq)*(4*horizon_term*En*(1 + sqrt_one_minus_a_sq) + 2*sqrt_one_minus_a_sq*(((-a)*En + Lz)^2 + (1 + sqrt_one_minus_a_sq)^2)) + 2*(10*a^3*En*Lz - 14*a*En*Lz*(1 + sqrt_one_minus_a_sq) + (4 + 16*En^2 + Lz^2)*(1 + sqrt_one_minus_a_sq) + a^2*(-7 + En^2*(-7 + sqrt_one_minus_a_sq) - 5*sqrt_one_minus_a_sq + 2*Lz^2*sqrt_one_minus_a_sq) + a^4*(3 + sqrt_one_minus_a_sq - En^2*(5 + sqrt_one_minus_a_sq))) + 
            4*I*(1 + sqrt_one_minus_a_sq)*(a^4*(1 + En^2) - 8*a*En*Lz*(1 + sqrt_one_minus_a_sq) + (4 + 16*En^2 + Lz^2)*(1 + sqrt_one_minus_a_sq) + a^3*En*Lz*(4 + sqrt_one_minus_a_sq) - a^2*(5 + 3*sqrt_one_minus_a_sq + En^2*(13 + 5*sqrt_one_minus_a_sq)))*omega)*(((-a)*En + Lz)^2*(a*m - 2*(1 + sqrt_one_minus_a_sq)*omega)*(swsh_piover2) + horizon_term*(2*a*En - 2*Lz)*((psptheta_piover2) + (m - a*omega)*(swsh_piover2)))) - 
            2*I*(4*((-a)*En + Lz)^2*(1 + sqrt_one_minus_a_sq)*(a*Lz + En*(a^4 + 4*(1 + sqrt_one_minus_a_sq) - a^2*(7 + 4*sqrt_one_minus_a_sq)))*(a*m - 2*(1 + sqrt_one_minus_a_sq)*omega)*(swsh_piover2) - 4*((-a)*En + Lz)^2*((1 + sqrt_one_minus_a_sq)*(a*Lz*(-1 + sqrt_one_minus_a_sq) + 2*En*(2 - 6*sqrt_one_minus_a_sq - 8*sqrt_one_minus_a_sq^2 + a^2*(-1 + 4*sqrt_one_minus_a_sq + sqrt_one_minus_a_sq^2)))*((-a)*m + 2*(1 + sqrt_one_minus_a_sq)*omega) + 
            2*sqrt_one_minus_a_sq*((-a)*Lz*(1 + sqrt_one_minus_a_sq) + En*(4 - 2*a^2 + 4*sqrt_one_minus_a_sq))*(I*sqrt_one_minus_a_sq + (1 + sqrt_one_minus_a_sq)^2*omega))*(swsh_piover2) + 4*horizon_term*(2*a*En - 2*Lz)*sqrt_one_minus_a_sq*(1 + sqrt_one_minus_a_sq)*(a*Lz - 8*En*(1 + sqrt_one_minus_a_sq) + a^2*En*(4 + sqrt_one_minus_a_sq))*((psptheta_piover2) + (m - a*omega)*(swsh_piover2)) + 
            ((2*a*En - 2*Lz)*sqrt_one_minus_a_sq*(1 + sqrt_one_minus_a_sq)*(a*Lz*(1 + sqrt_one_minus_a_sq) + 2*En*(a^2 - 2*(1 + sqrt_one_minus_a_sq)))*(8*horizon_term*En*(1 + sqrt_one_minus_a_sq) + 2*sqrt_one_minus_a_sq*(((-a)*En + Lz)^2 + (1 + sqrt_one_minus_a_sq)^2))*((psptheta_piover2) + (m - a*omega)*(swsh_piover2)))/horizon_term) - 
            (1/horizon_term)*(2*I*((-I)*a*m + 2*sqrt_one_minus_a_sq + 2*I*(1 + sqrt_one_minus_a_sq)*omega)*(2*horizon_term*((-a)*En + Lz)^2*((1 + sqrt_one_minus_a_sq)*(a*Lz - 8*En*(1 + sqrt_one_minus_a_sq) + a^2*En*(4 + sqrt_one_minus_a_sq))*((-a)*m + 2*(1 + sqrt_one_minus_a_sq)*omega) + 2*((-a)*Lz*(1 + sqrt_one_minus_a_sq) + En*(4 - 2*a^2 + 4*sqrt_one_minus_a_sq))*(I*sqrt_one_minus_a_sq + (1 + sqrt_one_minus_a_sq)^2*omega))*(swsh_piover2) + 
            (2*a*En - 2*Lz)*(1 + sqrt_one_minus_a_sq)*(-2*horizon_term^2*(a*Lz - 8*En*(1 + sqrt_one_minus_a_sq) + a^2*En*(4 + sqrt_one_minus_a_sq)) + (((-a)*En + Lz)^2*sqrt_one_minus_a_sq + 4*horizon_term*En*(1 + sqrt_one_minus_a_sq) + sqrt_one_minus_a_sq*(1 + sqrt_one_minus_a_sq)^2)*((-a)*Lz*(1 + sqrt_one_minus_a_sq) + En*(4 - 2*a^2 + 4*sqrt_one_minus_a_sq)))*((psptheta_piover2) + (m - a*omega)*(swsh_piover2)))) - 
            (1/(horizon_term^2*(1 + sqrt_one_minus_a_sq)^2))*(2*((-a)*Lz*(1 + sqrt_one_minus_a_sq) + En*(4 - 2*a^2 + 4*sqrt_one_minus_a_sq))*(2*horizon_term^2 + sqrt_one_minus_a_sq*(1 + sqrt_one_minus_a_sq)*(a^2*En^2 - 2*a*En*Lz + Lz^2 + (1 + sqrt_one_minus_a_sq)^2) + horizon_term*(2*a^2*En - a*Lz + En*(sqrt_one_minus_a_sq - sqrt_one_minus_a_sq^3)))*(((-a)*En + Lz)^2*(1 + sqrt_one_minus_a_sq)^2*(a*m - 2*(1 + sqrt_one_minus_a_sq)*omega)*((-a)*m + 2*omega + 2*sqrt_one_minus_a_sq*(I + omega))*(swsh_piover2) - 
            2*horizon_term*(2*a*En - 2*Lz)*(1 + sqrt_one_minus_a_sq)^2*(a*m - 2*(1 + sqrt_one_minus_a_sq)*omega)*((psptheta_piover2) + (m - a*omega)*(swsh_piover2)) - 4*horizon_term^2*((2 - a^2 + 2*sqrt_one_minus_a_sq)*(p2sptheta2_piover2) + (1 + sqrt_one_minus_a_sq)*(a^2*omega*(2*I + omega + sqrt_one_minus_a_sq*omega)*(swsh_piover2) - 2*a*(I + omega + sqrt_one_minus_a_sq*omega)*((psptheta_piover2) + m*(swsh_piover2)) + 
            (1 + sqrt_one_minus_a_sq)*(2*m*(psptheta_piover2) + (-2 + m^2)*(swsh_piover2)))))) + 
            (1/(1 + sqrt_one_minus_a_sq)^2)*(8*(horizon_term^2*((-a)*Lz*(1 + sqrt_one_minus_a_sq) + En*(4 - 2*a^2 + 4*sqrt_one_minus_a_sq))*(2*(p2sptheta2_piover2) - a^2*(p2sptheta2_piover2) + 2*sqrt_one_minus_a_sq*(p2sptheta2_piover2) + 2*(1 + sqrt_one_minus_a_sq)^2*(m - a*omega)*(psptheta_piover2) - 2*(1 + sqrt_one_minus_a_sq)^2*(swsh_piover2) + (1 + sqrt_one_minus_a_sq)^2*(m - a*omega)^2*(swsh_piover2)) - 
            ((1 + sqrt_one_minus_a_sq)*(((-a)*En + Lz)^2*sqrt_one_minus_a_sq + 4*horizon_term*En*(1 + sqrt_one_minus_a_sq) + sqrt_one_minus_a_sq*(1 + sqrt_one_minus_a_sq)^2)*((-a)*Lz*(1 + sqrt_one_minus_a_sq) + En*(4 - 2*a^2 + 4*sqrt_one_minus_a_sq)) - horizon_term^2*(2*a*Lz*(1 + sqrt_one_minus_a_sq) + a^2*En*(6 + 5*sqrt_one_minus_a_sq + sqrt_one_minus_a_sq^2) - 4*En*(3 + 5*sqrt_one_minus_a_sq + 2*sqrt_one_minus_a_sq^2)))*
            ((2 - a^2 + 2*sqrt_one_minus_a_sq)*(p2sptheta2_piover2) + (1 + sqrt_one_minus_a_sq)*(a^2*omega*(2*I + omega + sqrt_one_minus_a_sq*omega)*(swsh_piover2) - 2*a*(I + omega + sqrt_one_minus_a_sq*omega)*((psptheta_piover2) + m*(swsh_piover2)) + (1 + sqrt_one_minus_a_sq)*(2*m*(psptheta_piover2) + (-2 + m^2)*(swsh_piover2)))))))
        end
    else
        # Throw an error, this spin weight is not supported
        throw(DomainError(s, "Currently only spin weight s of +2 is supported"))
    end
end

function _nearhorizon_ansatz_zerothorder(s::Int, l::Int, m::Int, a, omega, En, Lz; swsh_piover2=nothing, psptheta_piover2=nothing, p2sptheta2_piover2=nothing, lambda=nothing)
    if isnothing(swsh_piover2) || isnothing(psptheta_piover2) || isnothing(p2sptheta2_piover2)
        swsh_piover2 = spin_weighted_spheroidal_harmonic(s, l, m, a*omega, pi/2, 0)
        psptheta_piover2 = spin_weighted_spheroidal_harmonic(s, l, m, a*omega, pi/2, 0; theta_derivative=1)
        p2sptheta2_piover2 = spin_weighted_spheroidal_harmonic(s, l, m, a*omega, pi/2, 0; theta_derivative=2)
        lambda = Teukolsky_lambda_const(a*omega, s, l, m)
    end

    # Some useful terms to pre-compute
    sqrt_one_minus_a_sq = sqrt(1 - a^2)
    horizon_term = (a*Lz - 2*En*(1 + sqrt_one_minus_a_sq))

    #=
    We have derived/shown the explicit expression for
    s = +2 ONLY
    =#
    if s == 2
        # s = +2
        coefficient_for_scriptA0 = begin
            (-(1/(1 + sqrt_one_minus_a_sq)^5))*(6*sqrt_one_minus_a_sq*(1 + sqrt_one_minus_a_sq)^3 + 24*(1 + sqrt_one_minus_a_sq)^3*(1 + 3*sqrt_one_minus_a_sq) + (4*(1 + sqrt_one_minus_a_sq)*(-19*a^4 - 48*(1 + sqrt_one_minus_a_sq) + a^2*(67 + 43*sqrt_one_minus_a_sq)))/sqrt_one_minus_a_sq + 
            (1/(2*horizon_term^2*(sqrt_one_minus_a_sq - a^2*sqrt_one_minus_a_sq)))*(1 + sqrt_one_minus_a_sq)^4*((-a)*(a*horizon_term^2*sqrt_one_minus_a_sq - 4*(-1 + a^2)*horizon_term*(a*En + Lz*sqrt_one_minus_a_sq) - 2*a*(1 - a^2)*(((-a)*En + Lz)^2*sqrt_one_minus_a_sq + 2*horizon_term*En*(1 + sqrt_one_minus_a_sq) + sqrt_one_minus_a_sq*(1 + sqrt_one_minus_a_sq)^2))*
            m^2 + 4*(2*a^7 + 8*a*(-1 + En^2)*(1 + sqrt_one_minus_a_sq) - 8*En*Lz*(1 + sqrt_one_minus_a_sq) - 6*a^5*(2 + sqrt_one_minus_a_sq) + a^6*En*Lz*(4 + sqrt_one_minus_a_sq) + 4*a^2*En*Lz*(4 + 3*sqrt_one_minus_a_sq) - 3*a^4*En*Lz*(4 + 3*sqrt_one_minus_a_sq) + 
            a^3*(18 + (14 + Lz^2)*sqrt_one_minus_a_sq - 4*En^2*(2 + sqrt_one_minus_a_sq)))*m*omega + 4*(1 + sqrt_one_minus_a_sq)*(-2*a^6*(1 + En^2) - 8*a*En*Lz*(1 + sqrt_one_minus_a_sq) + 2*(4 + 8*En^2 + Lz^2)*(1 + sqrt_one_minus_a_sq) - 2*a^5*En*Lz*(4 + sqrt_one_minus_a_sq) + 
            2*a^3*En*Lz*(8 + 7*sqrt_one_minus_a_sq) + a^4*(Lz^2 + 6*(2 + sqrt_one_minus_a_sq) + 2*En^2*(12 + 5*sqrt_one_minus_a_sq)) - a^2*(3*Lz^2*(1 + sqrt_one_minus_a_sq) + 2*(9 + 7*sqrt_one_minus_a_sq) + En^2*(38 + 30*sqrt_one_minus_a_sq)))*omega^2) - 
            (2*(1 + sqrt_one_minus_a_sq)^3*((-a)*m + 2*(1 + sqrt_one_minus_a_sq)*omega)*((-a)*m + 2*omega + 2*sqrt_one_minus_a_sq*(-2*I + omega)))/sqrt_one_minus_a_sq + (1/(horizon_term^2*sqrt_one_minus_a_sq))*3*I*(1 + sqrt_one_minus_a_sq)^3*
            (-2*a^6*(1 + En^2)*omega - a*(1 + sqrt_one_minus_a_sq)*((4 + 40*En^2 - Lz^2)*m + 128*En*Lz*omega) + 4*(1 + sqrt_one_minus_a_sq)*(-2*En*Lz*m + 64*En^2*omega + (4 + Lz^2)*omega) + 
            a^5*((-3 - sqrt_one_minus_a_sq + En^2*(5 + sqrt_one_minus_a_sq))*m - 2*En*Lz*(5 + sqrt_one_minus_a_sq)*omega) + 2*a^4*(-12*En*Lz*m + 3*En^2*(13 + 2*sqrt_one_minus_a_sq)*omega + (9 - 5*Lz^2 + 4*sqrt_one_minus_a_sq)*omega) + 
            a^3*((7 + 5*sqrt_one_minus_a_sq - 2*Lz^2*(1 + 3*sqrt_one_minus_a_sq) + En^2*(31 + 11*sqrt_one_minus_a_sq))*m + 2*En*Lz*(65 + 33*sqrt_one_minus_a_sq)*omega) - 2*a^2*(-2*En*Lz*(9 + 8*sqrt_one_minus_a_sq)*m + 4*(4 + 3*sqrt_one_minus_a_sq)*omega - Lz^2*(4 + 5*sqrt_one_minus_a_sq)*omega + 
            2*En^2*(81 + 49*sqrt_one_minus_a_sq)*omega)) + (1/(horizon_term^2*sqrt_one_minus_a_sq))*2*I*(1 + sqrt_one_minus_a_sq)^2*(4*a^6*(5 + sqrt_one_minus_a_sq + En^2*(7 + sqrt_one_minus_a_sq))*omega + 2*a^7*((-1 + En^2)*m - 2*En*Lz*omega) + 
            4*a*(1 + sqrt_one_minus_a_sq)*((4 + 44*En^2 - Lz^2)*m + 136*En*Lz*omega) - 16*(1 + sqrt_one_minus_a_sq)*(-2*En*Lz*m + 68*En^2*omega + (4 + Lz^2)*omega) + 
            a^5*((18 - 13*Lz^2 + 8*sqrt_one_minus_a_sq - 2*En^2*(-7 + 6*sqrt_one_minus_a_sq))*m + 4*En*Lz*(41 + 6*sqrt_one_minus_a_sq)*omega) - a^3*((8*(4 + 3*sqrt_one_minus_a_sq) + 8*En^2*(23 + 12*sqrt_one_minus_a_sq) - Lz^2*(19 + 17*sqrt_one_minus_a_sq))*m + 16*En*Lz*(43 + 26*sqrt_one_minus_a_sq)*omega) + 
            8*a^2*(-2*En*Lz*(10 + 9*sqrt_one_minus_a_sq)*m - Lz^2*(4 + 5*sqrt_one_minus_a_sq)*omega + 2*(9 + 7*sqrt_one_minus_a_sq)*omega + 2*En^2*(103 + 69*sqrt_one_minus_a_sq)*omega) - 
            2*a^4*(-2*En*Lz*(30 + 13*sqrt_one_minus_a_sq)*m + 2*En^2*(151 + 47*sqrt_one_minus_a_sq)*omega + (50 + 26*sqrt_one_minus_a_sq - 11*Lz^2*(2 + sqrt_one_minus_a_sq))*omega)) + 
            ((1 + sqrt_one_minus_a_sq)^4*(-8*omega^2 + 4*a^2*omega*(I + omega) + 2*a*m*(-I + omega + sqrt_one_minus_a_sq*omega) + sqrt_one_minus_a_sq*(lambda - 8*omega^2)))/sqrt_one_minus_a_sq + ((2*(1 + sqrt_one_minus_a_sq)^3)/sqrt_one_minus_a_sq - (1 + sqrt_one_minus_a_sq)^4/(2*(-1 + a^2)))*
            (a^2*m^2 + 4*(1 + sqrt_one_minus_a_sq)*omega*(omega + sqrt_one_minus_a_sq*(-I + omega)) - 2*a*m*(2*omega + sqrt_one_minus_a_sq*(-I + 2*omega))))
        end
        
        coefficient_for_scriptA1 = begin
            -((4*I*(-1 + a^2)*(a^3*m - 2*a*(1 + sqrt_one_minus_a_sq)*m + a^2*(2*I + sqrt_one_minus_a_sq*(I - 2*omega) - 6*omega) + 2*(1 + sqrt_one_minus_a_sq)*(-I + 4*omega)))/(sqrt_one_minus_a_sq*(1 + sqrt_one_minus_a_sq)^5))
        end

        return coefficient_for_scriptA0, coefficient_for_scriptA1
    else
        # Throw an error, this spin weight is not supported
        throw(DomainError(s, "Currently only spin weight s of +2 is supported"))
    end
end

function _nearhorizon_ansatz_firstorder(s::Int, l::Int, m::Int, a, omega, En, Lz; swsh_piover2=nothing, psptheta_piover2=nothing, p2sptheta2_piover2=nothing, lambda=nothing)
    if isnothing(swsh_piover2) || isnothing(psptheta_piover2) || isnothing(p2sptheta2_piover2)
        swsh_piover2 = spin_weighted_spheroidal_harmonic(s, l, m, a*omega, pi/2, 0)
        psptheta_piover2 = spin_weighted_spheroidal_harmonic(s, l, m, a*omega, pi/2, 0; theta_derivative=1)
        p2sptheta2_piover2 = spin_weighted_spheroidal_harmonic(s, l, m, a*omega, pi/2, 0; theta_derivative=2)
        lambda = Teukolsky_lambda_const(a*omega, s, l, m)
    end

    # Some useful terms to pre-compute
    sqrt_one_minus_a_sq = sqrt(1 - a^2)
    horizon_term = (a*Lz - 2*En*(1 + sqrt_one_minus_a_sq))

    #=
    We have derived/shown the explicit expression for
    s = +2 ONLY
    =#
    if s == 2
        # s = +2
        coefficient_for_scriptA0 = begin
            (-(1/(16*horizon_term^4*(1 + sqrt_one_minus_a_sq)^5*(sqrt_one_minus_a_sq - a^2*sqrt_one_minus_a_sq)^2)))*(72*(-1 + a^2)^2*horizon_term^4*sqrt_one_minus_a_sq*(1 + sqrt_one_minus_a_sq)^2*(1 + 3*sqrt_one_minus_a_sq) + 96*(1 - a^2)*(-1 + a^2)*horizon_term^4*(1 + sqrt_one_minus_a_sq)^2*(15*a^2 - 2*(8 + 5*sqrt_one_minus_a_sq)) + 
            16*(1 - a^2)*(-1 + a^2)*horizon_term^4*(101*a^4 + 324*(1 + sqrt_one_minus_a_sq) - 6*a^2*(70 + 43*sqrt_one_minus_a_sq)) - 4*horizon_term^2*(1 + sqrt_one_minus_a_sq)^3*(5 - 5*a^2 + sqrt_one_minus_a_sq)*
            (a*(a*horizon_term^2*sqrt_one_minus_a_sq - 4*(-1 + a^2)*horizon_term*(a*En + Lz*sqrt_one_minus_a_sq) - 2*a*(1 - a^2)*(((-a)*En + Lz)^2*sqrt_one_minus_a_sq + 2*horizon_term*En*(1 + sqrt_one_minus_a_sq) + sqrt_one_minus_a_sq*(1 + sqrt_one_minus_a_sq)^2))*m^2 - 
            4*(2*a^7 + 8*a*(-1 + En^2)*(1 + sqrt_one_minus_a_sq) - 8*En*Lz*(1 + sqrt_one_minus_a_sq) - 6*a^5*(2 + sqrt_one_minus_a_sq) + a^6*En*Lz*(4 + sqrt_one_minus_a_sq) + 4*a^2*En*Lz*(4 + 3*sqrt_one_minus_a_sq) - 3*a^4*En*Lz*(4 + 3*sqrt_one_minus_a_sq) + 
            a^3*(18 + (14 + Lz^2)*sqrt_one_minus_a_sq - 4*En^2*(2 + sqrt_one_minus_a_sq)))*m*omega - 4*(1 + sqrt_one_minus_a_sq)*(-2*a^6*(1 + En^2) - 8*a*En*Lz*(1 + sqrt_one_minus_a_sq) + 2*(4 + 8*En^2 + Lz^2)*(1 + sqrt_one_minus_a_sq) - 2*a^5*En*Lz*(4 + sqrt_one_minus_a_sq) + 
            2*a^3*En*Lz*(8 + 7*sqrt_one_minus_a_sq) + a^4*(Lz^2 + 6*(2 + sqrt_one_minus_a_sq) + 2*En^2*(12 + 5*sqrt_one_minus_a_sq)) - a^2*(3*Lz^2*(1 + sqrt_one_minus_a_sq) + 2*(9 + 7*sqrt_one_minus_a_sq) + En^2*(38 + 30*sqrt_one_minus_a_sq)))*omega^2) - 
            6*I*(1 + sqrt_one_minus_a_sq)^2*(-2*(-1 + a^2)*horizon_term^4*sqrt_one_minus_a_sq*(1 + 25*(1 - a^2) + 14*sqrt_one_minus_a_sq)*(a*m - 2*(1 + sqrt_one_minus_a_sq)*omega) - horizon_term^2*(1 + sqrt_one_minus_a_sq)*(3 - a^2 + 10*sqrt_one_minus_a_sq + sqrt_one_minus_a_sq^3 - 2*a^2*(1 + 5*sqrt_one_minus_a_sq))*
            (horizon_term*(4*a*En - 4*a^3*En + a*horizon_term*sqrt_one_minus_a_sq + 4*Lz*sqrt_one_minus_a_sq - 4*a^2*Lz*sqrt_one_minus_a_sq)*m - 2*a*(1 - a^2)*(((-a)*En + Lz)^2*sqrt_one_minus_a_sq + 2*horizon_term*En*(1 + sqrt_one_minus_a_sq) + sqrt_one_minus_a_sq*(1 + sqrt_one_minus_a_sq)^2)*m + 
            4*(1 - a^2)*(1 + sqrt_one_minus_a_sq)*(((-a)*En + Lz)^2*sqrt_one_minus_a_sq + 2*horizon_term*En*(1 + sqrt_one_minus_a_sq) + sqrt_one_minus_a_sq*(1 + sqrt_one_minus_a_sq)^2)*omega + 2*horizon_term*(-4*(-1 + a^2)*horizon_term - 4*a^4*En + 3*a^3*Lz - 2*En*sqrt_one_minus_a_sq*(1 + sqrt_one_minus_a_sq)^2 - a*Lz*(3 + sqrt_one_minus_a_sq) + 
            2*a^2*En*(2 + 2*(1 - a^2) + 2*sqrt_one_minus_a_sq + sqrt_one_minus_a_sq^3))*omega) - sqrt_one_minus_a_sq*(1 + sqrt_one_minus_a_sq)*((1 + sqrt_one_minus_a_sq)*((-horizon_term^3)*sqrt_one_minus_a_sq*(4*a*En - 4*a^3*En + a*horizon_term*sqrt_one_minus_a_sq) + 
            2*(1 - a^2)*horizon_term*(4*a*En - 4*a^3*En + a*horizon_term*sqrt_one_minus_a_sq + 4*Lz*sqrt_one_minus_a_sq - 4*a^2*Lz*sqrt_one_minus_a_sq)*(((-a)*En + Lz)^2*sqrt_one_minus_a_sq + 2*horizon_term*En*(1 + sqrt_one_minus_a_sq) + sqrt_one_minus_a_sq*(1 + sqrt_one_minus_a_sq)^2) - 
            2*a*(1 - a^2)*(3*(1 - a^2)*(((-a)*En + Lz)^2*sqrt_one_minus_a_sq + 2*horizon_term*En*(1 + sqrt_one_minus_a_sq) + sqrt_one_minus_a_sq*(1 + sqrt_one_minus_a_sq)^2)^2 + horizon_term^2*sqrt_one_minus_a_sq*(4 + (1 - a^2)*(2 - 12*En^2) + (5 - 12*En^2 + Lz^2)*sqrt_one_minus_a_sq + sqrt_one_minus_a_sq^3 + 
            a^2*(-4 + (-4 + 5*En^2)*sqrt_one_minus_a_sq))))*m + 2*(2*(1 - a^2)*horizon_term*(1 + sqrt_one_minus_a_sq)*(((-a)*En + Lz)^2*sqrt_one_minus_a_sq + 2*horizon_term*En*(1 + sqrt_one_minus_a_sq) + sqrt_one_minus_a_sq*(1 + sqrt_one_minus_a_sq)^2)*
            (-4*(-1 + a^2)*horizon_term - 4*a^4*En + 3*a^3*Lz - 2*En*sqrt_one_minus_a_sq*(1 + sqrt_one_minus_a_sq)^2 - a*Lz*(3 + sqrt_one_minus_a_sq) + 2*a^2*En*(2 + 2*(1 - a^2) + 2*sqrt_one_minus_a_sq + sqrt_one_minus_a_sq^3)) + 
            2*(1 - a^2)*(1 + sqrt_one_minus_a_sq)^2*(3*(1 - a^2)*(((-a)*En + Lz)^2*sqrt_one_minus_a_sq + 2*horizon_term*En*(1 + sqrt_one_minus_a_sq) + sqrt_one_minus_a_sq*(1 + sqrt_one_minus_a_sq)^2)^2 + horizon_term^2*sqrt_one_minus_a_sq*(4 + (1 - a^2)*(2 - 12*En^2) + (5 - 12*En^2 + Lz^2)*sqrt_one_minus_a_sq + sqrt_one_minus_a_sq^3 + 
            a^2*(-4 + (-4 + 5*En^2)*sqrt_one_minus_a_sq))) + horizon_term^3*(-4*(1 - a^2)*(-1 + a^2)*horizon_term - a*(1 - a^2)*(-8*a^3*En + 7*a^2*Lz - 4*Lz*(2 + sqrt_one_minus_a_sq) + 2*a*En*(5 + 3*sqrt_one_minus_a_sq)) + 
            4*(1 - a^2)*(-4*a^4*En + 3*a^3*Lz - 2*En*sqrt_one_minus_a_sq*(1 + sqrt_one_minus_a_sq)^2 - a*Lz*(3 + sqrt_one_minus_a_sq) + 2*a^2*En*(2 + 2*(1 - a^2) + 2*sqrt_one_minus_a_sq + sqrt_one_minus_a_sq^3))))*omega)) + 
            8*I*(1 + sqrt_one_minus_a_sq)*(-2*(-1 + a^2)*horizon_term^4*(6 + 34*(1 - a^2) + 25*sqrt_one_minus_a_sq + 20*sqrt_one_minus_a_sq^3 + 5*sqrt_one_minus_a_sq^4 - 6*(a + 2*a*sqrt_one_minus_a_sq)^2)*(a*m - 2*(1 + sqrt_one_minus_a_sq)*omega) - 
            horizon_term^2*(1 + sqrt_one_minus_a_sq)*(8 + 6*a^4 + 7*(1 - a^2) + 14*sqrt_one_minus_a_sq + sqrt_one_minus_a_sq^3 - 2*a^2*(7 + 3*(1 - a^2) + 7*sqrt_one_minus_a_sq))*(horizon_term*(4*a*En - 4*a^3*En + a*horizon_term*sqrt_one_minus_a_sq + 4*Lz*sqrt_one_minus_a_sq - 4*a^2*Lz*sqrt_one_minus_a_sq)*m - 
            2*a*(1 - a^2)*(((-a)*En + Lz)^2*sqrt_one_minus_a_sq + 2*horizon_term*En*(1 + sqrt_one_minus_a_sq) + sqrt_one_minus_a_sq*(1 + sqrt_one_minus_a_sq)^2)*m + 4*(1 - a^2)*(1 + sqrt_one_minus_a_sq)*(((-a)*En + Lz)^2*sqrt_one_minus_a_sq + 2*horizon_term*En*(1 + sqrt_one_minus_a_sq) + sqrt_one_minus_a_sq*(1 + sqrt_one_minus_a_sq)^2)*omega + 
            2*horizon_term*(-4*(-1 + a^2)*horizon_term - 4*a^4*En + 3*a^3*Lz - 2*En*sqrt_one_minus_a_sq*(1 + sqrt_one_minus_a_sq)^2 - a*Lz*(3 + sqrt_one_minus_a_sq) + 2*a^2*En*(2 + 2*(1 - a^2) + 2*sqrt_one_minus_a_sq + sqrt_one_minus_a_sq^3))*omega) - 
            (1 + sqrt_one_minus_a_sq)*(1 - a^2 + sqrt_one_minus_a_sq)*((1 + sqrt_one_minus_a_sq)*((-horizon_term^3)*sqrt_one_minus_a_sq*(4*a*En - 4*a^3*En + a*horizon_term*sqrt_one_minus_a_sq) + 2*(1 - a^2)*horizon_term*(4*a*En - 4*a^3*En + a*horizon_term*sqrt_one_minus_a_sq + 4*Lz*sqrt_one_minus_a_sq - 4*a^2*Lz*sqrt_one_minus_a_sq)*
            (((-a)*En + Lz)^2*sqrt_one_minus_a_sq + 2*horizon_term*En*(1 + sqrt_one_minus_a_sq) + sqrt_one_minus_a_sq*(1 + sqrt_one_minus_a_sq)^2) - 2*a*(1 - a^2)*(3*(1 - a^2)*(((-a)*En + Lz)^2*sqrt_one_minus_a_sq + 2*horizon_term*En*(1 + sqrt_one_minus_a_sq) + sqrt_one_minus_a_sq*(1 + sqrt_one_minus_a_sq)^2)^2 + 
            horizon_term^2*sqrt_one_minus_a_sq*(4 + (1 - a^2)*(2 - 12*En^2) + (5 - 12*En^2 + Lz^2)*sqrt_one_minus_a_sq + sqrt_one_minus_a_sq^3 + a^2*(-4 + (-4 + 5*En^2)*sqrt_one_minus_a_sq))))*m + 
            2*(2*(1 - a^2)*horizon_term*(1 + sqrt_one_minus_a_sq)*(((-a)*En + Lz)^2*sqrt_one_minus_a_sq + 2*horizon_term*En*(1 + sqrt_one_minus_a_sq) + sqrt_one_minus_a_sq*(1 + sqrt_one_minus_a_sq)^2)*(-4*(-1 + a^2)*horizon_term - 4*a^4*En + 3*a^3*Lz - 2*En*sqrt_one_minus_a_sq*(1 + sqrt_one_minus_a_sq)^2 - a*Lz*(3 + sqrt_one_minus_a_sq) + 
            2*a^2*En*(2 + 2*(1 - a^2) + 2*sqrt_one_minus_a_sq + sqrt_one_minus_a_sq^3)) + 2*(1 - a^2)*(1 + sqrt_one_minus_a_sq)^2*(3*(1 - a^2)*(((-a)*En + Lz)^2*sqrt_one_minus_a_sq + 2*horizon_term*En*(1 + sqrt_one_minus_a_sq) + sqrt_one_minus_a_sq*(1 + sqrt_one_minus_a_sq)^2)^2 + 
            horizon_term^2*sqrt_one_minus_a_sq*(4 + (1 - a^2)*(2 - 12*En^2) + (5 - 12*En^2 + Lz^2)*sqrt_one_minus_a_sq + sqrt_one_minus_a_sq^3 + a^2*(-4 + (-4 + 5*En^2)*sqrt_one_minus_a_sq))) + 
            horizon_term^3*(-4*(1 - a^2)*(-1 + a^2)*horizon_term - a*(1 - a^2)*(-8*a^3*En + 7*a^2*Lz - 4*Lz*(2 + sqrt_one_minus_a_sq) + 2*a*En*(5 + 3*sqrt_one_minus_a_sq)) + 4*(1 - a^2)*(-4*a^4*En + 3*a^3*Lz - 2*En*sqrt_one_minus_a_sq*(1 + sqrt_one_minus_a_sq)^2 - 
            a*Lz*(3 + sqrt_one_minus_a_sq) + 2*a^2*En*(2 + 2*(1 - a^2) + 2*sqrt_one_minus_a_sq + sqrt_one_minus_a_sq^3))))*omega)) + 
            (1 + sqrt_one_minus_a_sq)^3*(-((1 + sqrt_one_minus_a_sq)*(8*a*(1 - a^2)*horizon_term*(((-a)*En + Lz)^2*sqrt_one_minus_a_sq + 2*horizon_term*En*(1 + sqrt_one_minus_a_sq) + sqrt_one_minus_a_sq*(1 + sqrt_one_minus_a_sq)^2)*(a*horizon_term*sqrt_one_minus_a_sq - 4*(-1 + a^2)*(a*En + Lz*sqrt_one_minus_a_sq)) - 
            horizon_term^2*sqrt_one_minus_a_sq*(3*a^2*horizon_term^2*sqrt_one_minus_a_sq - 16*a*(-1 + a^2)*horizon_term*(a*En + Lz*sqrt_one_minus_a_sq) - 8*(-1 + a^2)*sqrt_one_minus_a_sq*(2*Lz^2 + a^2*(2*En^2 - 3*Lz^2) + 2*a*En*Lz*(1 + 3*sqrt_one_minus_a_sq))) - 
            4*a^2*(1 - a^2)*(4*(1 - a^2)*(((-a)*En + Lz)^2*sqrt_one_minus_a_sq + 2*horizon_term*En*(1 + sqrt_one_minus_a_sq) + sqrt_one_minus_a_sq*(1 + sqrt_one_minus_a_sq)^2)^2 + horizon_term^2*sqrt_one_minus_a_sq*(4 + (1 - a^2)*(2 - 12*En^2) + (5 - 12*En^2 + Lz^2)*sqrt_one_minus_a_sq + sqrt_one_minus_a_sq^3 + 
            a^2*(-4 + (-4 + 5*En^2)*sqrt_one_minus_a_sq))))*m^2) - 4*(1 - a^2)*(16*a^11*(-1 + En^2) - 128*En*Lz*(-4 + 4*En^2 - Lz^2)*(1 + sqrt_one_minus_a_sq) + 128*a*(4 + 4*En^4 + Lz^2 - 4*En^2*(2 + Lz^2))*(1 + sqrt_one_minus_a_sq) + 
            8*a^10*En*Lz*(-7 - sqrt_one_minus_a_sq + En^2*(9 + sqrt_one_minus_a_sq)) + 64*a^2*En*Lz*(-34 - 30*sqrt_one_minus_a_sq - Lz^2*(5 + 4*sqrt_one_minus_a_sq) + 4*En^2*(6 + 5*sqrt_one_minus_a_sq)) - 
            64*a^3*(4*En^4*(2 + sqrt_one_minus_a_sq) + Lz^2*(5 + 4*sqrt_one_minus_a_sq) + 4*(7 + 6*sqrt_one_minus_a_sq) - 2*En^2*(2 + Lz^2)*(13 + 11*sqrt_one_minus_a_sq)) - 4*a^9*(-80 + Lz^2 - 24*sqrt_one_minus_a_sq + 2*En^2*(58 + 5*Lz^2 + 16*sqrt_one_minus_a_sq)) - 
            16*a^4*En*Lz*((-Lz^2)*(13 + 6*sqrt_one_minus_a_sq) - 8*(27 + 20*sqrt_one_minus_a_sq) + 6*En^2*(29 + 23*sqrt_one_minus_a_sq)) - 2*a^8*En*Lz*(Lz^2*(7 + 3*sqrt_one_minus_a_sq) - 16*(23 + 8*sqrt_one_minus_a_sq) + 4*En^2*(93 + 35*sqrt_one_minus_a_sq)) + 
            2*a^6*En*Lz*((-Lz^2)*(13 + sqrt_one_minus_a_sq) + 8*En^2*(146 + 85*sqrt_one_minus_a_sq) - 4*(309 + 175*sqrt_one_minus_a_sq)) + 2*a^5*(-16*En^4*(-1 + sqrt_one_minus_a_sq) + Lz^4*(1 + sqrt_one_minus_a_sq) + 4*Lz^2*(31 + 17*sqrt_one_minus_a_sq) + 16*(73 + 51*sqrt_one_minus_a_sq) - 
            32*En^2*(63 + 43*sqrt_one_minus_a_sq + 2*Lz^2*(12 + 7*sqrt_one_minus_a_sq))) + a^7*(16*En^4 + Lz^4 - 4*Lz^2*(13 + 2*sqrt_one_minus_a_sq) - 16*(85 + 44*sqrt_one_minus_a_sq) + 16*En^2*(68*(2 + sqrt_one_minus_a_sq) + Lz^2*(31 + 12*sqrt_one_minus_a_sq))))*m*omega + 
            4*(1 + sqrt_one_minus_a_sq)*(8*(1 - a^2)*horizon_term*(1 + sqrt_one_minus_a_sq)*(((-a)*En + Lz)^2*sqrt_one_minus_a_sq + 2*horizon_term*En*(1 + sqrt_one_minus_a_sq) + sqrt_one_minus_a_sq*(1 + sqrt_one_minus_a_sq)^2)*(-4*(-1 + a^2)*horizon_term - 4*a^4*En + 3*a^3*Lz - 2*En*sqrt_one_minus_a_sq*(1 + sqrt_one_minus_a_sq)^2 - 
            a*Lz*(3 + sqrt_one_minus_a_sq) + 2*a^2*En*(2 + 2*(1 - a^2) + 2*sqrt_one_minus_a_sq + sqrt_one_minus_a_sq^3)) + 4*(1 - a^2)*(1 + sqrt_one_minus_a_sq)^2*(4*(1 - a^2)*(((-a)*En + Lz)^2*sqrt_one_minus_a_sq + 2*horizon_term*En*(1 + sqrt_one_minus_a_sq) + sqrt_one_minus_a_sq*(1 + sqrt_one_minus_a_sq)^2)^2 + 
            horizon_term^2*sqrt_one_minus_a_sq*(4 + (1 - a^2)*(2 - 12*En^2) + (5 - 12*En^2 + Lz^2)*sqrt_one_minus_a_sq + sqrt_one_minus_a_sq^3 + a^2*(-4 + (-4 + 5*En^2)*sqrt_one_minus_a_sq))) + 
            horizon_term^2*(-24*(1 - a^2)*(-1 + a^2)*horizon_term^2 - 2*a*(1 - a^2)*horizon_term*(-8*a^3*En + 7*a^2*Lz - 4*Lz*(2 + sqrt_one_minus_a_sq) + 2*a*En*(5 + 3*sqrt_one_minus_a_sq)) + 
            (4*a^4*En - 3*a^3*Lz + 2*En*sqrt_one_minus_a_sq*(1 + sqrt_one_minus_a_sq)^2 + a*Lz*(3 + sqrt_one_minus_a_sq) - 2*a^2*En*(2 + 2*(1 - a^2) + 2*sqrt_one_minus_a_sq + sqrt_one_minus_a_sq^3))^2 + 
            16*(1 - a^2)*horizon_term*(-4*a^4*En + 3*a^3*Lz - 2*En*sqrt_one_minus_a_sq*(1 + sqrt_one_minus_a_sq)^2 - a*Lz*(3 + sqrt_one_minus_a_sq) + 2*a^2*En*(2 + 2*(1 - a^2) + 2*sqrt_one_minus_a_sq + sqrt_one_minus_a_sq^3))))*omega^2 - 
            2*I*((1 + sqrt_one_minus_a_sq)*(-4*(-1 + a^2)*horizon_term^3*sqrt_one_minus_a_sq*((-1 + a^2)*Lz + a*En*sqrt_one_minus_a_sq) + horizon_term^3*((-a)*(-4 + 4*a^2 + 3*(1 - a^2))*horizon_term*sqrt_one_minus_a_sq - 4*(-1 + a^2)*Lz*sqrt_one_minus_a_sq^3) - 
            8*(1 - a^2)*(-1 + a^2)*horizon_term*sqrt_one_minus_a_sq*(a*En + Lz*sqrt_one_minus_a_sq)*(((-a)*En + Lz)^2*sqrt_one_minus_a_sq + 2*horizon_term*En*(1 + sqrt_one_minus_a_sq) + sqrt_one_minus_a_sq*(1 + sqrt_one_minus_a_sq)^2) - 8*(-1 + a^2)^2*horizon_term*sqrt_one_minus_a_sq*(a*En + Lz*sqrt_one_minus_a_sq)*
            (((-a)*En + Lz)^2*sqrt_one_minus_a_sq + 2*horizon_term*En*(1 + sqrt_one_minus_a_sq) + sqrt_one_minus_a_sq*(1 + sqrt_one_minus_a_sq)^2) - 2*a*(1 - a^2)*sqrt_one_minus_a_sq*(3*(1 - a^2)*(((-a)*En + Lz)^2*sqrt_one_minus_a_sq + 2*horizon_term*En*(1 + sqrt_one_minus_a_sq) + sqrt_one_minus_a_sq*(1 + sqrt_one_minus_a_sq)^2)^2 + 
            horizon_term^2*sqrt_one_minus_a_sq*(4 + (1 - a^2)*(2 - 12*En^2) + (5 - 12*En^2 + Lz^2)*sqrt_one_minus_a_sq + sqrt_one_minus_a_sq^3 + a^2*(-4 + (-4 + 5*En^2)*sqrt_one_minus_a_sq))) + 
            (1 - a^2)*sqrt_one_minus_a_sq*(-12*a*(-1 + a^2)*(2*horizon_term*En*(1 + sqrt_one_minus_a_sq) + sqrt_one_minus_a_sq*(a^2*En^2 - 2*a*En*Lz + Lz^2 + (1 + sqrt_one_minus_a_sq)^2))^2 + 
            horizon_term*(-4*a*(-1 + a^2)*horizon_term*(6 + 5*a^2*(-1 + En^2) + Lz^2 + 6*sqrt_one_minus_a_sq - 12*En^2*(1 + sqrt_one_minus_a_sq)) - (4*a*En - 4*a^3*En + a*horizon_term*sqrt_one_minus_a_sq + 4*Lz*sqrt_one_minus_a_sq - 4*a^2*Lz*sqrt_one_minus_a_sq)*
            (4*horizon_term*En*(1 + sqrt_one_minus_a_sq) + 2*sqrt_one_minus_a_sq*(((-a)*En + Lz)^2 + (1 + sqrt_one_minus_a_sq)^2)))))*m + 
            2*(1 - a^2)*sqrt_one_minus_a_sq*(8*(1 - a^2)*horizon_term^2*(1 + sqrt_one_minus_a_sq)*(((-a)*En + Lz)^2*sqrt_one_minus_a_sq + 2*horizon_term*En*(1 + sqrt_one_minus_a_sq) + sqrt_one_minus_a_sq*(1 + sqrt_one_minus_a_sq)^2) + 
            horizon_term^3*(4*(1 - a^2)*horizon_term - 8*a^4*En + 7*a^3*Lz - 4*a*Lz*(2 + sqrt_one_minus_a_sq) + 2*a^2*En*(5 + 3*sqrt_one_minus_a_sq)) - 4*horizon_term^2*(-2*a^6*(1 + En^2) + 2*(4 + Lz^2)*(1 + sqrt_one_minus_a_sq) - 2*a^5*En*Lz*(4 + sqrt_one_minus_a_sq) + 
            2*a^3*En*Lz*(4 + 3*sqrt_one_minus_a_sq) + a^4*(3*(4 + Lz^2 + 2*sqrt_one_minus_a_sq) + 2*En^2*(8 + 5*sqrt_one_minus_a_sq)) - a^2*(14*En^2*(1 + sqrt_one_minus_a_sq) + Lz^2*(5 + 3*sqrt_one_minus_a_sq) + 2*(9 + 7*sqrt_one_minus_a_sq))) + 
            2*(1 + sqrt_one_minus_a_sq)^2*(3*(1 - a^2)*(((-a)*En + Lz)^2*sqrt_one_minus_a_sq + 2*horizon_term*En*(1 + sqrt_one_minus_a_sq) + sqrt_one_minus_a_sq*(1 + sqrt_one_minus_a_sq)^2)^2 + horizon_term^2*sqrt_one_minus_a_sq*(4 + (1 - a^2)*(2 - 12*En^2) + (5 - 12*En^2 + Lz^2)*sqrt_one_minus_a_sq + sqrt_one_minus_a_sq^3 + 
            a^2*(-4 + (-4 + 5*En^2)*sqrt_one_minus_a_sq))) - (1 + sqrt_one_minus_a_sq)*(-12*(-1 + a^2)*(1 + sqrt_one_minus_a_sq)*(2*horizon_term*En*(1 + sqrt_one_minus_a_sq) + sqrt_one_minus_a_sq*(a^2*En^2 - 2*a*En*Lz + Lz^2 + (1 + sqrt_one_minus_a_sq)^2))^2 + 
            horizon_term*(-4*(-1 + a^2)*horizon_term*(1 + sqrt_one_minus_a_sq)*(6 + 5*a^2*(-1 + En^2) + Lz^2 + 6*sqrt_one_minus_a_sq - 12*En^2*(1 + sqrt_one_minus_a_sq)) + (-4*(-1 + a^2)*horizon_term - 4*a^4*En + 3*a^3*Lz - 2*En*sqrt_one_minus_a_sq*(1 + sqrt_one_minus_a_sq)^2 - a*Lz*(3 + sqrt_one_minus_a_sq) + 
            2*a^2*En*(2 + 2*(1 - a^2) + 2*sqrt_one_minus_a_sq + sqrt_one_minus_a_sq^3))*(4*horizon_term*En*(1 + sqrt_one_minus_a_sq) + 2*sqrt_one_minus_a_sq*(((-a)*En + Lz)^2 + (1 + sqrt_one_minus_a_sq)^2)))))*omega)) + 
            horizon_term^4*sqrt_one_minus_a_sq*(1 + sqrt_one_minus_a_sq)^2*(16 + 2*(1 - a^2) + 41*sqrt_one_minus_a_sq + sqrt_one_minus_a_sq^3 - 8*a^2*(2 + 5*sqrt_one_minus_a_sq))*(a^2*m^2 + 4*(1 + sqrt_one_minus_a_sq)*omega*(omega + sqrt_one_minus_a_sq*(-I + omega)) - 2*a*m*(2*omega + sqrt_one_minus_a_sq*(-I + 2*omega))) - 
            2*horizon_term^2*(10 - 10*a^2 + 3*(1 - a^2) + 3*sqrt_one_minus_a_sq)*(-12*(1 - a^2)*(-1 + a^2)*horizon_term^2*(1 + sqrt_one_minus_a_sq)^2 - 48*(-1 + a^2)*horizon_term^2*sqrt_one_minus_a_sq*(1 + sqrt_one_minus_a_sq)^2*(1 + 3*sqrt_one_minus_a_sq) + 
            8*(-1 + a^2)*horizon_term^2*(19*a^4 + 48*(1 + sqrt_one_minus_a_sq) - a^2*(67 + 43*sqrt_one_minus_a_sq)) - 
            (1 + sqrt_one_minus_a_sq)^3*(a*(a*horizon_term^2*sqrt_one_minus_a_sq - 4*(-1 + a^2)*horizon_term*(a*En + Lz*sqrt_one_minus_a_sq) - 2*a*(1 - a^2)*(((-a)*En + Lz)^2*sqrt_one_minus_a_sq + 2*horizon_term*En*(1 + sqrt_one_minus_a_sq) + sqrt_one_minus_a_sq*(1 + sqrt_one_minus_a_sq)^2))*m^2 - 
            4*(2*a^7 + 8*a*(-1 + En^2)*(1 + sqrt_one_minus_a_sq) - 8*En*Lz*(1 + sqrt_one_minus_a_sq) - 6*a^5*(2 + sqrt_one_minus_a_sq) + a^6*En*Lz*(4 + sqrt_one_minus_a_sq) + 4*a^2*En*Lz*(4 + 3*sqrt_one_minus_a_sq) - 3*a^4*En*Lz*(4 + 3*sqrt_one_minus_a_sq) + 
            a^3*(18 + (14 + Lz^2)*sqrt_one_minus_a_sq - 4*En^2*(2 + sqrt_one_minus_a_sq)))*m*omega - 4*(1 + sqrt_one_minus_a_sq)*(-2*a^6*(1 + En^2) - 8*a*En*Lz*(1 + sqrt_one_minus_a_sq) + 2*(4 + 8*En^2 + Lz^2)*(1 + sqrt_one_minus_a_sq) - 2*a^5*En*Lz*(4 + sqrt_one_minus_a_sq) + 
            2*a^3*En*Lz*(8 + 7*sqrt_one_minus_a_sq) + a^4*(Lz^2 + 6*(2 + sqrt_one_minus_a_sq) + 2*En^2*(12 + 5*sqrt_one_minus_a_sq)) - a^2*(3*Lz^2*(1 + sqrt_one_minus_a_sq) + 2*(9 + 7*sqrt_one_minus_a_sq) + En^2*(38 + 30*sqrt_one_minus_a_sq)))*omega^2) + 
            4*(-1 + a^2)*horizon_term^2*(1 + sqrt_one_minus_a_sq)^2*(a*m - 2*(1 + sqrt_one_minus_a_sq)*omega)*(a*m - 2*omega - 2*sqrt_one_minus_a_sq*(-2*I + omega)) + 6*I*(-1 + a^2)*(1 + sqrt_one_minus_a_sq)^2*
            (2*a^6*(1 + En^2)*omega + a*(1 + sqrt_one_minus_a_sq)*((4 + 40*En^2 - Lz^2)*m + 128*En*Lz*omega) - 4*(1 + sqrt_one_minus_a_sq)*(-2*En*Lz*m + 64*En^2*omega + (4 + Lz^2)*omega) + 
            a^5*((3 + sqrt_one_minus_a_sq - En^2*(5 + sqrt_one_minus_a_sq))*m + 2*En*Lz*(5 + sqrt_one_minus_a_sq)*omega) - 2*a^4*(-12*En*Lz*m + 3*En^2*(13 + 2*sqrt_one_minus_a_sq)*omega + (9 - 5*Lz^2 + 4*sqrt_one_minus_a_sq)*omega) - 
            a^3*((7 + 5*sqrt_one_minus_a_sq - 2*Lz^2*(1 + 3*sqrt_one_minus_a_sq) + En^2*(31 + 11*sqrt_one_minus_a_sq))*m + 2*En*Lz*(65 + 33*sqrt_one_minus_a_sq)*omega) + 2*a^2*(-2*En*Lz*(9 + 8*sqrt_one_minus_a_sq)*m + 4*(4 + 3*sqrt_one_minus_a_sq)*omega - Lz^2*(4 + 5*sqrt_one_minus_a_sq)*omega + 
            2*En^2*(81 + 49*sqrt_one_minus_a_sq)*omega)) - 4*I*(-1 + a^2)*(1 + sqrt_one_minus_a_sq)*(4*a^6*(5 + sqrt_one_minus_a_sq + En^2*(7 + sqrt_one_minus_a_sq))*omega + 2*a^7*((-1 + En^2)*m - 2*En*Lz*omega) + 
            4*a*(1 + sqrt_one_minus_a_sq)*((4 + 44*En^2 - Lz^2)*m + 136*En*Lz*omega) - 16*(1 + sqrt_one_minus_a_sq)*(-2*En*Lz*m + 68*En^2*omega + (4 + Lz^2)*omega) + 
            a^5*((18 - 13*Lz^2 + 8*sqrt_one_minus_a_sq - 2*En^2*(-7 + 6*sqrt_one_minus_a_sq))*m + 4*En*Lz*(41 + 6*sqrt_one_minus_a_sq)*omega) - a^3*((8*(4 + 3*sqrt_one_minus_a_sq) + 8*En^2*(23 + 12*sqrt_one_minus_a_sq) - Lz^2*(19 + 17*sqrt_one_minus_a_sq))*m + 16*En*Lz*(43 + 26*sqrt_one_minus_a_sq)*omega) + 
            8*a^2*(-2*En*Lz*(10 + 9*sqrt_one_minus_a_sq)*m - Lz^2*(4 + 5*sqrt_one_minus_a_sq)*omega + 2*(9 + 7*sqrt_one_minus_a_sq)*omega + 2*En^2*(103 + 69*sqrt_one_minus_a_sq)*omega) - 
            2*a^4*(-2*En*Lz*(30 + 13*sqrt_one_minus_a_sq)*m + 2*En^2*(151 + 47*sqrt_one_minus_a_sq)*omega + (50 + 26*sqrt_one_minus_a_sq - 11*Lz^2*(2 + sqrt_one_minus_a_sq))*omega)) - 
            2*(-1 + a^2)*horizon_term^2*(1 + sqrt_one_minus_a_sq)^3*(-8*omega^2 + 4*a^2*omega*(I + omega) + 2*a*m*(-I + omega + sqrt_one_minus_a_sq*omega) + sqrt_one_minus_a_sq*(lambda - 8*omega^2)) + horizon_term^2*(1 + sqrt_one_minus_a_sq)^2*(5 - 5*a^2 + sqrt_one_minus_a_sq)*
            (a^2*m^2 + 4*(1 + sqrt_one_minus_a_sq)*omega*(omega + sqrt_one_minus_a_sq*(-I + omega)) - 2*a*m*(2*omega + sqrt_one_minus_a_sq*(-I + 2*omega)))) + 4*(-1 + a^2)*horizon_term^4*sqrt_one_minus_a_sq*(1 + sqrt_one_minus_a_sq)^2*
            (4*a*m*(10*I + sqrt_one_minus_a_sq*(4*I - 15*omega) - 15*omega) + 4*a^4*omega*(11*I + (13 + sqrt_one_minus_a_sq)*omega) + 2*a^3*m*(-20*I + (30 + 9*sqrt_one_minus_a_sq)*omega) - 2*(1 + sqrt_one_minus_a_sq)*(5*lambda + 4*(3*I - 28*omega)*omega) + 
            a^2*(10*lambda - 4*omega*(5*I + 69*omega) + sqrt_one_minus_a_sq*(6*m^2 + 9*lambda - 4*omega*(2*I + 41*omega)))))
        end
        
        coefficient_for_scriptA1 = begin
            (1/(4*(1 - a^2)*horizon_term^2*sqrt_one_minus_a_sq^2*(1 + sqrt_one_minus_a_sq)^6))*(4*I*(-1 + a^2)*horizon_term^2*(10 - 10*a^2 + 3*(1 - a^2) + 3*sqrt_one_minus_a_sq)*(a^3*m - 2*a*(1 + sqrt_one_minus_a_sq)*m + a^2*(2*I + sqrt_one_minus_a_sq*(I - 2*omega) - 6*omega) + 
            2*(1 + sqrt_one_minus_a_sq)*(-I + 4*omega)) + sqrt_one_minus_a_sq*(1 + sqrt_one_minus_a_sq)*(-72*(1 - a^2)*(-1 + a^2)*horizon_term^2*(1 + sqrt_one_minus_a_sq) - 12*(-1 + a^2)*horizon_term^2*sqrt_one_minus_a_sq*(1 + sqrt_one_minus_a_sq)*(5 + 9*sqrt_one_minus_a_sq) - 
            16*(-1 + a^2)*horizon_term^2*(-9*(1 + sqrt_one_minus_a_sq) + a^2*(9 + 7*sqrt_one_minus_a_sq)) - 
            2*(1 + sqrt_one_minus_a_sq)^2*(a*(a*horizon_term^2*sqrt_one_minus_a_sq - 4*(-1 + a^2)*horizon_term*(a*En + Lz*sqrt_one_minus_a_sq) - 2*a*(1 - a^2)*(((-a)*En + Lz)^2*sqrt_one_minus_a_sq + 2*horizon_term*En*(1 + sqrt_one_minus_a_sq) + sqrt_one_minus_a_sq*(1 + sqrt_one_minus_a_sq)^2))*m^2 - 
            4*(2*a^7 + 8*a*(-1 + En^2)*(1 + sqrt_one_minus_a_sq) - 8*En*Lz*(1 + sqrt_one_minus_a_sq) - 6*a^5*(2 + sqrt_one_minus_a_sq) + a^6*En*Lz*(4 + sqrt_one_minus_a_sq) + 4*a^2*En*Lz*(4 + 3*sqrt_one_minus_a_sq) - 3*a^4*En*Lz*(4 + 3*sqrt_one_minus_a_sq) + 
            a^3*(18 + (14 + Lz^2)*sqrt_one_minus_a_sq - 4*En^2*(2 + sqrt_one_minus_a_sq)))*m*omega - 4*(1 + sqrt_one_minus_a_sq)*(-2*a^6*(1 + En^2) - 8*a*En*Lz*(1 + sqrt_one_minus_a_sq) + 2*(4 + 8*En^2 + Lz^2)*(1 + sqrt_one_minus_a_sq) - 2*a^5*En*Lz*(4 + sqrt_one_minus_a_sq) + 
            2*a^3*En*Lz*(8 + 7*sqrt_one_minus_a_sq) + a^4*(Lz^2 + 6*(2 + sqrt_one_minus_a_sq) + 2*En^2*(12 + 5*sqrt_one_minus_a_sq)) - a^2*(3*Lz^2*(1 + sqrt_one_minus_a_sq) + 2*(9 + 7*sqrt_one_minus_a_sq) + En^2*(38 + 30*sqrt_one_minus_a_sq)))*omega^2) - 
            horizon_term^2*(1 + sqrt_one_minus_a_sq)*(5 - 5*a^2 + sqrt_one_minus_a_sq)*(a*m - 2*(1 + sqrt_one_minus_a_sq)*omega)*(a*m - 2*omega - 2*sqrt_one_minus_a_sq*(-2*I + omega)) + 
            6*I*(-1 + a^2)*(1 + sqrt_one_minus_a_sq)*(4*a^6*(1 + En^2)*omega + 2*a*(1 + sqrt_one_minus_a_sq)*((4 + 32*En^2 - Lz^2)*m + 112*En*Lz*omega) - 8*(1 + sqrt_one_minus_a_sq)*(-2*En*Lz*m + 56*En^2*omega + (4 + Lz^2)*omega) + 
            a^5*(6*m + 2*sqrt_one_minus_a_sq*m - 2*En^2*(5 + sqrt_one_minus_a_sq)*m + 4*En*Lz*(5 + sqrt_one_minus_a_sq)*omega) - 2*a^4*(-18*En*Lz*m + 6*En^2*(11 + 2*sqrt_one_minus_a_sq)*omega + (18 - 7*Lz^2 + 8*sqrt_one_minus_a_sq)*omega) - 
            a^3*((2*(7 + 5*sqrt_one_minus_a_sq) + 2*En^2*(21 + 5*sqrt_one_minus_a_sq) - Lz^2*(5 + 9*sqrt_one_minus_a_sq))*m + 4*En*Lz*(55 + 27*sqrt_one_minus_a_sq)*omega) + 4*a^2*(-2*En*Lz*(8 + 7*sqrt_one_minus_a_sq)*m + 4*(4 + 3*sqrt_one_minus_a_sq)*omega - Lz^2*(3 + 4*sqrt_one_minus_a_sq)*omega + 
            2*En^2*(69 + 41*sqrt_one_minus_a_sq)*omega)) + 4*I*sqrt_one_minus_a_sq*(1 + sqrt_one_minus_a_sq)*(4*a^6*(4 + sqrt_one_minus_a_sq + En^2*(6 + sqrt_one_minus_a_sq))*omega + 2*a^7*((-1 + En^2)*m - 2*En*Lz*omega) + 2*a*(1 + sqrt_one_minus_a_sq)*((4 + 56*En^2 - Lz^2)*m + 160*En*Lz*omega) - 
            8*(1 + sqrt_one_minus_a_sq)*(-2*En*Lz*m + 80*En^2*omega + (4 + Lz^2)*omega) + a^5*((12 - 15*Lz^2 + En^2*(32 - 10*sqrt_one_minus_a_sq) + 6*sqrt_one_minus_a_sq)*m + 20*En*Lz*(8 + sqrt_one_minus_a_sq)*omega) + 
            a^4*(20*En*Lz*(4 + 3*sqrt_one_minus_a_sq)*m - 4*(16 + 9*sqrt_one_minus_a_sq)*omega + 2*Lz^2*(14 + 13*sqrt_one_minus_a_sq)*omega - 4*En^2*(124 + 45*sqrt_one_minus_a_sq)*omega) - 
            a^3*(((-Lz^2)*(17 + 5*sqrt_one_minus_a_sq) + 2*(9 + 7*sqrt_one_minus_a_sq) + 2*En^2*(73 + 45*sqrt_one_minus_a_sq))*m + 4*En*Lz*(119 + 79*sqrt_one_minus_a_sq)*omega) + 4*a^2*(-2*En*Lz*(12 + 11*sqrt_one_minus_a_sq)*m + 4*(5 + 4*sqrt_one_minus_a_sq)*omega - Lz^2*(5 + 6*sqrt_one_minus_a_sq)*omega + 
            2*En^2*(139 + 99*sqrt_one_minus_a_sq)*omega)) - 4*(-1 + a^2)*horizon_term^2*(1 + sqrt_one_minus_a_sq)^2*(-8*omega^2 + 4*a^2*omega*(I + omega) + 2*a*m*(-I + omega + sqrt_one_minus_a_sq*omega) + sqrt_one_minus_a_sq*(lambda - 8*omega^2)) + 
            horizon_term^2*(1 + sqrt_one_minus_a_sq)*(4 - 4*a^2 + 3*(1 - a^2) + 3*sqrt_one_minus_a_sq)*(a^2*m^2 + 4*(1 + sqrt_one_minus_a_sq)*omega*(omega + sqrt_one_minus_a_sq*(-I + omega)) - 2*a*m*(2*omega + sqrt_one_minus_a_sq*(-I + 2*omega)))))
        end

        return coefficient_for_scriptA0, coefficient_for_scriptA1
    else
        # Throw an error, this spin weight is not supported
        throw(DomainError(s, "Currently only spin weight s of +2 is supported"))
    end
end

function sourceterm_regularization_ansatz_coefficients(s::Int, l::Int, m::Int, a, omega, En, Lz)
    # Compute the necessary angular terms using SpinWeightedSpheroidalHarmonics.jl
    # There is a caching mechanism so only need to do spectral decomposition once
    swsh_piover2 = spin_weighted_spheroidal_harmonic(s, l, m, a*omega, pi/2, 0)
    psptheta_piover2 = spin_weighted_spheroidal_harmonic(s, l, m, a*omega, pi/2, 0; theta_derivative=1)
    p2sptheta2_piover2 = spin_weighted_spheroidal_harmonic(s, l, m, a*omega, pi/2, 0; theta_derivative=2)
    lambda = Teukolsky_lambda_const(a*omega, s, l, m)

    # These alpha's are the series expansion coefficients of the RHS/source term
    alpha0 = _sourceterm_without_phasing_circularorbit_nearhorizon_seriesexpansion_zerothorder(s, l, m, a, omega, En, Lz; swsh_piover2=swsh_piover2, psptheta_piover2=psptheta_piover2, p2sptheta2_piover2=p2sptheta2_piover2)
    alpha1 = _sourceterm_without_phasing_circularorbit_nearhorizon_seriesexpansion_firstorder(s, l, m, a, omega, En, Lz; swsh_piover2=swsh_piover2, psptheta_piover2=psptheta_piover2, p2sptheta2_piover2=p2sptheta2_piover2)
    
    # These beta's are the series expansion coefficients of the LHS when the ansatz is substituted
    beta00, beta01 = _nearhorizon_ansatz_zerothorder(s, l, m, a, omega, En, Lz; swsh_piover2=swsh_piover2, psptheta_piover2=psptheta_piover2, p2sptheta2_piover2=p2sptheta2_piover2, lambda=lambda)
    beta10, beta11 = _nearhorizon_ansatz_firstorder(s, l, m, a, omega, En, Lz; swsh_piover2=swsh_piover2, psptheta_piover2=psptheta_piover2, p2sptheta2_piover2=p2sptheta2_piover2, lambda=lambda)

    # Solve the simple linear algebra
    scriptA0 = (beta11*alpha0 - beta01*alpha1)/(beta00*beta11 - beta10*beta01)
    scriptA1 = (beta10*alpha0 - beta00*alpha1)/(beta10*beta01 - beta00*beta11)

    return scriptA0, scriptA1
end

function _lhs_without_phasing_ansatz(s::Int, l::Int, m::Int, a, omega, En, Lz, r, scriptA0, scriptA1)
    lambda = Teukolsky_lambda_const(a*omega, s, l, m)

    coefficient_for_scriptA0 = begin
        (-(1/(r^3*(a^2 + (-2 + r)*r)^3)))*(24*(-1 + r)^2*r^2 + 6*(-1 + r)*r*(a^2 + (-2 + r)*r) - 2*(a^4 + 4*a^2*(-2 + r)*r + 3*r^2*(8 - 12*r + 5*r^2)) + 
        (6*I*(-1 + r)*r^3*((-Lz)*(-2 + r)*m + En*r^3*omega + a^2*En*(2 + r)*omega - 2*a*(En*m + Lz*omega)))/sqrt(-((a^2 + (-2 + r)*r)*(((-a)*En + Lz)^2 + r^2)) + (a*Lz - En*(a^2 + r^2))^2) - 
        (2*I*r^2*(a^2 + r*(-6 + 5*r))*((-Lz)*(-2 + r)*m + En*r^3*omega + a^2*En*(2 + r)*omega - 2*a*(En*m + Lz*omega)))/sqrt(-((a^2 + (-2 + r)*r)*(((-a)*En + Lz)^2 + r^2)) + (a*Lz - En*(a^2 + r^2))^2) + 
        r^2*((a^2 + (-2 + r)*r)*lambda - 8*I*r*(a^2 + (-2 + r)*r)*omega - (a*m - (a^2 + r^2)*omega)^2 + 4*I*(-1 + r)*((-a)*m + (a^2 + r^2)*omega)) + 
        (1/((a^2 + (-2 + r)*r)*(((-a)*En + Lz)^2 + r^2) - (a*Lz - En*(a^2 + r^2))^2))*r^3*((-(2*a*En + Lz*(-2 + r))^2)*r*m^2 + 2*(2*a*En + Lz*(-2 + r))*r*(-2*a*Lz + En*r^3 + a^2*En*(2 + r))*m*omega - 
        r*(-2*a*Lz + En*r^3 + a^2*En*(2 + r))^2*omega^2 + (1/sqrt(-((a^2 + (-2 + r)*r)*(((-a)*En + Lz)^2 + r^2)) + (a*Lz - En*(a^2 + r^2))^2))*I*
        (Lz^3*(-2 + r)^2*(-1 + r)*r*m - Lz*(-2 + r)^2*r^3*(3 + 2*(-1 + En^2)*r)*m + En*r^6*(6 + (-3 + 2*En^2)*r)*omega + En*Lz^2*r^4*(10 - 7*r + r^2)*omega - a^6*En*(-r^2 + En^2*(2 + 3*r + r^2))*omega + 
        2*a^5*En^2*(En*m + 3*Lz*(1 + r)*omega) - 2*a^3*(3*En*(-1 + r)*(Lz^2 - r^2)*m + En^3*r*(-2 + r + 3*r^2)*m - 3*En^2*Lz*r*(2 - 3*r + r^2)*omega - Lz*(Lz^2 + 3*(-1 + r)*r^2)*omega) + 
        a^2*((-Lz^3)*(2 - 3*r + r^2)*m + Lz*r*(r*(6 - 5*r + r^2) - En^2*(12 - 12*r - 5*r^2 + r^3))*m + 2*En*Lz^2*r*(-6 + 12*r - 5*r^2 + r^3)*omega + 
        En*r^3*(-12 + 2*(8 + En^2)*r + 5*(-2 + En^2)*r^2 - (-1 + En^2)*r^3)*omega) - 2*a*r*(En^3*r^3*(-4 + 3*r)*m - 3*En*(Lz^2 + r^2)*(2 - 3*r + r^2)*m + 6*En^2*Lz*r^3*omega - 
        Lz*(-2 + r)*(3*(-1 + r)*r^2 + Lz^2*(-1 + 2*r))*omega) + a^4*((-Lz)*r^2*m + En^2*Lz*(-6 + 3*r + r^2)*m - 2*En^3*r*(2 - 2*r + r^3)*omega + 
        En*(Lz^2*(-6 - 3*r + r^2) + r^2*(6 - 7*r + 2*r^2))*omega))))
    end

    coefficient_for_scriptA1 = begin
        (1/(2*r^5*(a^2 + (-2 + r)*r)^2))*(24*(-1 + r)^2*r^2 + 36*(-1 + r)*r*(a^2 + (-2 + r)*r) - 4*(6*a^4 + a^2*r*(-30 + 17*r) + r^2*(40 - 48*r + 15*r^2)) + 
        (2*(2*a*En + Lz*(-2 + r))^2*r^4*m^2)/(-((a^2 + (-2 + r)*r)*(((-a)*En + Lz)^2 + r^2)) + (a*Lz - En*(a^2 + r^2))^2) - (4*(2*a*En + Lz*(-2 + r))*r^4*(-2*a*Lz + En*r^3 + a^2*En*(2 + r))*m*omega)/
        (-((a^2 + (-2 + r)*r)*(((-a)*En + Lz)^2 + r^2)) + (a*Lz - En*(a^2 + r^2))^2) + (2*r^4*(-2*a*Lz + En*r^3 + a^2*En*(2 + r))^2*omega^2)/
        (-((a^2 + (-2 + r)*r)*(((-a)*En + Lz)^2 + r^2)) + (a*Lz - En*(a^2 + r^2))^2) + (1/(-((a^2 + (-2 + r)*r)*(((-a)*En + Lz)^2 + r^2)) + (a*Lz - En*(a^2 + r^2))^2)^(3/2))*I*r^2*
        ((2*(2*a*En + Lz*(-2 + r))*r*(a^2 + (-2 + r)*r)*(-2*a*En*Lz - Lz^2*(-1 + r) + r^2*(3 + 2*(-1 + En^2)*r) + a^2*(-r + En^2*(1 + r))) + 
        2*(2*a*En + Lz*(-2 + r))*(2 - 2*r)*r*((a^2 + (-2 + r)*r)*(((-a)*En + Lz)^2 + r^2) - (a*Lz - En*(a^2 + r^2))^2) + 4*(a*En + Lz*(-1 + r))*(a^2 + (-2 + r)*r)*
        ((a^2 + (-2 + r)*r)*(((-a)*En + Lz)^2 + r^2) - (a*Lz - En*(a^2 + r^2))^2))*m + 
        (-2*r*(a^2 + (-2 + r)*r)*(-2*a*Lz + En*r^3 + a^2*En*(2 + r))*(-2*a*En*Lz - Lz^2*(-1 + r) + r^2*(3 + 2*(-1 + En^2)*r) + a^2*(-r + En^2*(1 + r))) - 
        4*(a^2 + (-2 + r)*r)*(-2*a*Lz + En*r^3 + a^2*En*(2 + r))*((a^2 + (-2 + r)*r)*(((-a)*En + Lz)^2 + r^2) - (a*Lz - En*(a^2 + r^2))^2) - 
        4*((-a^4)*En + a^3*Lz - 2*a^2*En*(-2 + r)*r - En*r^4 + a*Lz*r*(-4 + 3*r))*((a^2 + (-2 + r)*r)*(((-a)*En + Lz)^2 + r^2) - (a*Lz - En*(a^2 + r^2))^2))*omega) + 
        (12*I*(-1 + r)*r^3*((-Lz)*(-2 + r)*m + En*r^3*omega + a^2*En*(2 + r)*omega - 2*a*(En*m + Lz*omega)))/sqrt(-((a^2 + (-2 + r)*r)*(((-a)*En + Lz)^2 + r^2)) + (a*Lz - En*(a^2 + r^2))^2) - 
        (4*I*r^2*(3*a^2 + r*(-8 + 5*r))*((-Lz)*(-2 + r)*m + En*r^3*omega + a^2*En*(2 + r)*omega - 2*a*(En*m + Lz*omega)))/sqrt(-((a^2 + (-2 + r)*r)*(((-a)*En + Lz)^2 + r^2)) + (a*Lz - En*(a^2 + r^2))^2) + 
        2*r^2*((a^2 + (-2 + r)*r)*lambda - 8*I*r*(a^2 + (-2 + r)*r)*omega - (a*m - (a^2 + r^2)*omega)^2 + 4*I*(-1 + r)*((-a)*m + (a^2 + r^2)*omega)))
    end

    return coefficient_for_scriptA0*scriptA0 + coefficient_for_scriptA1+scriptA1
end