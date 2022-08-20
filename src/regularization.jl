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
            -((4*I*(-2*sqrt_one_minus_a_sq^4*(3*I + a*m - 6*omega) - (-1 + a^2)*(a*m - 2*omega) + 2*(-1 + a^2)*sqrt_one_minus_a_sq*(-2*I + I*a^2 - a*m + 3*omega) + sqrt_one_minus_a_sq^5*(-3*I + 4*omega) + sqrt_one_minus_a_sq^3*(-3*I - 4*a*m + 10*omega + 2*a^2*omega) - sqrt_one_minus_a_sq^2*(a*m + a^3*m + a^2*(4*I - 6*omega) + 2*(-2*I + omega))))/(sqrt_one_minus_a_sq*(1 + sqrt_one_minus_a_sq)^5))
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
    neg_horizon_term = -horizon_term
    am_term = a*m - 2*omega - 2*sqrt_one_minus_a_sq*omega
    long_term = 2*a^4*(1 + En^2)*omega - 2*a^2*(5 + 3*sqrt_one_minus_a_sq + En^2*(13 + 5*sqrt_one_minus_a_sq))*omega + 2*(1 + sqrt_one_minus_a_sq)*(-2*En*Lz*sqrt_one_minus_a_sq*m + 16*En^2*omega + (4 + Lz^2)*omega) + a*((-2 + 4*En^2 + Lz^2)*sqrt_one_minus_a_sq*m + 2*(-1 + 2*En^2)*sqrt_one_minus_a_sq^2*m - 16*En*Lz*omega - 16*En*Lz*sqrt_one_minus_a_sq*omega)

    #=
    We have derived/shown the explicit expression for
    s = +2 ONLY
    =#
    if s == 2
        # s = +2
        coefficient_for_scriptA0 = begin
            (1/(8*(1 + sqrt_one_minus_a_sq)^5))*(72*(1 + sqrt_one_minus_a_sq)^2 + (12*I*am_term*(1 + sqrt_one_minus_a_sq)^2)/(-1 + a^2) - (20*am_term*(am_term + 4*I*sqrt_one_minus_a_sq)*(1 + sqrt_one_minus_a_sq)^2)/sqrt_one_minus_a_sq^2 - (120*I*am_term*(1 + sqrt_one_minus_a_sq)^2)/sqrt_one_minus_a_sq + (2*am_term*(am_term + 4*I*sqrt_one_minus_a_sq)*(1 + sqrt_one_minus_a_sq)^2)/((-1 + a^2)*sqrt_one_minus_a_sq) + (24*(1 + sqrt_one_minus_a_sq)^2*(7 - 8*a^2 + sqrt_one_minus_a_sq))/(-1 + a^2) - 
            (36*I*am_term*(1 + sqrt_one_minus_a_sq)^2*(1 + 2*sqrt_one_minus_a_sq))/sqrt_one_minus_a_sq^2 + (216*(1 + sqrt_one_minus_a_sq)^2*(1 + 2*sqrt_one_minus_a_sq))/sqrt_one_minus_a_sq + (240*(1 + sqrt_one_minus_a_sq)^2*(1 + 4*sqrt_one_minus_a_sq))/sqrt_one_minus_a_sq - (48*(1 + sqrt_one_minus_a_sq)^2*(1 + 2*sqrt_one_minus_a_sq)*(1 + 4*sqrt_one_minus_a_sq))/sqrt_one_minus_a_sq^2 - (30*I*am_term*(1 + sqrt_one_minus_a_sq)^2*(1 + 7*sqrt_one_minus_a_sq))/(-1 + a^2) - 
            (5*am_term*(am_term + 4*I*sqrt_one_minus_a_sq)*(1 + sqrt_one_minus_a_sq)^2*(1 + 7*sqrt_one_minus_a_sq))/((-1 + a^2)*sqrt_one_minus_a_sq) + (120*sqrt_one_minus_a_sq*(1 + sqrt_one_minus_a_sq)^2*(1 + 7*sqrt_one_minus_a_sq))/(-1 + a^2) + (40*I*am_term*(1 + sqrt_one_minus_a_sq)*(1 - a^2 + sqrt_one_minus_a_sq)*(1 + 7*sqrt_one_minus_a_sq))/((-1 + a^2)*sqrt_one_minus_a_sq) + (120*(2 - a^2 + 2*sqrt_one_minus_a_sq)*(1 + 7*sqrt_one_minus_a_sq))/sqrt_one_minus_a_sq + 
            (24*(1 + sqrt_one_minus_a_sq)^3*(1 + sqrt_one_minus_a_sq + sqrt_one_minus_a_sq^3/(-1 + a^2)))/sqrt_one_minus_a_sq^2 + (40*I*am_term*(1 + sqrt_one_minus_a_sq)*(-9*a^2 + 11*(1 + sqrt_one_minus_a_sq)))/sqrt_one_minus_a_sq^2 + 6*I*am_term*(1 + sqrt_one_minus_a_sq)^3*(1/(1 - a^2) - (1 + sqrt_one_minus_a_sq)/sqrt_one_minus_a_sq^3) - (8*I*am_term*(1 + sqrt_one_minus_a_sq)*(-10*(1 + sqrt_one_minus_a_sq) + a^2*(11 + 10*sqrt_one_minus_a_sq)))/((-1 + a^2)*sqrt_one_minus_a_sq) + 
            (1824*(1 + sqrt_one_minus_a_sq) - 16*a^2*(98 + 41*sqrt_one_minus_a_sq))/sqrt_one_minus_a_sq - (80*(19*a^4 + 48*(1 + sqrt_one_minus_a_sq) - a^2*(67 + 43*sqrt_one_minus_a_sq)))/sqrt_one_minus_a_sq^2 + (20*(1 + sqrt_one_minus_a_sq)^2*(2*I*am_term*sqrt_one_minus_a_sq + (a*m - 2*(1 + sqrt_one_minus_a_sq)*omega)^2))/sqrt_one_minus_a_sq^2 - (2*(1 + sqrt_one_minus_a_sq)^2*(2*I*am_term*sqrt_one_minus_a_sq + (a*m - 2*(1 + sqrt_one_minus_a_sq)*omega)^2))/((-1 + a^2)*sqrt_one_minus_a_sq) + 
            (5*(1 + sqrt_one_minus_a_sq)^2*(1 + 7*sqrt_one_minus_a_sq)*(2*I*am_term*sqrt_one_minus_a_sq + (a*m - 2*(1 + sqrt_one_minus_a_sq)*omega)^2))/((-1 + a^2)*sqrt_one_minus_a_sq) - 
            (1/((1 - a^2)*neg_horizon_term^4*sqrt_one_minus_a_sq^2))*(6*I*(1 + sqrt_one_minus_a_sq)^3*(-((-1 + a^2)*am_term*neg_horizon_term^4) + m*(1 + sqrt_one_minus_a_sq)*(-2*a*(-1 + a^2)*En*neg_horizon_term^3 - 4*(1 - a^2)*neg_horizon_term*sqrt_one_minus_a_sq*(a*En + Lz*sqrt_one_minus_a_sq)*(((-a)*En + Lz)^2*sqrt_one_minus_a_sq - 2*En*neg_horizon_term*(1 + sqrt_one_minus_a_sq) + sqrt_one_minus_a_sq*(1 + sqrt_one_minus_a_sq)^2) - 
            a*sqrt_one_minus_a_sq*(3*(1 - a^2)*(((-a)*En + Lz)^2*sqrt_one_minus_a_sq - 2*En*neg_horizon_term*(1 + sqrt_one_minus_a_sq) + sqrt_one_minus_a_sq*(1 + sqrt_one_minus_a_sq)^2)^2 + 2*neg_horizon_term^2*sqrt_one_minus_a_sq*(En*neg_horizon_term - 2*(-1 + a^2)*(1 + sqrt_one_minus_a_sq) - 2*En^2*sqrt_one_minus_a_sq*(1 + sqrt_one_minus_a_sq)^2))) + 
            2*sqrt_one_minus_a_sq*((-neg_horizon_term^3)*(4*a^3*Lz + 2*neg_horizon_term*sqrt_one_minus_a_sq + a*Lz*(-4 + sqrt_one_minus_a_sq + 3*sqrt_one_minus_a_sq^2) - 4*(-1 + a^2)*En*(-2 + 2*a^2 - 2*sqrt_one_minus_a_sq - sqrt_one_minus_a_sq^2 - sqrt_one_minus_a_sq^3)) + 2*(1 - a^2)*neg_horizon_term*(1 + sqrt_one_minus_a_sq)*(((-a)*En + Lz)^2*sqrt_one_minus_a_sq - 2*En*neg_horizon_term*(1 + sqrt_one_minus_a_sq) + sqrt_one_minus_a_sq*(1 + sqrt_one_minus_a_sq)^2)*
            (a*Lz + 2*neg_horizon_term + En*(2 - 2*a^2 + 3*sqrt_one_minus_a_sq + 2*sqrt_one_minus_a_sq^2 + sqrt_one_minus_a_sq^3)) + (1 + sqrt_one_minus_a_sq)^2*(3*(1 - a^2)*(((-a)*En + Lz)^2*sqrt_one_minus_a_sq - 2*En*neg_horizon_term*(1 + sqrt_one_minus_a_sq) + sqrt_one_minus_a_sq*(1 + sqrt_one_minus_a_sq)^2)^2 + 2*neg_horizon_term^2*sqrt_one_minus_a_sq*(En*neg_horizon_term - 2*(-1 + a^2)*(1 + sqrt_one_minus_a_sq) - 2*En^2*sqrt_one_minus_a_sq*(1 + sqrt_one_minus_a_sq)^2)))*omega - 
            2*(1 - a^2)*neg_horizon_term^2*sqrt_one_minus_a_sq*(m*sqrt_one_minus_a_sq*(a^3*(-1 + En^2) + 4*En*Lz*(1 + sqrt_one_minus_a_sq) - a*(-2 + Lz^2 - 2*sqrt_one_minus_a_sq + 4*En^2*(1 + sqrt_one_minus_a_sq))) - 2*(a^4*(1 + En^2) - 8*a*En*Lz*(1 + sqrt_one_minus_a_sq) + (4 + 16*En^2 + Lz^2)*(1 + sqrt_one_minus_a_sq) + a^3*En*Lz*(4 + sqrt_one_minus_a_sq) - a^2*(5 + 3*sqrt_one_minus_a_sq + En^2*(13 + 5*sqrt_one_minus_a_sq)))*omega))) + 
            (1/((1 - a^2)*neg_horizon_term^4*sqrt_one_minus_a_sq^3))*(8*I*(1 + sqrt_one_minus_a_sq)^2*(1 - a^2 + sqrt_one_minus_a_sq)*(-((-1 + a^2)*am_term*neg_horizon_term^4) + m*(1 + sqrt_one_minus_a_sq)*(-2*a*(-1 + a^2)*En*neg_horizon_term^3 - 4*(1 - a^2)*neg_horizon_term*sqrt_one_minus_a_sq*(a*En + Lz*sqrt_one_minus_a_sq)*(((-a)*En + Lz)^2*sqrt_one_minus_a_sq - 2*En*neg_horizon_term*(1 + sqrt_one_minus_a_sq) + sqrt_one_minus_a_sq*(1 + sqrt_one_minus_a_sq)^2) - 
            a*sqrt_one_minus_a_sq*(3*(1 - a^2)*(((-a)*En + Lz)^2*sqrt_one_minus_a_sq - 2*En*neg_horizon_term*(1 + sqrt_one_minus_a_sq) + sqrt_one_minus_a_sq*(1 + sqrt_one_minus_a_sq)^2)^2 + 2*neg_horizon_term^2*sqrt_one_minus_a_sq*(En*neg_horizon_term - 2*(-1 + a^2)*(1 + sqrt_one_minus_a_sq) - 2*En^2*sqrt_one_minus_a_sq*(1 + sqrt_one_minus_a_sq)^2))) + 
            2*sqrt_one_minus_a_sq*((-neg_horizon_term^3)*(4*a^3*Lz + 2*neg_horizon_term*sqrt_one_minus_a_sq + a*Lz*(-4 + sqrt_one_minus_a_sq + 3*sqrt_one_minus_a_sq^2) - 4*(-1 + a^2)*En*(-2 + 2*a^2 - 2*sqrt_one_minus_a_sq - sqrt_one_minus_a_sq^2 - sqrt_one_minus_a_sq^3)) + 2*(1 - a^2)*neg_horizon_term*(1 + sqrt_one_minus_a_sq)*(((-a)*En + Lz)^2*sqrt_one_minus_a_sq - 2*En*neg_horizon_term*(1 + sqrt_one_minus_a_sq) + sqrt_one_minus_a_sq*(1 + sqrt_one_minus_a_sq)^2)*
            (a*Lz + 2*neg_horizon_term + En*(2 - 2*a^2 + 3*sqrt_one_minus_a_sq + 2*sqrt_one_minus_a_sq^2 + sqrt_one_minus_a_sq^3)) + (1 + sqrt_one_minus_a_sq)^2*(3*(1 - a^2)*(((-a)*En + Lz)^2*sqrt_one_minus_a_sq - 2*En*neg_horizon_term*(1 + sqrt_one_minus_a_sq) + sqrt_one_minus_a_sq*(1 + sqrt_one_minus_a_sq)^2)^2 + 2*neg_horizon_term^2*sqrt_one_minus_a_sq*(En*neg_horizon_term - 2*(-1 + a^2)*(1 + sqrt_one_minus_a_sq) - 2*En^2*sqrt_one_minus_a_sq*(1 + sqrt_one_minus_a_sq)^2)))*omega - 
            2*(1 - a^2)*neg_horizon_term^2*sqrt_one_minus_a_sq*(m*sqrt_one_minus_a_sq*(a^3*(-1 + En^2) + 4*En*Lz*(1 + sqrt_one_minus_a_sq) - a*(-2 + Lz^2 - 2*sqrt_one_minus_a_sq + 4*En^2*(1 + sqrt_one_minus_a_sq))) - 2*(a^4*(1 + En^2) - 8*a*En*Lz*(1 + sqrt_one_minus_a_sq) + (4 + 16*En^2 + Lz^2)*(1 + sqrt_one_minus_a_sq) + a^3*En*Lz*(4 + sqrt_one_minus_a_sq) - a^2*(5 + 3*sqrt_one_minus_a_sq + En^2*(13 + 5*sqrt_one_minus_a_sq)))*omega))) + 
            (12*(1 + sqrt_one_minus_a_sq)^2*(-am_term^2 - 4*I*am_term*sqrt_one_minus_a_sq + (1 + sqrt_one_minus_a_sq)*(sqrt_one_minus_a_sq*lambda - 8*omega^2 - 8*sqrt_one_minus_a_sq*omega^2 + 4*a^2*omega*(I + omega) + 2*a*m*(-I + omega + sqrt_one_minus_a_sq*omega))))/sqrt_one_minus_a_sq^2 + 
            (2*(1 + sqrt_one_minus_a_sq)^2*(am_term^2 + 4*I*am_term*sqrt_one_minus_a_sq + 2*(1 + sqrt_one_minus_a_sq)*(-2*omega*(3*I + omega) - 2*sqrt_one_minus_a_sq*omega*(3*I + 10*omega) + 2*a^2*omega*(2*I + 6*I*sqrt_one_minus_a_sq + 2*omega + 7*sqrt_one_minus_a_sq*omega + sqrt_one_minus_a_sq^2*omega) + 2*sqrt_one_minus_a_sq^2*(lambda - 9*omega^2) + a*m*(I - omega + 4*sqrt_one_minus_a_sq^2*omega + 3*sqrt_one_minus_a_sq*(-I + omega)))))/((-1 + a^2)*sqrt_one_minus_a_sq) + 
            (36*I*(1 + sqrt_one_minus_a_sq)^2*((-am_term)*neg_horizon_term^2 + (1 + sqrt_one_minus_a_sq)*(long_term + a^3*(8*En*Lz*omega + sqrt_one_minus_a_sq*(m - En^2*m + 2*En*Lz*omega)))))/(neg_horizon_term^2*sqrt_one_minus_a_sq) - (80*I*(1 + sqrt_one_minus_a_sq)*(1 - a^2 + sqrt_one_minus_a_sq)*((-am_term)*neg_horizon_term^2 + (1 + sqrt_one_minus_a_sq)*(long_term + a^3*(8*En*Lz*omega + sqrt_one_minus_a_sq*(m - En^2*m + 2*En*Lz*omega)))))/(neg_horizon_term^2*sqrt_one_minus_a_sq^2) - 
            (12*I*(1 + sqrt_one_minus_a_sq)^2*(1 + 2*sqrt_one_minus_a_sq)*((-am_term)*neg_horizon_term^2 + (1 + sqrt_one_minus_a_sq)*(long_term + a^3*(8*En*Lz*omega + sqrt_one_minus_a_sq*(m - En^2*m + 2*En*Lz*omega)))))/(neg_horizon_term^2*sqrt_one_minus_a_sq^2) + (8*I*(1 + sqrt_one_minus_a_sq)*(-9*a^2 + 11*(1 + sqrt_one_minus_a_sq))*((-am_term)*neg_horizon_term^2 + (1 + sqrt_one_minus_a_sq)*(long_term + a^3*(8*En*Lz*omega + sqrt_one_minus_a_sq*(m - En^2*m + 2*En*Lz*omega)))))/(neg_horizon_term^2*sqrt_one_minus_a_sq^2) - 
            (1/sqrt_one_minus_a_sq^3)*(2*(1 + sqrt_one_minus_a_sq)^2*((sqrt_one_minus_a_sq^2*(2*I*am_term*sqrt_one_minus_a_sq + (a*m - 2*(1 + sqrt_one_minus_a_sq)*omega)^2))/(-1 + a^2) + (1/((1 - a^2)*neg_horizon_term^4))*((1 + sqrt_one_minus_a_sq)*(2*m^2*(1 + sqrt_one_minus_a_sq)*(4*a*(1 - a^2)*neg_horizon_term*sqrt_one_minus_a_sq*(a*En + Lz*sqrt_one_minus_a_sq)*(((-a)*En + Lz)^2*sqrt_one_minus_a_sq - 2*En*neg_horizon_term*(1 + sqrt_one_minus_a_sq) + sqrt_one_minus_a_sq*(1 + sqrt_one_minus_a_sq)^2) - 
            (1 - a^2)*neg_horizon_term^2*(a^2*En*neg_horizon_term - 2*sqrt_one_minus_a_sq*(a*En + Lz*sqrt_one_minus_a_sq)^2) + a^2*sqrt_one_minus_a_sq*(2*(1 - a^2)*(((-a)*En + Lz)^2*sqrt_one_minus_a_sq - 2*En*neg_horizon_term*(1 + sqrt_one_minus_a_sq) + sqrt_one_minus_a_sq*(1 + sqrt_one_minus_a_sq)^2)^2 + neg_horizon_term^2*sqrt_one_minus_a_sq*(En*neg_horizon_term - 2*(-1 + a^2)*(1 + sqrt_one_minus_a_sq) - 2*En^2*sqrt_one_minus_a_sq*(1 + sqrt_one_minus_a_sq)^2))) - 
            2*m*sqrt_one_minus_a_sq^2*(4*a^10*En*(-1 + En^2)*Lz + 64*En*Lz*(4 + Lz^2)*(1 + sqrt_one_minus_a_sq) - 64*a*(-4 + 4*En^4 - Lz^2 + En^2*(4 + 6*Lz^2))*(1 + sqrt_one_minus_a_sq) - 8*a^9*(-6 - sqrt_one_minus_a_sq + En^2*(8 + sqrt_one_minus_a_sq)) + 16*a^2*En*Lz*(-68 - 60*sqrt_one_minus_a_sq + 36*En^2*(1 + sqrt_one_minus_a_sq) - Lz^2*(7 + 5*sqrt_one_minus_a_sq)) - 
            4*a^8*En*Lz*(-34 + Lz^2 - 7*sqrt_one_minus_a_sq + En^2*(34 + 9*sqrt_one_minus_a_sq)) - 4*a^4*En*Lz*(36*En^2*(7 + 5*sqrt_one_minus_a_sq) - Lz^2*(13 + 5*sqrt_one_minus_a_sq) - 8*(47 + 33*sqrt_one_minus_a_sq)) + 4*a^6*En*Lz*(-201 - 95*sqrt_one_minus_a_sq + Lz^2*(2 + 3*sqrt_one_minus_a_sq) + En^2*(149 + 77*sqrt_one_minus_a_sq)) + 
            4*a^7*(-88 - 38*sqrt_one_minus_a_sq + Lz^2*(1 + sqrt_one_minus_a_sq) + 2*En^2*(62 + 27*sqrt_one_minus_a_sq + 3*Lz^2*(5 + sqrt_one_minus_a_sq))) + 16*a^3*(4*En^4*(5 + 3*sqrt_one_minus_a_sq) - 8*(6 + 5*sqrt_one_minus_a_sq) - Lz^2*(7 + 5*sqrt_one_minus_a_sq) + 4*En^2*(3*Lz^2*(4 + 3*sqrt_one_minus_a_sq) + 2*(7 + 6*sqrt_one_minus_a_sq))) - 
            a^5*(Lz^4*(1 + sqrt_one_minus_a_sq) + 16*En^4*(5 + sqrt_one_minus_a_sq) - 4*Lz^2*(11 + 3*sqrt_one_minus_a_sq) - 48*(17 + 11*sqrt_one_minus_a_sq) + 16*En^2*(67 + 45*sqrt_one_minus_a_sq + 3*Lz^2*(11 + 6*sqrt_one_minus_a_sq))))*omega + 4*sqrt_one_minus_a_sq*(1 + sqrt_one_minus_a_sq)*(4*(1 - a^2)*neg_horizon_term*(1 + sqrt_one_minus_a_sq)*(((-a)*En + Lz)^2*sqrt_one_minus_a_sq - 2*En*neg_horizon_term*(1 + sqrt_one_minus_a_sq) + sqrt_one_minus_a_sq*(1 + sqrt_one_minus_a_sq)^2)*
            (a*Lz + 2*neg_horizon_term + 4*En*(1 + sqrt_one_minus_a_sq) - a^2*En*(4 + sqrt_one_minus_a_sq)) + 2*(1 + sqrt_one_minus_a_sq)^2*(2*(1 - a^2)*(((-a)*En + Lz)^2*sqrt_one_minus_a_sq - 2*En*neg_horizon_term*(1 + sqrt_one_minus_a_sq) + sqrt_one_minus_a_sq*(1 + sqrt_one_minus_a_sq)^2)^2 + neg_horizon_term^2*sqrt_one_minus_a_sq*(En*neg_horizon_term - 2*(-1 + a^2)*(1 + sqrt_one_minus_a_sq) - 2*En^2*sqrt_one_minus_a_sq*(1 + sqrt_one_minus_a_sq)^2)) - 
            neg_horizon_term^2*(2*neg_horizon_term^2*(1 - 2*sqrt_one_minus_a_sq)*sqrt_one_minus_a_sq + 8*neg_horizon_term*sqrt_one_minus_a_sq^2*((-a)*Lz - 4*En*(1 + sqrt_one_minus_a_sq) + a^2*En*(4 + sqrt_one_minus_a_sq)) + neg_horizon_term*sqrt_one_minus_a_sq*(-4*(-1 + a^2)*En*(1 + sqrt_one_minus_a_sq) + a*Lz*(1 + 3*sqrt_one_minus_a_sq)) - (1 - a^2)*(a*Lz + En*(2 - 2*a^2 + 3*sqrt_one_minus_a_sq + 2*sqrt_one_minus_a_sq^2 + sqrt_one_minus_a_sq^3))^2))*omega^2 - 
            I*(m*(1 + sqrt_one_minus_a_sq)*(-4*(-1 + a^2)*Lz*neg_horizon_term^3*sqrt_one_minus_a_sq^2 + 4*(1 - a^2)*neg_horizon_term*sqrt_one_minus_a_sq^2*(a*En + Lz*sqrt_one_minus_a_sq)*(((-a)*En + Lz)^2*sqrt_one_minus_a_sq - 2*En*neg_horizon_term*(1 + sqrt_one_minus_a_sq) + sqrt_one_minus_a_sq*(1 + sqrt_one_minus_a_sq)^2) - 2*(1 - a^2)*neg_horizon_term*sqrt_one_minus_a_sq*(((-a)*En + Lz)^2*sqrt_one_minus_a_sq - 2*En*neg_horizon_term*(1 + sqrt_one_minus_a_sq) + sqrt_one_minus_a_sq*(1 + sqrt_one_minus_a_sq)^2)*
            (a*neg_horizon_term + 2*sqrt_one_minus_a_sq*(a*En + Lz*sqrt_one_minus_a_sq)) + neg_horizon_term^3*((-a)*(-1 + a^2)*neg_horizon_term + 2*a*En*sqrt_one_minus_a_sq^3 + 4*(-1 + a^2)*sqrt_one_minus_a_sq*(a*En + Lz*sqrt_one_minus_a_sq)) - a*sqrt_one_minus_a_sq^2*(3*(1 - a^2)*(((-a)*En + Lz)^2*sqrt_one_minus_a_sq - 2*En*neg_horizon_term*(1 + sqrt_one_minus_a_sq) + sqrt_one_minus_a_sq*(1 + sqrt_one_minus_a_sq)^2)^2 + 2*neg_horizon_term^2*sqrt_one_minus_a_sq*
            (En*neg_horizon_term - 2*(-1 + a^2)*(1 + sqrt_one_minus_a_sq) - 2*En^2*sqrt_one_minus_a_sq*(1 + sqrt_one_minus_a_sq)^2)) - (1 - a^2)*sqrt_one_minus_a_sq^2*(-2*a*neg_horizon_term^2*(6 + 5*a^2*(-1 + En^2) + Lz^2 + 6*sqrt_one_minus_a_sq - 12*En^2*(1 + sqrt_one_minus_a_sq)) - 6*a*(a^2*En^2*sqrt_one_minus_a_sq - 2*En*(neg_horizon_term + a*Lz*sqrt_one_minus_a_sq + neg_horizon_term*sqrt_one_minus_a_sq) + sqrt_one_minus_a_sq*(Lz^2 + (1 + sqrt_one_minus_a_sq)^2))^2 + 2*neg_horizon_term*(a*En + Lz*sqrt_one_minus_a_sq)*
            (4*En*neg_horizon_term*(1 + sqrt_one_minus_a_sq) - 2*sqrt_one_minus_a_sq*(((-a)*En + Lz)^2 + (1 + sqrt_one_minus_a_sq)^2)))) - 2*sqrt_one_minus_a_sq*(-2*(1 - a^2)*neg_horizon_term*(1 + sqrt_one_minus_a_sq)*(((-a)*En + Lz)^2*sqrt_one_minus_a_sq - 2*En*neg_horizon_term*(1 + sqrt_one_minus_a_sq) + sqrt_one_minus_a_sq*(1 + sqrt_one_minus_a_sq)^2)*(2*neg_horizon_term*sqrt_one_minus_a_sq - a*Lz*(1 + sqrt_one_minus_a_sq) + En*(4 - 2*a^2 + 4*sqrt_one_minus_a_sq)) - 
            sqrt_one_minus_a_sq*(1 + sqrt_one_minus_a_sq)^2*(3*(1 - a^2)*(((-a)*En + Lz)^2*sqrt_one_minus_a_sq - 2*En*neg_horizon_term*(1 + sqrt_one_minus_a_sq) + sqrt_one_minus_a_sq*(1 + sqrt_one_minus_a_sq)^2)^2 + 2*neg_horizon_term^2*sqrt_one_minus_a_sq*(En*neg_horizon_term - 2*(-1 + a^2)*(1 + sqrt_one_minus_a_sq) - 2*En^2*sqrt_one_minus_a_sq*(1 + sqrt_one_minus_a_sq)^2)) + neg_horizon_term^3*sqrt_one_minus_a_sq*(4*a^3*Lz + 2*neg_horizon_term*sqrt_one_minus_a_sq + a*Lz*(-5 + sqrt_one_minus_a_sq + 4*sqrt_one_minus_a_sq^2) - 2*En*
            (-2 + 2*a^4 + 6*sqrt_one_minus_a_sq + 8*sqrt_one_minus_a_sq^2 - a^2*(1 + 5*sqrt_one_minus_a_sq))) - 4*(-1 + a^2)*neg_horizon_term^2*sqrt_one_minus_a_sq*(a^4*(1 + En^2) - 4*a*En*Lz*(1 + sqrt_one_minus_a_sq) + (4 + 8*En^2 + Lz^2)*(1 + sqrt_one_minus_a_sq) + a^3*En*Lz*(4 + sqrt_one_minus_a_sq) - a^2*(5 + Lz^2 + 3*sqrt_one_minus_a_sq + En^2*(9 + 5*sqrt_one_minus_a_sq))) - 
            (1 - a^2)*sqrt_one_minus_a_sq*(1 + sqrt_one_minus_a_sq)*(-6*(1 + sqrt_one_minus_a_sq)*(a^2*En^2*sqrt_one_minus_a_sq - 2*En*(neg_horizon_term + a*Lz*sqrt_one_minus_a_sq + neg_horizon_term*sqrt_one_minus_a_sq) + sqrt_one_minus_a_sq*(Lz^2 + (1 + sqrt_one_minus_a_sq)^2))^2 + neg_horizon_term*(-2*neg_horizon_term*(1 + sqrt_one_minus_a_sq)*(6 + 5*a^2*(-1 + En^2) + Lz^2 + 6*sqrt_one_minus_a_sq - 12*En^2*(1 + sqrt_one_minus_a_sq)) + (a*Lz + 2*neg_horizon_term + En*(2 - 2*a^2 + 3*sqrt_one_minus_a_sq + 2*sqrt_one_minus_a_sq^2 + sqrt_one_minus_a_sq^3))*
            (4*En*neg_horizon_term*(1 + sqrt_one_minus_a_sq) - 2*sqrt_one_minus_a_sq*(((-a)*En + Lz)^2 + (1 + sqrt_one_minus_a_sq)^2)))))*omega))) + (1/(neg_horizon_term^2*(a^2*Lz^2 - 4*a*En*Lz*(1 + sqrt_one_minus_a_sq) + En^2*(8 - 4*a^2 + 8*sqrt_one_minus_a_sq))))*
            (4*sqrt_one_minus_a_sq*(1 + sqrt_one_minus_a_sq)*(a*m^2*sqrt_one_minus_a_sq*(a^2*Lz^2 - 4*a*En*Lz*(1 + sqrt_one_minus_a_sq) + En^2*(8 - 4*a^2 + 8*sqrt_one_minus_a_sq))*(a^3*(-1 + En^2) + 4*En*Lz*(1 + sqrt_one_minus_a_sq) - a*(-2 + Lz^2 - 2*sqrt_one_minus_a_sq + 4*En^2*(1 + sqrt_one_minus_a_sq))) - 
            2*m*(2*a^5 + a^3*(-10 - 4*En^2 + Lz^2 - 6*sqrt_one_minus_a_sq) + 8*a*(1 + En^2)*(1 + sqrt_one_minus_a_sq) + 8*En*Lz*(1 + sqrt_one_minus_a_sq) + a^4*En*Lz*(4 + sqrt_one_minus_a_sq) - 4*a^2*En*Lz*(4 + 3*sqrt_one_minus_a_sq))*(a^2*Lz^2 - 4*a*En*Lz*(1 + sqrt_one_minus_a_sq) + En^2*(8 - 4*a^2 + 8*sqrt_one_minus_a_sq))*omega + 
            4*(1 + sqrt_one_minus_a_sq)*(a^2*Lz^2 - 4*a*En*Lz*(1 + sqrt_one_minus_a_sq) + En^2*(8 - 4*a^2 + 8*sqrt_one_minus_a_sq))*(a^4*(1 + En^2) - 8*a*En*Lz*(1 + sqrt_one_minus_a_sq) + (4 + 16*En^2 + Lz^2)*(1 + sqrt_one_minus_a_sq) + a^3*En*Lz*(4 + sqrt_one_minus_a_sq) - a^2*(5 + 3*sqrt_one_minus_a_sq + En^2*(13 + 5*sqrt_one_minus_a_sq)))*omega^2 + 
            I*(-512*En^4*(1 + sqrt_one_minus_a_sq)*omega + 128*a*En^3*(1 + sqrt_one_minus_a_sq)*(En*m + 4*Lz*omega) + a^5*(16*En^4*m - 24*En^2*Lz^2*m + Lz^4*m + 64*En^3*Lz*omega - 16*En*Lz^3*omega) - 16*a^3*En*(-3*En*Lz^2*m*(1 + sqrt_one_minus_a_sq) + 4*En^3*m*(2 + sqrt_one_minus_a_sq) - 2*Lz^3*(1 + sqrt_one_minus_a_sq)*omega + 16*En^2*Lz*(2 + sqrt_one_minus_a_sq)*omega) - 
            2*a^4*(4*En*Lz^3*m*(1 + sqrt_one_minus_a_sq) - 16*En^3*Lz*m*(3 + sqrt_one_minus_a_sq) + Lz^4*(1 + sqrt_one_minus_a_sq)*omega - 24*En^2*Lz^2*(3 + sqrt_one_minus_a_sq)*omega + 16*En^4*(5 + sqrt_one_minus_a_sq)*omega) + 64*a^2*En^2*(-2*En*Lz*m*(1 + sqrt_one_minus_a_sq) - 3*Lz^2*(1 + sqrt_one_minus_a_sq)*omega + 2*En^2*(5 + 3*sqrt_one_minus_a_sq)*omega)))))) + 
            (1/sqrt_one_minus_a_sq)*(12*(1 + sqrt_one_minus_a_sq)^2*(2*I*am_term + (a*m - 2*(1 + sqrt_one_minus_a_sq)*omega)^2/sqrt_one_minus_a_sq + (1/neg_horizon_term^2)*((1 + sqrt_one_minus_a_sq)*(a*m^2*(a^3*(-1 + En^2) + 4*En*Lz*(1 + sqrt_one_minus_a_sq) - a*(-2 + Lz^2 - 2*sqrt_one_minus_a_sq + 4*En^2*(1 + sqrt_one_minus_a_sq))) - 
            (2*m*(2*a^5 + a^3*(-10 - 4*En^2 + Lz^2 - 6*sqrt_one_minus_a_sq) + 8*a*(1 + En^2)*(1 + sqrt_one_minus_a_sq) + 8*En*Lz*(1 + sqrt_one_minus_a_sq) + a^4*En*Lz*(4 + sqrt_one_minus_a_sq) - 4*a^2*En*Lz*(4 + 3*sqrt_one_minus_a_sq))*omega)/sqrt_one_minus_a_sq + 
            (4*(1 + sqrt_one_minus_a_sq)*(a^4*(1 + En^2) - 8*a*En*Lz*(1 + sqrt_one_minus_a_sq) + (4 + 16*En^2 + Lz^2)*(1 + sqrt_one_minus_a_sq) + a^3*En*Lz*(4 + sqrt_one_minus_a_sq) - a^2*(5 + 3*sqrt_one_minus_a_sq + En^2*(13 + 5*sqrt_one_minus_a_sq)))*omega^2)/sqrt_one_minus_a_sq + (1/(sqrt_one_minus_a_sq*(a^2*Lz^2 - 4*a*En*Lz*(1 + sqrt_one_minus_a_sq) + En^2*(8 - 4*a^2 + 8*sqrt_one_minus_a_sq))))*
            (I*(-512*En^4*(1 + sqrt_one_minus_a_sq)*omega + 128*a*En^3*(1 + sqrt_one_minus_a_sq)*(En*m + 4*Lz*omega) + a^5*(16*En^4*m - 24*En^2*Lz^2*m + Lz^4*m + 64*En^3*Lz*omega - 16*En*Lz^3*omega) - 16*a^3*En*(-3*En*Lz^2*m*(1 + sqrt_one_minus_a_sq) + 4*En^3*m*(2 + sqrt_one_minus_a_sq) - 2*Lz^3*(1 + sqrt_one_minus_a_sq)*omega + 16*En^2*Lz*(2 + sqrt_one_minus_a_sq)*omega) - 
            2*a^4*(4*En*Lz^3*m*(1 + sqrt_one_minus_a_sq) - 16*En^3*Lz*m*(3 + sqrt_one_minus_a_sq) + Lz^4*(1 + sqrt_one_minus_a_sq)*omega - 24*En^2*Lz^2*(3 + sqrt_one_minus_a_sq)*omega + 16*En^4*(5 + sqrt_one_minus_a_sq)*omega) + 64*a^2*En^2*(-2*En*Lz*m*(1 + sqrt_one_minus_a_sq) - 3*Lz^2*(1 + sqrt_one_minus_a_sq)*omega + 2*En^2*(5 + 3*sqrt_one_minus_a_sq)*omega))))))))
        end
        
        coefficient_for_scriptA1 = begin
            -((1/(2*(1 + sqrt_one_minus_a_sq)^6))*(30*I*am_term*(1 + sqrt_one_minus_a_sq)^2 + (20*I*(-1 + a^2)*am_term*(1 + sqrt_one_minus_a_sq)^2)/sqrt_one_minus_a_sq^2 + (5*am_term*(am_term + 4*I*sqrt_one_minus_a_sq)*(1 + sqrt_one_minus_a_sq)^2)/sqrt_one_minus_a_sq - 24*sqrt_one_minus_a_sq*(1 + sqrt_one_minus_a_sq)^2 + 24*(1 + sqrt_one_minus_a_sq)^2*(1 + 2*sqrt_one_minus_a_sq) - (6*I*am_term*(1 + sqrt_one_minus_a_sq)^2*(1 + 2*sqrt_one_minus_a_sq))/sqrt_one_minus_a_sq + (4*I*am_term*(1 + sqrt_one_minus_a_sq)^2*(1 + 5*sqrt_one_minus_a_sq))/sqrt_one_minus_a_sq + 
            (40*(-1 + a^2)*(a^2 - 2*(1 + sqrt_one_minus_a_sq)))/sqrt_one_minus_a_sq + (4*(1 + sqrt_one_minus_a_sq)*(-16*(1 + sqrt_one_minus_a_sq) + a^2*(16 + 13*sqrt_one_minus_a_sq)))/sqrt_one_minus_a_sq - (5*(1 + sqrt_one_minus_a_sq)^2*(2*I*am_term*sqrt_one_minus_a_sq + (a*m - 2*(1 + sqrt_one_minus_a_sq)*omega)^2))/sqrt_one_minus_a_sq + 
            (2*(1 + sqrt_one_minus_a_sq)^2*(-am_term^2 - 4*I*am_term*sqrt_one_minus_a_sq + (1 + sqrt_one_minus_a_sq)*(sqrt_one_minus_a_sq*lambda - 8*omega^2 - 8*sqrt_one_minus_a_sq*omega^2 + 4*a^2*omega*(I + omega) + 2*a*m*(-I + omega + sqrt_one_minus_a_sq*omega))))/sqrt_one_minus_a_sq + 
            (1/(neg_horizon_term^2*sqrt_one_minus_a_sq^2))*(4*I*(-1 + a^2)*(1 + sqrt_one_minus_a_sq)^2*((-am_term)*neg_horizon_term^2 + (1 + sqrt_one_minus_a_sq)*(2*a^4*(1 + En^2)*omega - 2*a^2*(5 + 3*sqrt_one_minus_a_sq + En^2*(13 + 5*sqrt_one_minus_a_sq))*omega + 2*(1 + sqrt_one_minus_a_sq)*(-2*En*Lz*m*sqrt_one_minus_a_sq + 16*En^2*omega + (4 + Lz^2)*omega) + 
            a*((-2 + 4*En^2 + Lz^2)*m*sqrt_one_minus_a_sq + 2*(-1 + 2*En^2)*m*sqrt_one_minus_a_sq^2 - 16*En*Lz*omega - 16*En*Lz*sqrt_one_minus_a_sq*omega) + a^3*(8*En*Lz*omega + sqrt_one_minus_a_sq*(m - En^2*m + 2*En*Lz*omega))))) - 
            6*I*(1 + sqrt_one_minus_a_sq)^2*(am_term - ((1 + sqrt_one_minus_a_sq)*(2*a^4*(1 + En^2)*omega - 2*a^2*(5 + 3*sqrt_one_minus_a_sq + En^2*(13 + 5*sqrt_one_minus_a_sq))*omega + 2*(1 + sqrt_one_minus_a_sq)*(-2*En*Lz*m*sqrt_one_minus_a_sq + 16*En^2*omega + (4 + Lz^2)*omega) + a*((-2 + 4*En^2 + Lz^2)*m*sqrt_one_minus_a_sq + 2*(-1 + 2*En^2)*m*sqrt_one_minus_a_sq^2 - 16*En*Lz*omega - 16*En*Lz*sqrt_one_minus_a_sq*omega) + 
            a^3*(8*En*Lz*omega + sqrt_one_minus_a_sq*(m - En^2*m + 2*En*Lz*omega))))/neg_horizon_term^2) + 2*(1 + sqrt_one_minus_a_sq)^2*(2*I*am_term + (a*m - 2*(1 + sqrt_one_minus_a_sq)*omega)^2/sqrt_one_minus_a_sq + (1/neg_horizon_term^2)*((1 + sqrt_one_minus_a_sq)*(a*m^2*(a^3*(-1 + En^2) + 4*En*Lz*(1 + sqrt_one_minus_a_sq) - a*(-2 + Lz^2 - 2*sqrt_one_minus_a_sq + 4*En^2*(1 + sqrt_one_minus_a_sq))) - 
            (2*m*(2*a^5 + a^3*(-10 - 4*En^2 + Lz^2 - 6*sqrt_one_minus_a_sq) + 8*a*(1 + En^2)*(1 + sqrt_one_minus_a_sq) + 8*En*Lz*(1 + sqrt_one_minus_a_sq) + a^4*En*Lz*(4 + sqrt_one_minus_a_sq) - 4*a^2*En*Lz*(4 + 3*sqrt_one_minus_a_sq))*omega)/sqrt_one_minus_a_sq + 
            (4*(1 + sqrt_one_minus_a_sq)*(a^4*(1 + En^2) - 8*a*En*Lz*(1 + sqrt_one_minus_a_sq) + (4 + 16*En^2 + Lz^2)*(1 + sqrt_one_minus_a_sq) + a^3*En*Lz*(4 + sqrt_one_minus_a_sq) - a^2*(5 + 3*sqrt_one_minus_a_sq + En^2*(13 + 5*sqrt_one_minus_a_sq)))*omega^2)/sqrt_one_minus_a_sq + (1/(sqrt_one_minus_a_sq*(a^2*Lz^2 - 4*a*En*Lz*(1 + sqrt_one_minus_a_sq) + En^2*(8 - 4*a^2 + 8*sqrt_one_minus_a_sq))))*
            (I*(-512*En^4*(1 + sqrt_one_minus_a_sq)*omega + 128*a*En^3*(1 + sqrt_one_minus_a_sq)*(En*m + 4*Lz*omega) + a^5*(16*En^4*m - 24*En^2*Lz^2*m + Lz^4*m + 64*En^3*Lz*omega - 16*En*Lz^3*omega) - 16*a^3*En*(-3*En*Lz^2*m*(1 + sqrt_one_minus_a_sq) + 4*En^3*m*(2 + sqrt_one_minus_a_sq) - 2*Lz^3*(1 + sqrt_one_minus_a_sq)*omega + 16*En^2*Lz*(2 + sqrt_one_minus_a_sq)*omega) - 
            2*a^4*(4*En*Lz^3*m*(1 + sqrt_one_minus_a_sq) - 16*En^3*Lz*m*(3 + sqrt_one_minus_a_sq) + Lz^4*(1 + sqrt_one_minus_a_sq)*omega - 24*En^2*Lz^2*(3 + sqrt_one_minus_a_sq)*omega + 16*En^4*(5 + sqrt_one_minus_a_sq)*omega) + 64*a^2*En^2*(-2*En*Lz*m*(1 + sqrt_one_minus_a_sq) - 3*Lz^2*(1 + sqrt_one_minus_a_sq)*omega + 2*En^2*(5 + 3*sqrt_one_minus_a_sq)*omega))))))))
        end

        return coefficient_for_scriptA0, coefficient_for_scriptA1
    else
        # Throw an error, this spin weight is not supported
        throw(DomainError(s, "Currently only spin weight s of +2 is supported"))
    end
end

function sourceterm_regularization_ansatz_coefficients(s::Int, l::Int, m::Int, a, omega, En, Lz, lambda)
    # Compute the necessary angular terms using SpinWeightedSpheroidalHarmonics.jl
    # There is a caching mechanism so only need to do spectral decomposition once
    swsh_piover2 = spin_weighted_spheroidal_harmonic(s, l, m, a*omega, pi/2, 0)
    psptheta_piover2 = spin_weighted_spheroidal_harmonic(s, l, m, a*omega, pi/2, 0; theta_derivative=1)
    p2sptheta2_piover2 = spin_weighted_spheroidal_harmonic(s, l, m, a*omega, pi/2, 0; theta_derivative=2)
    lambda = Teukolsky_lambda_const(a*omega, s, l, m)

    #=
    We have derived/shown the explicit expression for
    s = +2 ONLY
    =#
    if s == 2
        # s = +2
        # These alpha's are the series expansion coefficients of the RHS/source term
        alpha0 = _sourceterm_without_phasing_circularorbit_nearhorizon_seriesexpansion_zerothorder(s, l, m, a, omega, En, Lz; swsh_piover2=swsh_piover2, psptheta_piover2=psptheta_piover2, p2sptheta2_piover2=p2sptheta2_piover2)
        alpha1 = _sourceterm_without_phasing_circularorbit_nearhorizon_seriesexpansion_firstorder(s, l, m, a, omega, En, Lz; swsh_piover2=swsh_piover2, psptheta_piover2=psptheta_piover2, p2sptheta2_piover2=p2sptheta2_piover2)
    
        # These beta's are the series expansion coefficients of the LHS when the ansatz is substituted
        beta00, beta01 = _nearhorizon_ansatz_zerothorder(s, l, m, a, omega, En, Lz; swsh_piover2=swsh_piover2, psptheta_piover2=psptheta_piover2, p2sptheta2_piover2=p2sptheta2_piover2, lambda=lambda)
        beta10, beta11 = _nearhorizon_ansatz_firstorder(s, l, m, a, omega, En, Lz; swsh_piover2=swsh_piover2, psptheta_piover2=psptheta_piover2, p2sptheta2_piover2=p2sptheta2_piover2, lambda=lambda)

        # Solve the simple system of linear equations
        scriptA0 = (alpha1*beta01 - alpha0*beta11)/(beta01*beta10 - beta00*beta11)
        scriptA1 = (-alpha1*beta00 + alpha0*beta10)/(beta01*beta10 - beta00*beta11)

        return scriptA0, scriptA1
    else
        # Throw an error, this spin weight is not supported
        throw(DomainError(s, "Currently only spin weight s of +2 is supported"))
    end
end

function _lhs_without_phasing_ansatz(s::Int, l::Int, m::Int, a, omega, En, Lz, lambda, r, scriptA0, scriptA1)
    #=
    We have derived/shown the explicit expression for
    s = +2 ONLY
    =#
    if s == 2
        sqrt_term = sqrt(-((a^2 + (-2 + r)*r)*(((-a)*En + Lz)^2 + r^2)) + (a*Lz - En*(a^2 + r^2))^2)
        _Delta = Delta(a, r)
        _K = K(m, a, omega, r)
        denominator_term = (a*Lz - En*(a^2 + r^2))^2 - (((-a)*En + Lz)^2 + r^2)*_Delta
        VT_related_term = -_K^2 + 4*I*_K*(-1 + r) + _Delta*lambda - 8*I*r*_Delta*omega
        three_im_term = (m*(a*(a*Lz - 2*En*r) - Lz*_Delta))/(sqrt_term*_Delta) + (r^2*(En + (2*(a^2*En - a*Lz + En*r^2))/(r*_Delta))*omega)/sqrt_term
        coefficient_m_term = begin
            (-2*a*En - Lz*(-2 + 2*r))/(sqrt_term*_Delta) - ((-2 + 2*r)*(a*(a*Lz - 2*En*r) - Lz*_Delta))/(sqrt_term*_Delta^2) - 
            ((a*(a*Lz - 2*En*r) - Lz*_Delta)*(-((-2 + 2*r)*(((-a)*En + Lz)^2 + r^2)) - 4*En*r*(a*Lz - En*(a^2 + r^2)) - 2*r*_Delta))/(2*denominator_term^(3/2)*_Delta)
        end
        coefficient_omega_term = begin
            -((r^2*(-((2*(-2 + 2*r)*(a^2*En - a*Lz + En*r^2))/(r*_Delta^2)) + (4*En)/_Delta - (2*(a^2*En - a*Lz + En*r^2))/(r^2*_Delta)))/sqrt_term) - (2*r*(En + (2*(a^2*En - a*Lz + En*r^2))/(r*_Delta)))/sqrt_term + 
            (r^2*(En + (2*(a^2*En - a*Lz + En*r^2))/(r*_Delta))*(-((-2 + 2*r)*(((-a)*En + Lz)^2 + r^2)) - 4*En*r*(a*Lz - En*(a^2 + r^2)) - 2*r*_Delta))/(2*denominator_term^(3/2))
        end
        r_Delta_term = begin
            (m^2*(a*(a*Lz - 2*En*r) - Lz*_Delta)^2)/(denominator_term*_Delta^2) + (2*m*r^2*(En + (2*(a^2*En - a*Lz + En*r^2))/(r*_Delta))*(a*(a*Lz - 2*En*r) - Lz*_Delta)*omega)/(denominator_term*_Delta) + 
            (r^4*(En + (2*(a^2*En - a*Lz + En*r^2))/(r*_Delta))^2*omega^2)/denominator_term + I*(coefficient_m_term*m - coefficient_omega_term*omega)
        end

        coefficient_for_scriptA0 = begin
            -((-12*(2 - 2*r)*(-1 + r)*r^4 - 2*r^2*(a^4 + 4*a^2*(-2 + r)*r + 3*r^2*(8 - 12*r + 5*r^2)) + r^4*VT_related_term - 3*(2 - 2*r)*r^3*_Delta - 3*I*three_im_term*(2 - 2*r)*r^4*_Delta - 
            2*I*three_im_term*r^3*(a^2 + r*(-6 + 5*r))*_Delta + r^4*r_Delta_term*_Delta^2)/(r^5*_Delta^3))
        end

        coefficient_for_scriptA1 = begin
            -((-6*(2 - 2*r)*(-1 + r)*r^2*_Delta - 2*I*three_im_term*r*(3*a^4 + 2*a^2*r*(-7 + 4*r) + r^2*(16 - 18*r + 5*r^2))*_Delta - 2*(6*a^4 + a^2*r*(-30 + 17*r) + r^2*(40 - 48*r + 15*r^2))*_Delta + 
            r^2*VT_related_term*_Delta - 9*(2 - 2*r)*r*_Delta^2 - 3*I*three_im_term*(2 - 2*r)*r^2*_Delta^2 + r^2*r_Delta_term*_Delta^3)/(r^5*_Delta^3))
        end

        return coefficient_for_scriptA0*scriptA0 + coefficient_for_scriptA1*scriptA1
    else
        # Throw an error, this spin weight is not supported
        throw(DomainError(s, "Currently only spin weight s of +2 is supported"))
    end
end
