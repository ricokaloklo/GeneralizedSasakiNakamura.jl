using SpinWeightedSpheroidalHarmonics

function sourceterm_without_phasing_circularorbit_nearhorizon_seriesexpansion_zerothorder(s::Int, l::Int, m::Int, a, omega, En, Lz, r)
    # Compute the necessary angular terms using SpinWeightedSpheroidalHarmonics.jl
    # There is a caching mechanism so only need to do spectral decomposition once
    swsh_piover2 = spin_weighted_spheroidal_harmonic(s, l, m, a*omega, pi/2, 0)
    psptheta_piover2 = spin_weighted_spheroidal_harmonic(s, l, m, a*omega, pi/2, 0; theta_derivative=1)
    p2sptheta2_piover2 = spin_weighted_spheroidal_harmonic(s, l, m, a*omega, pi/2, 0; theta_derivative=2)

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

function sourceterm_without_phasing_circularorbit_nearhorizon_seriesexpansion_firstorder(s::Int, l::Int, m::Int, a, omega, En, Lz, r)
    # Compute the necessary angular terms using SpinWeightedSpheroidalHarmonics.jl
    # There is a caching mechanism so only need to do spectral decomposition once
    swsh_piover2 = spin_weighted_spheroidal_harmonic(s, l, m, a*omega, pi/2, 0)
    psptheta_piover2 = spin_weighted_spheroidal_harmonic(s, l, m, a*omega, pi/2, 0; theta_derivative=1)
    p2sptheta2_piover2 = spin_weighted_spheroidal_harmonic(s, l, m, a*omega, pi/2, 0; theta_derivative=2)

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