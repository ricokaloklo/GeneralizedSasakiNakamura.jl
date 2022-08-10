using SpinWeightedSpheroidalHarmonics

# Mathematica being Mathematica
const Pi = pi
const I = 1im

function sourceterm_without_phasing_circularorbit(s::Int, l::Int, m::Int, a, omega, En, Lz, r)
    # Compute the necessary angular terms using SpinWeightedSpheroidalHarmonics.jl
    # There is a caching mechanism so only need to do spectral decomposition once
    swsh_piover2 = spin_weighted_spheroidal_harmonic(s, l, m, a*omega, pi/2, 0)
    psptheta_piover2 = spin_weighted_spheroidal_harmonic(s, l, m, a*omega, pi/2, 0; theta_derivative=1)
    p2sptheta2_piover2 = spin_weighted_spheroidal_harmonic(s, l, m, a*omega, pi/2, 0; theta_derivative=2)

    #=
    We have derived/shown the explicit expression for
    s = \pm 2
    =#
    if s == 2
        # s = +2
        return begin
            -((2*((-a)*En + Lz)^2*r^2*(En + (2*(a^2*En - a*Lz + En*r^2))/(r*(a^2 + (-2 + r)*r)))*((2*((-4*a*Lz + 3*a^2*En*r + 8*En*r^3 + a^2*En*(2 + 3*r))/(a^2 + (-2 + r)*r) - ((-2 + 2*r)*(a^4*En - 4*a*Lz*r + 2*En*r^4 + a^2*En*r*(2 + 3*r)))/(a^2 + (-2 + r)*r)^2)^2)/
            ((-a^2)*En + (a^4*En - 4*a*Lz*r + 2*En*r^4 + a^2*En*r*(2 + 3*r))/(a^2 + (-2 + r)*r))^3 - ((6*a^2*En + 24*En*r^2)/(a^2 + (-2 + r)*r) - (2*(-2 + 2*r)*(-4*a*Lz + 3*a^2*En*r + 8*En*r^3 + a^2*En*(2 + 3*r)))/(a^2 + (-2 + r)*r)^2 + 
            (2*(-2 + 2*r)^2*(a^4*En - 4*a*Lz*r + 2*En*r^4 + a^2*En*r*(2 + 3*r)))/(a^2 + (-2 + r)*r)^3 - (2*(a^4*En - 4*a*Lz*r + 2*En*r^4 + a^2*En*r*(2 + 3*r)))/(a^2 + (-2 + r)*r)^2)/((-a^2)*En + (a^4*En - 4*a*Lz*r + 2*En*r^4 + a^2*En*r*(2 + 3*r))/(a^2 + (-2 + r)*r))^2)*swsh_piover2)/
            (Pi*sqrt(-((a^2 + (-2 + r)*r)*(((-a)*En + Lz)^2 + r^2)) + (a*Lz - En*(a^2 + r^2))^2))) - (4*((-a)*En + Lz)^2*((-4*a*Lz + 3*a^2*En*r + 8*En*r^3 + a^2*En*(2 + 3*r))/(a^2 + (-2 + r)*r) - ((-2 + 2*r)*(a^4*En - 4*a*Lz*r + 2*En*r^4 + a^2*En*r*(2 + 3*r)))/(a^2 + (-2 + r)*r)^2)*
            ((r^2*(En + (2*(a^2*En - a*Lz + En*r^2))/(r*(a^2 + (-2 + r)*r)))*(-2*r*(a^2 + (-2 + r)*r) - (-2 + 2*r)*(((-a)*En + Lz)^2 + r^2) - 4*En*r*(a*Lz - En*(a^2 + r^2))))/(2*(-((a^2 + (-2 + r)*r)*(((-a)*En + Lz)^2 + r^2)) + (a*Lz - En*(a^2 + r^2))^2)^(3/2)) + 
            (I*m*r^2*(a*(a*Lz - 2*En*r) - Lz*(a^2 + (-2 + r)*r))*(En + (2*(a^2*En - a*Lz + En*r^2))/(r*(a^2 + (-2 + r)*r))))/((a^2 + (-2 + r)*r)*(-((a^2 + (-2 + r)*r)*(((-a)*En + Lz)^2 + r^2)) + (a*Lz - En*(a^2 + r^2))^2)) - 
            (r^2*((4*En)/(a^2 + (-2 + r)*r) - (2*(-2 + 2*r)*(a^2*En - a*Lz + En*r^2))/(r*(a^2 + (-2 + r)*r)^2) - (2*(a^2*En - a*Lz + En*r^2))/(r^2*(a^2 + (-2 + r)*r))))/sqrt(-((a^2 + (-2 + r)*r)*(((-a)*En + Lz)^2 + r^2)) + (a*Lz - En*(a^2 + r^2))^2) - 
            (2*r*(En + (2*(a^2*En - a*Lz + En*r^2))/(r*(a^2 + (-2 + r)*r))))/sqrt(-((a^2 + (-2 + r)*r)*(((-a)*En + Lz)^2 + r^2)) + (a*Lz - En*(a^2 + r^2))^2) + (I*r^4*(En + (2*(a^2*En - a*Lz + En*r^2))/(r*(a^2 + (-2 + r)*r)))^2*omega)/(-((a^2 + (-2 + r)*r)*(((-a)*En + Lz)^2 + r^2)) + (a*Lz - En*(a^2 + r^2))^2))*
            swsh_piover2)/(Pi*((-a^2)*En + (a^4*En - 4*a*Lz*r + 2*En*r^4 + a^2*En*r*(2 + 3*r))/(a^2 + (-2 + r)*r))^2) + 
            (2*((-a)*En + Lz)^2*((m^2*r^2*(a*(a*Lz - 2*En*r) - Lz*(a^2 + (-2 + r)*r))^2*(En + (2*(a^2*En - a*Lz + En*r^2))/(r*(a^2 + (-2 + r)*r))))/((a^2 + (-2 + r)*r)^2*(-((a^2 + (-2 + r)*r)*(((-a)*En + Lz)^2 + r^2)) + (a*Lz - En*(a^2 + r^2))^2)^(3/2)) + 
            ((-2*r*(a^2 + (-2 + r)*r) - (-2 + 2*r)*(((-a)*En + Lz)^2 + r^2) - 4*En*r*(a*Lz - En*(a^2 + r^2)))*(r^2*((4*En)/(a^2 + (-2 + r)*r) - (2*(-2 + 2*r)*(a^2*En - a*Lz + En*r^2))/(r*(a^2 + (-2 + r)*r)^2) - (2*(a^2*En - a*Lz + En*r^2))/(r^2*(a^2 + (-2 + r)*r))) + 
            2*r*(En + (2*(a^2*En - a*Lz + En*r^2))/(r*(a^2 + (-2 + r)*r)))))/(-((a^2 + (-2 + r)*r)*(((-a)*En + Lz)^2 + r^2)) + (a*Lz - En*(a^2 + r^2))^2)^(3/2) + r^2*(En + (2*(a^2*En - a*Lz + En*r^2))/(r*(a^2 + (-2 + r)*r)))*
            (-((3*(-2*r*(a^2 + (-2 + r)*r) - (-2 + 2*r)*(((-a)*En + Lz)^2 + r^2) - 4*En*r*(a*Lz - En*(a^2 + r^2)))^2)/(4*(-((a^2 + (-2 + r)*r)*(((-a)*En + Lz)^2 + r^2)) + (a*Lz - En*(a^2 + r^2))^2)^(5/2))) + 
            (8*En^2*r^2 - 4*r*(-2 + 2*r) - 2*(a^2 + (-2 + r)*r) - 2*(((-a)*En + Lz)^2 + r^2) - 4*En*(a*Lz - En*(a^2 + r^2)))/(2*(-((a^2 + (-2 + r)*r)*(((-a)*En + Lz)^2 + r^2)) + (a*Lz - En*(a^2 + r^2))^2)^(3/2))) + 
            (I*m*r^2*(En + (2*(a^2*En - a*Lz + En*r^2))/(r*(a^2 + (-2 + r)*r)))*(-(((a*(a*Lz - 2*En*r) - Lz*(a^2 + (-2 + r)*r))*(-2*r*(a^2 + (-2 + r)*r) - (-2 + 2*r)*(((-a)*En + Lz)^2 + r^2) - 4*En*r*(a*Lz - En*(a^2 + r^2))))/
            (2*(a^2 + (-2 + r)*r)*(-((a^2 + (-2 + r)*r)*(((-a)*En + Lz)^2 + r^2)) + (a*Lz - En*(a^2 + r^2))^2)^(3/2))) + (-2*a*En - Lz*(-2 + 2*r))/((a^2 + (-2 + r)*r)*sqrt(-((a^2 + (-2 + r)*r)*(((-a)*En + Lz)^2 + r^2)) + (a*Lz - En*(a^2 + r^2))^2)) - 
            ((-2 + 2*r)*(a*(a*Lz - 2*En*r) - Lz*(a^2 + (-2 + r)*r)))/((a^2 + (-2 + r)*r)^2*sqrt(-((a^2 + (-2 + r)*r)*(((-a)*En + Lz)^2 + r^2)) + (a*Lz - En*(a^2 + r^2))^2))))/sqrt(-((a^2 + (-2 + r)*r)*(((-a)*En + Lz)^2 + r^2)) + (a*Lz - En*(a^2 + r^2))^2) - 
            (2*I*m*(a*(a*Lz - 2*En*r) - Lz*(a^2 + (-2 + r)*r))*((r^2*(En + (2*(a^2*En - a*Lz + En*r^2))/(r*(a^2 + (-2 + r)*r)))*(-2*r*(a^2 + (-2 + r)*r) - (-2 + 2*r)*(((-a)*En + Lz)^2 + r^2) - 4*En*r*(a*Lz - En*(a^2 + r^2))))/(2*(-((a^2 + (-2 + r)*r)*(((-a)*En + Lz)^2 + r^2)) + (a*Lz - En*(a^2 + r^2))^2)^(3/2)) - 
            (r^2*((4*En)/(a^2 + (-2 + r)*r) - (2*(-2 + 2*r)*(a^2*En - a*Lz + En*r^2))/(r*(a^2 + (-2 + r)*r)^2) - (2*(a^2*En - a*Lz + En*r^2))/(r^2*(a^2 + (-2 + r)*r))))/sqrt(-((a^2 + (-2 + r)*r)*(((-a)*En + Lz)^2 + r^2)) + (a*Lz - En*(a^2 + r^2))^2) - 
            (2*r*(En + (2*(a^2*En - a*Lz + En*r^2))/(r*(a^2 + (-2 + r)*r))))/sqrt(-((a^2 + (-2 + r)*r)*(((-a)*En + Lz)^2 + r^2)) + (a*Lz - En*(a^2 + r^2))^2)))/((a^2 + (-2 + r)*r)*sqrt(-((a^2 + (-2 + r)*r)*(((-a)*En + Lz)^2 + r^2)) + (a*Lz - En*(a^2 + r^2))^2)) - 
            (4*r*((4*En)/(a^2 + (-2 + r)*r) - (2*(-2 + 2*r)*(a^2*En - a*Lz + En*r^2))/(r*(a^2 + (-2 + r)*r)^2) - (2*(a^2*En - a*Lz + En*r^2))/(r^2*(a^2 + (-2 + r)*r))) + 2*(En + (2*(a^2*En - a*Lz + En*r^2))/(r*(a^2 + (-2 + r)*r))) + 
            r^2*((4*r*(-((3*En)/r^2) + (3*(a^2*En - a*Lz + En*r^2))/r^4))/(a^2 + (-2 + r)*r) + 2*(-((4*r*(-2 + 2*r))/(a^2 + (-2 + r)*r)^2) + 4/(a^2 + (-2 + r)*r))*(En/r - (a^2*En - a*Lz + En*r^2)/r^3) + 
            ((a^2*En - a*Lz + En*r^2)*(-((8*(-2 + 2*r))/(a^2 + (-2 + r)*r)^2) + 4*r*((2*(-2 + 2*r)^2)/(a^2 + (-2 + r)*r)^3 - 2/(a^2 + (-2 + r)*r)^2)))/(2*r^2)))/sqrt(-((a^2 + (-2 + r)*r)*(((-a)*En + Lz)^2 + r^2)) + (a*Lz - En*(a^2 + r^2))^2) + 
            (2*m*r^4*(a*(a*Lz - 2*En*r) - Lz*(a^2 + (-2 + r)*r))*(En + (2*(a^2*En - a*Lz + En*r^2))/(r*(a^2 + (-2 + r)*r)))^2*omega)/((a^2 + (-2 + r)*r)*(-((a^2 + (-2 + r)*r)*(((-a)*En + Lz)^2 + r^2)) + (a*Lz - En*(a^2 + r^2))^2)^(3/2)) - 
            (3*I*r^2*(En + (2*(a^2*En - a*Lz + En*r^2))/(r*(a^2 + (-2 + r)*r)))*((r^2*(En + (2*(a^2*En - a*Lz + En*r^2))/(r*(a^2 + (-2 + r)*r)))*(-2*r*(a^2 + (-2 + r)*r) - (-2 + 2*r)*(((-a)*En + Lz)^2 + r^2) - 4*En*r*(a*Lz - En*(a^2 + r^2))))/
            (2*(-((a^2 + (-2 + r)*r)*(((-a)*En + Lz)^2 + r^2)) + (a*Lz - En*(a^2 + r^2))^2)^(3/2)) - (r^2*((4*En)/(a^2 + (-2 + r)*r) - (2*(-2 + 2*r)*(a^2*En - a*Lz + En*r^2))/(r*(a^2 + (-2 + r)*r)^2) - (2*(a^2*En - a*Lz + En*r^2))/(r^2*(a^2 + (-2 + r)*r))))/
            sqrt(-((a^2 + (-2 + r)*r)*(((-a)*En + Lz)^2 + r^2)) + (a*Lz - En*(a^2 + r^2))^2) - (2*r*(En + (2*(a^2*En - a*Lz + En*r^2))/(r*(a^2 + (-2 + r)*r))))/sqrt(-((a^2 + (-2 + r)*r)*(((-a)*En + Lz)^2 + r^2)) + (a*Lz - En*(a^2 + r^2))^2))*omega)/
            sqrt(-((a^2 + (-2 + r)*r)*(((-a)*En + Lz)^2 + r^2)) + (a*Lz - En*(a^2 + r^2))^2) + (r^6*(En + (2*(a^2*En - a*Lz + En*r^2))/(r*(a^2 + (-2 + r)*r)))^3*omega^2)/(-((a^2 + (-2 + r)*r)*(((-a)*En + Lz)^2 + r^2)) + (a*Lz - En*(a^2 + r^2))^2)^(3/2))*swsh_piover2)/
            (Pi*((-a^2)*En + (a^4*En - 4*a*Lz*r + 2*En*r^4 + a^2*En*r*(2 + 3*r))/(a^2 + (-2 + r)*r))) - 
            (r^2*(En + (2*(a^2*En - a*Lz + En*r^2))/(r*(a^2 + (-2 + r)*r)))*(-((4*((-a)*En + Lz)^2*(1/r^2 + (2*I*r*omega)/(a^2 - 2*r + r^2) - (I*(-2 + 2*r)*((-a)*m + (a^2 + r^2)*omega))/(a^2 - 2*r + r^2)^2)*swsh_piover2)/(Pi*((-a^2)*En + (a^4*En - 4*a*Lz*r + 2*En*r^4 + a^2*En*r*(2 + 3*r))/(a^2 + (-2 + r)*r)))) + 
            (4*((-a)*En + Lz)^2*((-4*a*Lz + 3*a^2*En*r + 8*En*r^3 + a^2*En*(2 + 3*r))/(a^2 + (-2 + r)*r) - ((-2 + 2*r)*(a^4*En - 4*a*Lz*r + 2*En*r^4 + a^2*En*r*(2 + 3*r)))/(a^2 + (-2 + r)*r)^2)*(-(1/r) + (I*((-a)*m + (a^2 + r^2)*omega))/(a^2 - 2*r + r^2))*swsh_piover2)/
            (Pi*((-a^2)*En + (a^4*En - 4*a*Lz*r + 2*En*r^4 + a^2*En*r*(2 + 3*r))/(a^2 + (-2 + r)*r))^2) + (2*I*(-2*a*En + 2*Lz)*(2*En*r + (-2*r*(a^2 + (-2 + r)*r) - (-2 + 2*r)*(((-a)*En + Lz)^2 + r^2) - 4*En*r*(a*Lz - En*(a^2 + r^2)))/
            (2*sqrt(-((a^2 + (-2 + r)*r)*(((-a)*En + Lz)^2 + r^2)) + (a*Lz - En*(a^2 + r^2))^2)))*(psptheta_piover2 + (m - a*omega)*swsh_piover2))/(Pi*(a^4*En - 4*a*Lz*r + 2*En*r^4 + a^2*En*r*(2 + 3*r) - a^2*En*(a^2 + (-2 + r)*r))) - 
            (2*I*(-2*a*En + 2*Lz)*(-4*a*Lz + 3*a^2*En*r + 8*En*r^3 - a^2*En*(-2 + 2*r) + a^2*En*(2 + 3*r))*(a^2*En - a*Lz + En*r^2 + sqrt(-((a^2 + (-2 + r)*r)*(((-a)*En + Lz)^2 + r^2)) + (a*Lz - En*(a^2 + r^2))^2))*(psptheta_piover2 + (m - a*omega)*swsh_piover2))/
            (Pi*(a^4*En - 4*a*Lz*r + 2*En*r^4 + a^2*En*r*(2 + 3*r) - a^2*En*(a^2 + (-2 + r)*r))^2)))/sqrt(-((a^2 + (-2 + r)*r)*(((-a)*En + Lz)^2 + r^2)) + (a*Lz - En*(a^2 + r^2))^2) + 
            ((r^2*(En + (2*(a^2*En - a*Lz + En*r^2))/(r*(a^2 + (-2 + r)*r)))*(-2*r*(a^2 + (-2 + r)*r) - (-2 + 2*r)*(((-a)*En + Lz)^2 + r^2) - 4*En*r*(a*Lz - En*(a^2 + r^2))))/(2*(-((a^2 + (-2 + r)*r)*(((-a)*En + Lz)^2 + r^2)) + (a*Lz - En*(a^2 + r^2))^2)^(3/2)) + 
            (I*m*r^2*(a*(a*Lz - 2*En*r) - Lz*(a^2 + (-2 + r)*r))*(En + (2*(a^2*En - a*Lz + En*r^2))/(r*(a^2 + (-2 + r)*r))))/((a^2 + (-2 + r)*r)*(-((a^2 + (-2 + r)*r)*(((-a)*En + Lz)^2 + r^2)) + (a*Lz - En*(a^2 + r^2))^2)) - 
            (r^2*((4*En)/(a^2 + (-2 + r)*r) - (2*(-2 + 2*r)*(a^2*En - a*Lz + En*r^2))/(r*(a^2 + (-2 + r)*r)^2) - (2*(a^2*En - a*Lz + En*r^2))/(r^2*(a^2 + (-2 + r)*r))))/sqrt(-((a^2 + (-2 + r)*r)*(((-a)*En + Lz)^2 + r^2)) + (a*Lz - En*(a^2 + r^2))^2) - 
            (2*r*(En + (2*(a^2*En - a*Lz + En*r^2))/(r*(a^2 + (-2 + r)*r))))/sqrt(-((a^2 + (-2 + r)*r)*(((-a)*En + Lz)^2 + r^2)) + (a*Lz - En*(a^2 + r^2))^2) + (I*r^4*(En + (2*(a^2*En - a*Lz + En*r^2))/(r*(a^2 + (-2 + r)*r)))^2*omega)/(-((a^2 + (-2 + r)*r)*(((-a)*En + Lz)^2 + r^2)) + (a*Lz - En*(a^2 + r^2))^2))*
            (-((4*((-a)*En + Lz)^2*(-(1/r) + (I*((-a)*m + (a^2 + r^2)*omega))/(a^2 - 2*r + r^2))*swsh_piover2)/(Pi*((-a^2)*En + (a^4*En - 4*a*Lz*r + 2*En*r^4 + a^2*En*r*(2 + 3*r))/(a^2 + (-2 + r)*r)))) + 
            (2*I*(-2*a*En + 2*Lz)*(a^2*En - a*Lz + En*r^2 + sqrt(-((a^2 + (-2 + r)*r)*(((-a)*En + Lz)^2 + r^2)) + (a*Lz - En*(a^2 + r^2))^2))*(psptheta_piover2 + (m - a*omega)*swsh_piover2))/(Pi*(a^4*En - 4*a*Lz*r + 2*En*r^4 + a^2*En*r*(2 + 3*r) - a^2*En*(a^2 + (-2 + r)*r)))) - 
            (r^2*(En + (2*(a^2*En - a*Lz + En*r^2))/(r*(a^2 + (-2 + r)*r)))*((2*((-a)*En + Lz)^2*(-((2*I*((-a)*m + (a^2 + r^2)*omega))/(r*(a^2 - 2*r + r^2))) - ((-a)*m + (a^2 + r^2)*omega)^2/(a^2 - 2*r + r^2)^2 + I*((2*r*omega)/(a^2 - 2*r + r^2) - ((-2 + 2*r)*((-a)*m + (a^2 + r^2)*omega))/(a^2 - 2*r + r^2)^2))*swsh_piover2)/
            (Pi*((-a^2)*En + (a^4*En - 4*a*Lz*r + 2*En*r^4 + a^2*En*r*(2 + 3*r))/(a^2 + (-2 + r)*r))) + (2*I*(-2*a*En + 2*Lz)*(a^2*En - a*Lz + En*r^2 + sqrt(-((a^2 + (-2 + r)*r)*(((-a)*En + Lz)^2 + r^2)) + (a*Lz - En*(a^2 + r^2))^2))*(2/r - (I*((-a)*m + (a^2 + r^2)*omega))/(a^2 - 2*r + r^2))*
            (psptheta_piover2 + (m - a*omega)*swsh_piover2))/(Pi*(a^4*En - 4*a*Lz*r + 2*En*r^4 + a^2*En*r*(2 + 3*r) - a^2*En*(a^2 + (-2 + r)*r))) + (2*(a^2*En - a*Lz + En*r^2 + sqrt(-((a^2 + (-2 + r)*r)*(((-a)*En + Lz)^2 + r^2)) + (a*Lz - En*(a^2 + r^2))^2))^2*
            ((-psptheta_piover2)*r*(m - a*omega) - r*(m - a*omega)^2*swsh_piover2 + (2*r - I*a*(m - a*omega))*swsh_piover2 + r^4*((6*I*a*psptheta_piover2)/r^4 - p2sptheta2_piover2/r^3 + (12*a^2*swsh_piover2)/r^5) + 4*I*a*r^3*(-(psptheta_piover2/r^3) + (3*I*a*swsh_piover2)/r^4) + r^4*(m - a*omega)*(-(psptheta_piover2/r^3) + (3*I*a*swsh_piover2)/r^4)))/
            (Pi*r*(a^2 + (-2 + r)*r)*(a^4*En - 4*a*Lz*r + 2*En*r^4 + a^2*En*r*(2 + 3*r) - a^2*En*(a^2 + (-2 + r)*r)))))/sqrt(-((a^2 + (-2 + r)*r)*(((-a)*En + Lz)^2 + r^2)) + (a*Lz - En*(a^2 + r^2))^2)
        end
    elseif s == -2
        # s = -2
        return nothing
    else
        # Throw an error, this spin weight is not supported
        throw(DomainError(s, "Currently only spin weight s of +/-2 are supported"))
    end
end