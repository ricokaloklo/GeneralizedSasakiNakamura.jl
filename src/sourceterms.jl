using SpinWeightedSpheroidalHarmonics

include("kerr.jl")

# Mathematica being Mathematica
const Pi = pi
const I = 1im

function sourceterm_without_phasing_circularorbit(s::Int, l::Int, m::Int, a, omega, En, Lz, r; swsh_piover2=nothing, psptheta_piover2=nothing, p2sptheta2_piover2=nothing)
    if isnothing(swsh_piover2) || isnothing(psptheta_piover2) || isnothing(p2sptheta2_piover2)
        swsh_piover2 = spin_weighted_spheroidal_harmonic(s, l, m, a*omega, pi/2, 0)
        psptheta_piover2 = spin_weighted_spheroidal_harmonic(s, l, m, a*omega, pi/2, 0; theta_derivative=1)
        p2sptheta2_piover2 = spin_weighted_spheroidal_harmonic(s, l, m, a*omega, pi/2, 0; theta_derivative=2)
    end

    #=
    We have derived/shown the explicit expression for
    s = \pm 2
    =#
    if s == 2
        # s = +2
        _Delta = Delta(a, r)
        sqrtterm = sqrt((a*Lz - En*(a^2 + r^2))^2 - (((-a)*En + Lz)^2 + r^2)*_Delta)
        a4Enterm = a^4*En - 4*a*Lz*r + 2*En*r^4 + a^2*En*r*(2 + 3*r)
        a2Enterm = a^2*En - a*Lz + En*r^2

        return begin
            -((2*((-a)*En + Lz)^2*r^2*swsh_piover2*(-(((2*a4Enterm*(-2 + 2*r)^2)/_Delta^3 - (2*a4Enterm)/_Delta^2 - (2*(-2 + 2*r)*(-4*a*Lz + 3*a^2*En*r + 8*En*r^3 + a^2*En*(2 + 3*r)))/_Delta^2 + (6*a^2*En + 24*En*r^2)/_Delta)/((-a^2)*En + a4Enterm/_Delta)^2) + 
            (2*(-((a4Enterm*(-2 + 2*r))/_Delta^2) + (-4*a*Lz + 3*a^2*En*r + 8*En*r^3 + a^2*En*(2 + 3*r))/_Delta)^2)/((-a^2)*En + a4Enterm/_Delta)^3)*(En + (2*a2Enterm)/(r*_Delta)))/(Pi*sqrtterm)) - 
            (4*((-a)*En + Lz)^2*swsh_piover2*(-((a4Enterm*(-2 + 2*r))/_Delta^2) + (-4*a*Lz + 3*a^2*En*r + 8*En*r^3 + a^2*En*(2 + 3*r))/_Delta)*(-((r^2*(-((2*a2Enterm*(-2 + 2*r))/(r*_Delta^2)) + (4*En)/_Delta - (2*a2Enterm)/(r^2*_Delta)))/sqrtterm) - (2*r*(En + (2*a2Enterm)/(r*_Delta)))/sqrtterm + 
            (r^2*(En + (2*a2Enterm)/(r*_Delta))*(-((-2 + 2*r)*(((-a)*En + Lz)^2 + r^2)) - 4*En*r*(a*Lz - En*(a^2 + r^2)) - 2*r*_Delta))/(2*sqrtterm^3) + (I*m*r^2*(En + (2*a2Enterm)/(r*_Delta))*(a*(a*Lz - 2*En*r) - Lz*_Delta))/(_Delta*((a*Lz - En*(a^2 + r^2))^2 - (((-a)*En + Lz)^2 + r^2)*_Delta)) + 
            (I*r^4*(En + (2*a2Enterm)/(r*_Delta))^2*omega)/((a*Lz - En*(a^2 + r^2))^2 - (((-a)*En + Lz)^2 + r^2)*_Delta)))/(Pi*((-a^2)*En + a4Enterm/_Delta)^2) + 
            (1/(Pi*((-a^2)*En + a4Enterm/_Delta)))*(2*((-a)*En + Lz)^2*swsh_piover2*(-((4*r*(-((2*a2Enterm*(-2 + 2*r))/(r*_Delta^2)) + (4*En)/_Delta - (2*a2Enterm)/(r^2*_Delta)) + 2*(En + (2*a2Enterm)/(r*_Delta)) + r^2*((a2Enterm*(4*r*((2*(-2 + 2*r)^2)/_Delta^3 - 2/_Delta^2) - (8*(-2 + 2*r))/_Delta^2))/(2*r^2) + 
            2*(-(a2Enterm/r^3) + En/r)*(-((4*r*(-2 + 2*r))/_Delta^2) + 4/_Delta) + (4*((3*a2Enterm)/r^4 - (3*En)/r^2)*r)/_Delta))/sqrtterm) + (m^2*r^2*(En + (2*a2Enterm)/(r*_Delta))*(a*(a*Lz - 2*En*r) - Lz*_Delta)^2)/(sqrtterm^3*_Delta^2) + 
            ((r^2*(-((2*a2Enterm*(-2 + 2*r))/(r*_Delta^2)) + (4*En)/_Delta - (2*a2Enterm)/(r^2*_Delta)) + 2*r*(En + (2*a2Enterm)/(r*_Delta)))*(-((-2 + 2*r)*(((-a)*En + Lz)^2 + r^2)) - 4*En*r*(a*Lz - En*(a^2 + r^2)) - 2*r*_Delta))/sqrtterm^3 - 
            (2*I*m*(a*(a*Lz - 2*En*r) - Lz*_Delta)*(-((r^2*(-((2*a2Enterm*(-2 + 2*r))/(r*_Delta^2)) + (4*En)/_Delta - (2*a2Enterm)/(r^2*_Delta)))/sqrtterm) - (2*r*(En + (2*a2Enterm)/(r*_Delta)))/sqrtterm + (r^2*(En + (2*a2Enterm)/(r*_Delta))*(-((-2 + 2*r)*(((-a)*En + Lz)^2 + r^2)) - 4*En*r*(a*Lz - En*(a^2 + r^2)) - 2*r*_Delta))/(2*sqrtterm^3)))/(sqrtterm*_Delta) + 
            (I*m*r^2*(En + (2*a2Enterm)/(r*_Delta))*((-2*a*En - Lz*(-2 + 2*r))/(sqrtterm*_Delta) - ((-2 + 2*r)*(a*(a*Lz - 2*En*r) - Lz*_Delta))/(sqrtterm*_Delta^2) - ((a*(a*Lz - 2*En*r) - Lz*_Delta)*(-((-2 + 2*r)*(((-a)*En + Lz)^2 + r^2)) - 4*En*r*(a*Lz - En*(a^2 + r^2)) - 2*r*_Delta))/(2*sqrtterm^3*_Delta)))/sqrtterm + 
            r^2*(En + (2*a2Enterm)/(r*_Delta))*((8*En^2*r^2 - 4*r*(-2 + 2*r) - 2*(((-a)*En + Lz)^2 + r^2) - 4*En*(a*Lz - En*(a^2 + r^2)) - 2*_Delta)/(2*sqrtterm^3) - (3*(-((-2 + 2*r)*(((-a)*En + Lz)^2 + r^2)) - 4*En*r*(a*Lz - En*(a^2 + r^2)) - 2*r*_Delta)^2)/(4*((a*Lz - En*(a^2 + r^2))^2 - (((-a)*En + Lz)^2 + r^2)*_Delta)^(5/2))) + 
            (2*m*r^4*(En + (2*a2Enterm)/(r*_Delta))^2*(a*(a*Lz - 2*En*r) - Lz*_Delta)*omega)/(sqrtterm^3*_Delta) - (3*I*r^2*(En + (2*a2Enterm)/(r*_Delta))*(-((r^2*(-((2*a2Enterm*(-2 + 2*r))/(r*_Delta^2)) + (4*En)/_Delta - (2*a2Enterm)/(r^2*_Delta)))/sqrtterm) - (2*r*(En + (2*a2Enterm)/(r*_Delta)))/sqrtterm + 
            (r^2*(En + (2*a2Enterm)/(r*_Delta))*(-((-2 + 2*r)*(((-a)*En + Lz)^2 + r^2)) - 4*En*r*(a*Lz - En*(a^2 + r^2)) - 2*r*_Delta))/(2*sqrtterm^3))*omega)/sqrtterm + (r^6*(En + (2*a2Enterm)/(r*_Delta))^3*omega^2)/sqrtterm^3)) + 
            (-((r^2*(-((2*a2Enterm*(-2 + 2*r))/(r*_Delta^2)) + (4*En)/_Delta - (2*a2Enterm)/(r^2*_Delta)))/sqrtterm) - (2*r*(En + (2*a2Enterm)/(r*_Delta)))/sqrtterm + (r^2*(En + (2*a2Enterm)/(r*_Delta))*(-((-2 + 2*r)*(((-a)*En + Lz)^2 + r^2)) - 4*En*r*(a*Lz - En*(a^2 + r^2)) - 2*r*_Delta))/(2*sqrtterm^3) + 
            (I*m*r^2*(En + (2*a2Enterm)/(r*_Delta))*(a*(a*Lz - 2*En*r) - Lz*_Delta))/(_Delta*((a*Lz - En*(a^2 + r^2))^2 - (((-a)*En + Lz)^2 + r^2)*_Delta)) + (I*r^4*(En + (2*a2Enterm)/(r*_Delta))^2*omega)/((a*Lz - En*(a^2 + r^2))^2 - (((-a)*En + Lz)^2 + r^2)*_Delta))*
            ((2*I*(-2*a*En + 2*Lz)*(a2Enterm + sqrtterm)*(psptheta_piover2 + swsh_piover2*(m - a*omega)))/(Pi*(a4Enterm - a^2*En*_Delta)) - (4*((-a)*En + Lz)^2*swsh_piover2*(-(1/r) + (I*((-a)*m + (a^2 + r^2)*omega))/(a^2 - 2*r + r^2)))/(Pi*((-a^2)*En + a4Enterm/_Delta))) - 
            (r^2*(En + (2*a2Enterm)/(r*_Delta))*(-((2*I*(-2*a*En + 2*Lz)*(-4*a*Lz + 3*a^2*En*r + 8*En*r^3 - a^2*En*(-2 + 2*r) + a^2*En*(2 + 3*r))*(a2Enterm + sqrtterm)*(psptheta_piover2 + swsh_piover2*(m - a*omega)))/(Pi*(a4Enterm - a^2*En*_Delta)^2)) + 
            (2*I*(-2*a*En + 2*Lz)*(2*En*r + (-((-2 + 2*r)*(((-a)*En + Lz)^2 + r^2)) - 4*En*r*(a*Lz - En*(a^2 + r^2)) - 2*r*_Delta)/(2*sqrtterm))*(psptheta_piover2 + swsh_piover2*(m - a*omega)))/(Pi*(a4Enterm - a^2*En*_Delta)) - 
            (4*((-a)*En + Lz)^2*swsh_piover2*(1/r^2 + (2*I*r*omega)/(a^2 - 2*r + r^2) - (I*(-2 + 2*r)*((-a)*m + (a^2 + r^2)*omega))/(a^2 - 2*r + r^2)^2))/(Pi*((-a^2)*En + a4Enterm/_Delta)) + 
            (4*((-a)*En + Lz)^2*swsh_piover2*(-((a4Enterm*(-2 + 2*r))/_Delta^2) + (-4*a*Lz + 3*a^2*En*r + 8*En*r^3 + a^2*En*(2 + 3*r))/_Delta)*(-(1/r) + (I*((-a)*m + (a^2 + r^2)*omega))/(a^2 - 2*r + r^2)))/(Pi*((-a^2)*En + a4Enterm/_Delta)^2)))/sqrtterm - 
            (r^2*(En + (2*a2Enterm)/(r*_Delta))*((2*I*(-2*a*En + 2*Lz)*(a2Enterm + sqrtterm)*(psptheta_piover2 + swsh_piover2*(m - a*omega))*(2/r - (I*((-a)*m + (a^2 + r^2)*omega))/(a^2 - 2*r + r^2)))/(Pi*(a4Enterm - a^2*En*_Delta)) + 
            (2*(a2Enterm + sqrtterm)^2*(r^4*((6*I*a*psptheta_piover2)/r^4 - p2sptheta2_piover2/r^3 + (12*a^2*swsh_piover2)/r^5) + 4*I*a*r^3*(-(psptheta_piover2/r^3) + (3*I*a*swsh_piover2)/r^4) - psptheta_piover2*r*(m - a*omega) + r^4*(-(psptheta_piover2/r^3) + (3*I*a*swsh_piover2)/r^4)*(m - a*omega) - r*swsh_piover2*(m - a*omega)^2 + swsh_piover2*(2*r - I*a*(m - a*omega))))/(Pi*r*_Delta*(a4Enterm - a^2*En*_Delta)) + 
            (2*((-a)*En + Lz)^2*swsh_piover2*(-((2*I*((-a)*m + (a^2 + r^2)*omega))/(r*(a^2 - 2*r + r^2))) - ((-a)*m + (a^2 + r^2)*omega)^2/(a^2 - 2*r + r^2)^2 + I*((2*r*omega)/(a^2 - 2*r + r^2) - ((-2 + 2*r)*((-a)*m + (a^2 + r^2)*omega))/(a^2 - 2*r + r^2)^2)))/(Pi*((-a^2)*En + a4Enterm/_Delta))))/sqrtterm
        end
    elseif s == -2
        # s = -2
        return nothing
    else
        # Throw an error, this spin weight is not supported
        throw(DomainError(s, "Currently only spin weight s of +/-2 are supported"))
    end
end