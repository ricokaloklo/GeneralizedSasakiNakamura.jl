using ForwardDiff
using StaticArrays
using DifferentialEquations
using ..GeneralizedSasakiNakamura: Potentials, InitialConditions,
    AsymptoticExpansionCoefficients, ExtremalRadial, Solutions

struct QNMFrequencyTangent end
_frequency_dual(z, dz=one(z)) = complex(
    ForwardDiff.Dual{QNMFrequencyTangent}(real(z),real(dz)),
    ForwardDiff.Dual{QNMFrequencyTangent}(imag(z),imag(dz)))
_primal(z::Complex) = complex(ForwardDiff.value(real(z)),ForwardDiff.value(imag(z)))
_primal(z::Real) = ForwardDiff.value(z)
_tangent_component(x::Real) = zero(x)
_tangent_component(x::ForwardDiff.Dual) = ForwardDiff.partials(x)[1]
_tangent(z::Complex) = complex(_tangent_component(real(z)),_tangent_component(imag(z)))
_tangent(z::Real) = _tangent_component(z)

# The angular spectral matrix is quadratic in c. Its symmetric (not Hermitian)
# eigenpair therefore gives an analytic derivative on the selected simple branch.
function _qnm_lambda_tangent(mode, a, omega)
    pair=angular_branch(mode,a,omega)
    v=pair.mixing
    c=a*omega
    n=length(v)
    Mplus=_big_angular_spectral_matrix(c+1,mode.s,mode.m,n)
    Mminus=_big_angular_spectral_matrix(c-1,mode.s,mode.m,n)
    dA=sum(v.*((Mplus-Mminus)*v))/(2sum(v.*v))
    return pair.lambda,a*dA+2a*c-2mode.m*a
end

function _qnm_complex_rstar(a,r)
    if abs(a)==1
        h=r-1
        return 1+h-log(4)+2*(log(h)-1/h)
    end
    k=sqrt(1-a*a)
    rp,rm=1+k,1-k
    return r+rp/k*log((r-rp)/2)-rm/k*log((r-rm)/2)
end

# Request the whole recurrence once; coefficient_sequence otherwise rebuilds
# the Taylor expansion at every requested order.
function _qnm_endpoint_coefficients(f,s,m,a,w,lambda,order)
    workspace=Ref{Any}(nothing)
    f(s,m,a,w,lambda,order;data_type=typeof(w),workspace)
    return workspace[].expansion_coeffs
end

function _qnm_infinity_jet(s,m,a,w,lambda,r,order,sign)
    A=AsymptoticExpansionCoefficients
    f=sign==1 ? A.outgoing_coefficient_at_inf : A.ingoing_coefficient_at_inf
    coeff=_qnm_endpoint_coefficients(f,s,m,a,w,lambda,order)
    # An isolated zero coefficient is not convergence (at a=0, c[3]=0 but
    # c[4] is nonzero). Use the same fixed polynomial for value and tangent.
    value=sum(coeff[n+1]/(w*r)^n for n in 0:order)
    deriv=sum(-n*coeff[n+1]/(w^n*r^(n+1)) for n in 1:order)
    q=(r*r-2r+a*a)/(r*r+a*a)
    phase=exp(im*sign*w*_qnm_complex_rstar(a,r))
    X=phase*value
    V=phase*(deriv+im*sign*w/q*value)
    return SVector(_primal(X),_primal(V),_tangent(X),_tangent(V))
end

function _qnm_horizon_jet(s,m,a,w,lambda,r,order)
    q=(r*r-2r+a*a)/(r*r+a*a)
    if abs(a)==1
        # Same irregular-horizon normalization as ExtremalRadial, whose
        # leading conversion coefficient is available analytically.
        raw=ExtremalRadial._raw_extremal_gsn_horizon_state(s,m,a,w,lambda,r,:ingoing;order)
        scale=inv(ConversionFactors.Btrans(s,m,a,w,lambda))
        X,V=raw[1]/scale,raw[2]/(q*scale)
    else
        coeff=_qnm_endpoint_coefficients(AsymptoticExpansionCoefficients.ingoing_coefficient_at_hor,
            s,m,a,w,lambda,order)
        value=InitialConditions.gansatz(i->coeff[i+1],a,r;order)
        deriv=InitialConditions.dgansatz_dr(i->coeff[i+1],a,r;order)
        p=w-m*a/(2*(1+sqrt(1-a*a)))
        phase=exp(-im*p*_qnm_complex_rstar(a,r))
        X=phase*value
        V=phase*(deriv-im*p/q*value)
    end
    return SVector(_primal(X),_primal(V),_tangent(X),_tangent(V))
end

@inline function _qnm_tangent_rhs(u,p,t)
    s,m,a,w,lambda,r0,dr=p
    r=r0+t*dr
    # Factor Delta before dividing near a horizon; expanded cancellation can
    # make the adaptive integrator resolve arithmetic noise instead of the ODE.
    kappa=sqrt(1-a*a)
    delta=(r-1-kappa)*(r-1+kappa)
    q=delta/(r*r+a*a)
    dq=((2r-2)*(r*r+a*a)-2r*delta)/(r*r+a*a)^2
    A=Potentials.sF(s,m,a,w,lambda,r)/q-dq/q
    B=Potentials.sU(s,m,a,w,lambda,r)/q^2
    av,bv,ad,bd=_primal(A),_primal(B),_tangent(A),_tangent(B)
    return SVector{length(u)}(ntuple(length(u)) do i
        k=4*((i-1)÷4)
        j=i-k
        dr*(j==1 ? u[k+2] : j==2 ? av*u[k+2]+bv*u[k+1] :
            j==3 ? u[k+4] : av*u[k+4]+bv*u[k+3]+ad*u[k+2]+bd*u[k+1])
    end)
end

function _qnm_tangent_propagate(u,s,m,a,w,lambda,r0,r1,tol; tortoise=false)
    parameters=(s,m,a,w,lambda,r0,r1-r0)
    # The sensitivity is already an explicit state. Do not let SciML also
    # promote that state to dual numbers from an AD-valued parameter tuple.
    rhs=if tortoise
        distance=_qnm_complex_rstar(a,r1)-_qnm_complex_rstar(a,r0)
        function(u,p,t)
            r=u[5]
            dr=(r*r-2r+a*a)/(r*r+a*a)*distance
            state=SVector(u[1],u[2],u[3],u[4])
            v=_qnm_tangent_rhs(state,(s,m,a,w,lambda,r,dr),0.)
            SVector(v[1],v[2],v[3],v[4],dr)
        end
    else
        (u,p,t)->_qnm_tangent_rhs(u,parameters,t)
    end
    tortoise && (u=vcat(u,SVector(r0)))
    prob=ODEProblem(rhs,u,(0.0,1.0))
    sol=solve(prob,Vern9();reltol=tol,abstol=tol,save_everystep=false,dense=false,maxiters=10^7)
    DifferentialEquations.SciMLBase.successful_retcode(sol) || error("tangent propagation: $(sol.retcode)")
    if tortoise
        abs(sol.u[end][5]-r1)<=1e-8*max(1,abs(r1)) ||
            error("horizon tortoise contour reached a different coordinate branch")
        return SVector{4}(sol.u[end][1:4]),sol.destats.nf
    end
    return sol.u[end],sol.destats.nf
end

"""
One frequency: propagate the IN state and its tangent, and the two infinity
bases and their tangents, then differentiate M*c = y at a fixed matching point.
The contour and endpoint coordinates are held fixed while differentiating.
Public dispatch uses a bounded mode family and checks refinement separately.
"""
function _qnm_variational_excitation(root;tol=1e-11,horizon_order=16,
        infinity_order=32,infinity_radius=100.0,match_radius=3.0,horizon_scale=0.15,
        infinity_coordinate=:tortoise,extremal_scale=nothing,
        extremal_radius=0.03,extremal_order=64,horizon_ray=false,
        extremal_precision=53,horizon_tortoise=false)
    mode=root.mode
    s,m,a,w0=mode.s,mode.m,Float64(root.a),ComplexF64(root.omega)
    s in -2:2 || throw(ArgumentError("variational GSN supports integer spin -2 through 2"))
    abs(a)==1 && w0==m*a/2 && throw(DomainError(w0,"synchronous branch-point endpoint"))
    t0=time_ns()
    lambda0,dlambda=_qnm_lambda_tangent(mode,a,w0)
    w,lambda=_frequency_dual(w0),_frequency_dual(lambda0,dlambda)
    angular_s=(time_ns()-t0)/1e9
    rp=1+sqrt(1-a*a)
    rh=if abs(a)==1
        alpha=im*(2w0-m*a)
        distance=extremal_scale === nothing ? extremal_radius :
            min(extremal_radius,abs(alpha)/extremal_scale)
        1+distance*cis(angle(alpha)-copysign(pi/2,imag(alpha)))
    else
        complex(rp+horizon_scale*sqrt(1-a*a))
    end
    rm=complex(match_radius)
    if abs(a)<1 && horizon_tortoise
        p=w0-m*a/(2rp)
        x0=_qnm_complex_rstar(a,rh)
        xm=_qnm_complex_rstar(a,rm)
        target=xm+real(x0-xm)*(1-im*imag(p)/real(p))
        k=sqrt(1-a*a)
        rh=rp+(rh-rp)*cis(imag(target)*k/rp)
        for _ in 1:60
            step=(_qnm_complex_rstar(a,rh)-target)*(rh*rh-2rh+a*a)/(rh*rh+a*a)
            cap=0.5abs(rh-rp)
            abs(step)>cap && (step*=cap/abs(step))
            rh-=step
            abs(step)<1e-14 && break
        end
        abs(_qnm_complex_rstar(a,rh)-target)<1e-8 ||
            error("horizon tortoise endpoint outside the principal coordinate branch")
    end
    radius=max(infinity_radius,(infinity_order+8)/abs(w0))
    ri=rm+radius*sign(real(w0))*cis(-angle(w0))
    if infinity_coordinate === :tortoise
        # Balance the two wave carriers in r*, not r: the logarithmic phase
        # otherwise amplifies tangent roundoff at higher damping.
        target=_qnm_complex_rstar(a,rm)+radius*sign(real(w0))*cis(-angle(w0))
        for _ in 1:12
            step=(_qnm_complex_rstar(a,ri)-target)*(ri*ri-2ri+a*a)/(ri*ri+a*a)
            ri-=step
            abs(step)<=1e-13*abs(ri) && break
        end
        abs(_qnm_complex_rstar(a,ri)-target)<=1e-11*max(1,abs(target)) ||
            error("complex tortoise endpoint did not converge")
    elseif infinity_coordinate !== :r
        throw(ArgumentError("infinity_coordinate must be :r or :tortoise"))
    end
    t0=time_ns()
    h=if abs(a)==1 && extremal_precision>53
        setprecision(BigFloat,extremal_precision) do
            wb=_frequency_dual(Complex{BigFloat}(w0))
            lb=_frequency_dual(Complex{BigFloat}(lambda0),Complex{BigFloat}(dlambda))
            jet=_qnm_horizon_jet(s,m,BigFloat(a),wb,lb,Complex{BigFloat}(rh),extremal_order)
            SVector{4,ComplexF64}(jet)
        end
    else
        _qnm_horizon_jet(s,m,a,w,lambda,rh,abs(a)==1 ? extremal_order : horizon_order)
    end
    incoming=_qnm_infinity_jet(s,m,a,w,lambda,ri,infinity_order,-1)
    outgoing=_qnm_infinity_jet(s,m,a,w,lambda,ri,infinity_order,1)
    endpoint_s=(time_ns()-t0)/1e9
    t0=time_ns()
    hin,nh=if abs(a)==1 && horizon_ray
        turn=1+max(0.5,abs(2w0-m*a))*cis(angle(rh-1))
        first,n1=_qnm_tangent_propagate(h,s,m,a,w,lambda,rh,turn,tol)
        last,n2=_qnm_tangent_propagate(first,s,m,a,w,lambda,turn,rm,tol)
        last,n1+n2
    else
        _qnm_tangent_propagate(h,s,m,a,w,lambda,rh,rm,tol;tortoise=horizon_tortoise)
    end
    both,ni=_qnm_tangent_propagate(vcat(incoming,outgoing),s,m,a,w,lambda,ri,rm,tol)
    propagation_s=(time_ns()-t0)/1e9
    M=SMatrix{2,2}(both[1],both[2],both[5],both[6])
    Mw=SMatrix{2,2}(both[3],both[4],both[7],both[8])
    coefficients=M\SVector(hin[1],hin[2])
    derivative=M\(SVector(hin[3],hin[4])-Mw*coefficients)
    ain,aout=coefficients
    alpha=derivative[1]
    eta(r)=Transformation.eta(s,m,a,w0,lambda0,r)
    q(r)=(r*r-2r+a*a)/(r*r+a*a)
    win=q(ri)*(incoming[1]*outgoing[2]-outgoing[1]*incoming[2])/eta(ri)
    wmatch=q(rm)*det(M)/eta(rm)
    expected=2im*w0/Transformation.eta_coefficient(s,m,a,w0,lambda0,0)
    wronskian_drift=max(abs(wmatch/win-1),abs(win/expected-1))
    return (incidence=ain,reflection=aout,alpha=alpha,B_gsn=aout/(2w0*alpha),
        lambda=lambda0,dlambda=dlambda,wronskian_drift,angular_s,endpoint_s,propagation_s,
        rhs_calls=nh+ni,horizon_rhs_calls=nh,infinity_rhs_calls=ni)
end
