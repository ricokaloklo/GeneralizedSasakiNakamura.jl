module ConvolutionIntegrals

using LinearAlgebra
using SpinWeightedSpheroidalHarmonics
using KerrGeodesics

using ..GeneralizedSasakiNakamura: BoundaryCondition, IN, UP
using ..Coordinates: rstar_from_r
using ..SolutionsY: Y_solution
using ..GridSampling

export convolution_integral_trapezoidal, convolution_integral_levin

function trapezoidal_1d_integral(f::Vector{<:Number})
    # Get grid size: M = number of sampling points (1D)
    M = length(f)
    M < 2 && return 0.0  # Need at least 2 points for integration
    
    # Step size: integration range [0, π] for 1D phase (qr or qθ)
    Δq = π / (M - 1)  # Uniform step between adjacent points
    
    # --------------------------
    # Assign weights based on position
    # --------------------------
    # Initialize weight vector with all 1.0 (internal points default)
    weight = ones(Float64, M)
    
    # Edge points (first and last) have weight 0.5 (trapezoidal rule)
    weight[1] = 0.5    # First point
    weight[end] = 0.5  # Last point
    
    # --------------------------
    # Compute 1D trapezoidal integral
    # --------------------------
    # Element-wise product of function values and weights
    weighted_f = f .* weight
    # Sum all weighted values and multiply by step size
    integral_value = sum(weighted_f) * Δq
    
    return integral_value
end

function trapezoidal_2d_integral(f::Matrix{<:Number})
    # Get grid size: N×K (N=radial points, K=polar points)
    N, K = size(f)
    
    # Calculate step sizes: integration range [0, π] for both qr and qθ
    Δqr = π / (N - 1)  # Step in radial direction (qr)
    Δqθ = π / (K - 1)  # Step in polar direction (qθ)
    
    # Initialize weight matrix with all 1.0 (internal points default)
    weight = ones(Float64, N, K)
    
    # --------------------------
    # Assign weights based on position
    # --------------------------
    # 1. Vertex points (4 corners): weight = 1/4
    weight[1, 1] = 0.25    # (1,1): top-left vertex
    weight[1, K] = 0.25    # (1,K): top-right vertex
    weight[N, 1] = 0.25    # (N,1): bottom-left vertex
    weight[N, K] = 0.25    # (N,K): bottom-right vertex
    
    # 2. Edge points (non-vertex): weight = 1/2
    # Top edge (row=1, columns 2~K-1)
    weight[1, 2:K-1] .= 0.5
    # Bottom edge (row=N, columns 2~K-1)
    weight[N, 2:K-1] .= 0.5
    # Left edge (column=1, rows 2~N-1)
    weight[2:N-1, 1] .= 0.5
    # Right edge (column=K, rows 2~N-1)
    weight[2:N-1, K] .= 0.5
    
    # --------------------------
    # Compute 2D trapezoidal integral
    # --------------------------
    # Element-wise product of function values and weights
    weighted_f = f .* weight
    # Sum all weighted values, multiply by step product (Δqr*Δqθ)
    integral_value = sum(weighted_f) * Δqr * Δqθ
    
    return integral_value
end

function cheb_nodes_and_D(n::Int)
    # Chebyshev extremal nodes on reference interval [-1, 1]
    x_ref = [cos(π * k / (n - 1)) for k in 0:n-1]
    
    # Construct Chebyshev differentiation matrix
    D_ref = zeros(Float64, n, n)
    for i in 1:n, j in 1:n
        if i != j
            # Weight factors (2 for endpoints, 1 otherwise)
            c_i = (i == 1 || i == n) ? 2.0 : 1.0
            c_j = (j == 1 || j == n) ? 2.0 : 1.0
            # Off-diagonal elements (standard Chebyshev differentiation formula)
            D_ref[i,j] = (-1)^(i+j) * c_i / (c_j * (x_ref[i] - x_ref[j]))
        end
    end
    # Diagonal elements enforce row sum = 0 (property of differentiation matrices)
    for i in 1:n
        D_ref[i,i] = -sum(D_ref[i,:])
    end
    return x_ref, D_ref
end

function levin_1d_integral(f_vals::Vector{<:Number}, gprime_vals::Vector{<:Number}, g0::Number, gπ::Number)
    k = length(f_vals)
    @assert length(gprime_vals) == k "gprime length must match f ($(length(gprime_vals)) ≠ $k)"
    
    # Fixed integration interval [0, π]
    a, b = 0.0, π
    scale = 2.0 / (b - a)  # Scaling factor for differentiation matrix
    
    # Generate Chebyshev nodes (reference [-1,1] and physical [0,π]) and differentiation matrix
    x_ref, D_ref = cheb_nodes_and_D(k)
    x_phys = (b - a) * (x_ref .+ 1.0) ./ 2.0 .+ a  # Map [-1,1] → [0,π]
    D = scale * D_ref  # Physical-space differentiation matrix
    
    # Prepare complex-valued inputs for numerical stability
    rhs = convert(Vector{ComplexF64}, f_vals)
    gprime_hat = convert(Vector{ComplexF64}, gprime_vals)
    
    # Construct system matrix: A = D + i·diag(g') (core of Levin method)
    A = D .+ 1im .* Diagonal(gprime_hat)
    
    # Solve linear system using pivoted QR decomposition (numerically robust)
    F = qr(A, Val(true))  # Pivoted QR factorization
    R = F.R               # Upper triangular factor
    piv = F.p             # Pivot indices
    diagR = abs.(diag(R)) # Magnitudes of diagonal elements (singularity measure)
    Anorm = maximum(diagR)# Estimate of matrix norm
    tol_sing = Anorm * 1e-12  # Singularity tolerance
    
    # Compute solution p(x) at Chebyshev nodes
    pvals = zeros(ComplexF64, k)
    if !all(diagR .< tol_sing)
        l = count(≥(tol_sing), diagR)  # Number of significant singular values
        y = F.Q[:, 1:l]' * rhs          # Project RHS onto significant subspace
        ytop = R[1:l, 1:l] \ y          # Solve reduced system
        colp = zeros(ComplexF64, k)
        colp[piv[1:l]] = ytop           # Map solution back to original indices
        pvals .= colp
    end
    
    # Identify boundary indices using physical coordinates (most reliable method)
    idx_x0 = findmin(x_phys)[2]  # Index of x=0 in physical nodes
    idx_xπ = findmax(x_phys)[2]  # Index of x=π in physical nodes
    
    # Core integral result from Levin's method: p(b)e^{ig(b)} - p(a)e^{ig(a)}
    term_b = pvals[idx_xπ] * exp(1im * gπ)
    term_a = pvals[idx_x0] * exp(1im * g0)
    return term_b - term_a
end

function levin_2d_integral(f_vals::Matrix{<:Number}, g1_prime_vals::Vector{<:Number}, g2_prime_vals::Vector{<:Number}, g1_0::Number, g1_π::Number, g2_0::Number, g2_π::Number)
    # Get grid dimensions (n×k, non-square allowed)
    n, k = size(f_vals)
    
    # Validate input dimensions
    @assert length(g1_prime_vals) == n "g1_prime_vals length must match x-dimension n ($(length(g1_prime_vals)) ≠ $n)"
    @assert length(g2_prime_vals) == k "g2_prime_vals length must match y-dimension k ($(length(g2_prime_vals)) ≠ $k)"
    
    # ----------------------
    # Step 1: Solve y-direction Levin equation for each x (n rows)
    # ----------------------
    # Generate y-direction Chebyshev nodes and differentiation matrix
    y_ref, D_ref_y = cheb_nodes_and_D(k)
    y_phys = π .* (y_ref .+ 1.0) ./ 2.0  # y ∈ [0,π]
    scale_y = 2.0 / π  # Scale factor for y-direction
    D_y = scale_y * D_ref_y  # y-differentiation matrix
    
    # Storage for y-solutions: p(x,y) has size n×k
    p_vals = zeros(ComplexF64, n, k)
    
    for i in 1:n  # Iterate over x-nodes (each row of f_vals)
        # Extract f(xᵢ, y) for fixed xᵢ (i-th row)
        f_row = f_vals[i, :]
        rhs_y = convert(Vector{ComplexF64}, f_row)
        
        # Construct y-matrix: Aᵧ = Dᵧ + i·diag(∂g₂/∂y)
        A_y = D_y .+ 1im .* Diagonal(convert(Vector{ComplexF64}, g2_prime_vals))
        
        # Solve Aᵧ·p = f using pivoted QR
        F_y = qr(A_y, Val(true))
        R_y = F_y.R
        piv_y = F_y.p
        diagR_y = abs.(diag(R_y))
        Anorm_y = maximum(diagR_y)
        tol_sing_y = Anorm_y * 1e-16
        
        p_row = zeros(ComplexF64, k)
        if !all(diagR_y .< tol_sing_y)
            l_y = count(≥(tol_sing_y), diagR_y)
            y_proj = F_y.Q[:, 1:l_y]' * rhs_y
            y_sol = R_y[1:l_y, 1:l_y] \ y_proj
            colp = zeros(ComplexF64, k)
            colp[piv_y[1:l_y]] = y_sol
            p_row .= colp
        end
        p_vals[i, :] = p_row
    end
    
    # ----------------------
    # Step 2: Reduce to x-direction problem using y-boundaries
    # ----------------------
    # Find y-boundary indices (y=0 and y=π)
    idx_y0 = findmin(y_phys)[2]  # Index of y=0 in y_phys
    idx_yπ = findmax(y_phys)[2]  # Index of y=π in y_phys
    
    # Compute h(x) = p(x, π)e^{ig₂(π)} - p(x, 0)e^{ig₂(0)} (length n)
    h_vals = [p_vals[i, idx_yπ] * exp(1im * g2_π) - p_vals[i, idx_y0] * exp(1im * g2_0) for i in 1:n]
    
    # ----------------------
    # Step 3: Solve x-direction Levin equation
    # ----------------------
    # Generate x-direction Chebyshev nodes and differentiation matrix
    x_ref, D_ref_x = cheb_nodes_and_D(n)
    x_phys = π .* (x_ref .+ 1.0) ./ 2.0  # x ∈ [0,π]
    scale_x = 2.0 / π  # Scale factor for x-direction
    D_x = scale_x * D_ref_x  # x-differentiation matrix
    
    # Construct x-matrix: Aₓ = Dₓ + i·diag(∂g₁/∂x)
    A_x = D_x .+ 1im .* Diagonal(convert(Vector{ComplexF64}, g1_prime_vals))
    rhs_x = convert(Vector{ComplexF64}, h_vals)
    
    # Solve Aₓ·q = h using pivoted QR
    F_x = qr(A_x, Val(true))
    R_x = F_x.R
    piv_x = F_x.p
    diagR_x = abs.(diag(R_x))
    Anorm_x = maximum(diagR_x)
    tol_sing_x = Anorm_x * 1e-16
    
    q_vals = zeros(ComplexF64, n)
    if !all(diagR_x .< tol_sing_x)
        l_x = count(≥(tol_sing_x), diagR_x)
        x_proj = F_x.Q[:, 1:l_x]' * rhs_x
        x_sol = R_x[1:l_x, 1:l_x] \ x_proj
        colp_x = zeros(ComplexF64, n)
        colp_x[piv_x[1:l_x]] = x_sol
        q_vals .= colp_x
    end
    
    # ----------------------
    # Step 4: Compute final 2D integral from x-boundaries
    # ----------------------
    idx_x0 = findmin(x_phys)[2]  # Index of x=0 in x_phys
    idx_xπ = findmax(x_phys)[2]  # Index of x=π in x_phys
    integral_val = q_vals[idx_xπ] * exp(1im * g1_π) - q_vals[idx_x0] * exp(1im * g1_0)
    
    return integral_val
end

function convolution_integral_generic_trapezoidal(a, p, e, x, s, l, m, n, k, N_sample, K_sample)

    KG = Kerr_Geodesics(a, p, e, x)

    Frequencies = KG["Frequencies"]
    Γ = Frequencies["ϒt"]
    ϒr = Frequencies["ϒr"]
    ϒθ = Frequencies["ϒθ"]
    ϒφ = Frequencies["ϒϕ"]
    omega = (m * ϒφ + n * ϒr + k * ϒθ) / Γ
    rsin = min(-50.0, 50.0 * log10(1 - a))
    rsout = max(500.0, 10pi / abs(omega))

    SH = spin_weighted_spheroidal_harmonic(s, l, m, a * omega)
    λ = SH.lambda
    Yin_soln = Y_solution(s, l, m, a, omega, IN; rsin = rsin, rsout = rsout)
    Binc = Yin_soln.asymptotic.Binc

    KG_samp = kerr_geo_generic_sample(KG, N_sample, K_sample)
    Yin_samp = y_sample(Yin_soln, KG_samp)
    SH_samp = swsh_sample(SH, KG_samp)

    inte_in_samp = integrand_generic_sample(KG_samp, Yin_samp, SH_samp, n, k)
    integral = trapezoidal_2d_integral(inte_in_samp) / Binc
    carter_samp = carter_ingredients_sample(KG_samp, a, m, omega)
    carter_factor = trapezoidal_1d_integral(carter_samp)

    return Dict(
        "Amplitude" => integral,
        "omega" => omega,
        "EnergyFlux" => abs2(integral)/(4.0pi*omega^2),
        "AngularMomentumFlux" => m * abs2(integral) / (4.0pi * omega^3),
        "CarterConstantFlux" => abs2(integral) * (carter_factor + k * ϒθ) / (2.0pi*omega^3),
        "Trajectory" => KG,
        "YSolution" => Yin_soln,
        "SWSH" => SH
    )
end

function convolution_integral_generic_levin(a, p, e, x, s, l, m, n, k, N_sample, K_sample)

    KG = Kerr_Geodesics(a, p, e, x)

    Frequencies = KG["Frequencies"]
    Γ = Frequencies["ϒt"]
    ϒr = Frequencies["ϒr"]
    ϒθ = Frequencies["ϒθ"]
    ϒφ = Frequencies["ϒϕ"]
    omega = (m * ϒφ + n * ϒr + k * ϒθ) / Γ
    rsin = min(-50.0, 50.0 * log10(1 - a))
    rsout = max(500.0, 10pi / abs(omega))

    SH = spin_weighted_spheroidal_harmonic(s, l, m, a * omega)
    λ = SH.lambda
    Yin_soln = Y_solution(s, l, m, a, omega, IN; rsin = rsin, rsout = rsout)
    Binc = Yin_soln.asymptotic.Binc

    KG_samp = kerr_geo_generic_sample_cheby(KG, N_sample, K_sample)
    Yin_samp = y_sample(Yin_soln, KG_samp)
    SH_samp = swsh_sample(SH, KG_samp)

    Jpp_grid_in, Jpm_grid_in, Jmp_grid_in, Jmm_grid_in, drphase_in, dθphase_in, rphaseL_in, rphaseR_in, θphaseL_in, θphaseR_in, prefactor_in = integrand_generic_sample_cheby(KG_samp, Yin_samp, SH_samp, n, k)

    integralpp_in = levin_2d_integral(Jpp_grid_in, drphase_in, dθphase_in, rphaseL_in, rphaseR_in, θphaseL_in, θphaseR_in)
    integralpm_in = levin_2d_integral(Jpm_grid_in, drphase_in, - dθphase_in, rphaseL_in, rphaseR_in, - θphaseL_in, - θphaseR_in)
    integralmp_in = levin_2d_integral(Jmp_grid_in, - drphase_in, dθphase_in, - rphaseL_in, - rphaseR_in, θphaseL_in, θphaseR_in)
    integralmm_in = levin_2d_integral(Jmm_grid_in, - drphase_in, - dθphase_in, - rphaseL_in, - rphaseR_in, - θphaseL_in, - θphaseR_in)
    integral = (integralpp_in + integralpm_in + integralmp_in + integralmm_in) * prefactor_in / Binc
    KG_trap = kerr_geo_generic_sample(KG, N_sample, K_sample)
    carter_samp = carter_ingredients_sample(KG_trap, a, m, omega)
    carter_factor = trapezoidal_1d_integral(carter_samp)

    return Dict(
        "Amplitude" => integral,
        "omega" => omega,
        "EnergyFlux" => abs2(integral)/(4.0pi*omega^2),
        "AngularMomentumFlux" => m * abs2(integral) / (4.0pi * omega^3),
        "CarterConstantFlux" => abs2(integral) * (carter_factor + k * ϒθ) / (2.0pi*omega^3),
        "Trajectory" => KG,
        "YSolution" => Yin_soln,
        "SWSH" => SH
    )
end

function convolution_integral_eccentric_trapezoidal(a, p, e, s, l, m, n, N_sample)
    # --------------------------
    # 1. Generate geodesic data (equatorial: θ fixed)
    # --------------------------
    KG = Kerr_Geodesics(a, p, e, 1.0)  # Assume this handles equatorial case
    
    # Extract frequencies and compute omega (radial-only mode)
    Frequencies = KG["Frequencies"]
    Γ = Frequencies["ϒt"]
    ϒr = Frequencies["ϒr"]
    ϒφ = Frequencies["ϒϕ"]
    omega = (m * ϒφ + n * ϒr) / Γ  # No θ contribution (k omitted)
    
    # Radial solution parameters
    rsin = min(-50.0, 50.0 * log10(1 - a))
    rsout = max(500.0, 10pi / abs(omega))
    
    # --------------------------
    # 2. Solve radial function and spherical harmonics
    # --------------------------
    Y_soln = Y_solution(s, l, m, a, omega, IN; rsin = rsin, rsout = rsout)
    Binc = Y_soln.asymptotic.Binc
    SH = spin_weighted_spheroidal_harmonic(s, l, m, a * omega)  # θ-fixed harmonics
    
    # --------------------------
    # 3. 1D sampling (radial only)
    # --------------------------
    KG_samp = kerr_geo_eccentric_sample(KG, N_sample)  # Equatorial sampling (θ fixed)
    Y_samp = y_sample(Y_soln, KG_samp)                 # Radial function sampling 
    S0 = SH(π/2, 0.0)      
    S1 = SH(π/2, 0.0; theta_derivative=1)
    S2 = (m^2 + (a * omega)^2 - s - 2 * a * omega * m - SH.lambda) * S0

    # --------------------------
    # 4. 1D integrand and integral
    # --------------------------
    inte_samp = integrand_eccentric_sample(KG_samp, Y_samp, (S0, S1, S2), n)  # No k
    integral = trapezoidal_1d_integral(inte_samp) / Binc                        # 1D trapezoidal
    
    return Dict(
        "Amplitude" => integral,
        "omega" => omega,
        "EnergyFlux" => abs2(integral)/(4.0pi*omega^2),
        "AngularMomentumFlux" => m * abs2(integral) / (4.0pi * omega^3),
        "CarterConstantFlux" => 0.0,
        "Trajectory" => KG,
        "YSolution" => Y_soln,
        "SWSH" => SH
    )
end

function convolution_integral_eccentric_levin(a, p, e, x, s, l, m, n, N_sample)

    KG = Kerr_Geodesics(a, p, e, x)

    Frequencies = KG["Frequencies"]
    Γ = Frequencies["ϒt"]
    ϒr = Frequencies["ϒr"]
    ϒφ = Frequencies["ϒϕ"]
    omega = (m * ϒφ + n * ϒr) / Γ
    rsin = min(-50.0, 50.0 * log10(1 - a))
    rsout = max(500.0, 10pi / abs(omega))

    Y_soln = Y_solution(s, l, m, a, omega, IN; rsin = rsin, rsout = rsout)
    Binc = Y_soln.asymptotic.Binc
    SH = spin_weighted_spheroidal_harmonic(s, l, m, a * omega)

    KG_samp = kerr_geo_eccentric_sample_cheby(KG, N_sample)
    Y_samp = y_sample(Y_soln, KG_samp)                 # Radial function sampling 
    S0 = SH(π/2, 0.0)      
    S1 = SH(π/2, 0.0; theta_derivative=1)
    S2 = (m^2 + (a * omega)^2 - s - 2 * a * omega * m - SH.lambda) * S0

    Jp, Jm, drphase, rphaseL, rphaseR, prefactor = integrand_eccentric_sample_cheby(KG_samp, Y_samp, (S0, S1, S2), n)
    
    integralp = levin_1d_integral(Jp, drphase, rphaseL, rphaseR)
    integralm = levin_1d_integral(Jm, - drphase, - rphaseL, - rphaseR)
    integral = (integralp + integralm) * prefactor / Binc

    return Dict(
        "Amplitude" => integral,
        "omega" => omega,
        "EnergyFlux" => abs2(integral)/(4.0pi*omega^2),
        "AngularMomentumFlux" => m * abs2(integral) / (4.0pi * omega^3),
        "CarterConstantFlux" => 0.0,
        "Trajectory" => KG,
        "YSolution" => Y_soln,
        "SWSH" => SH
    )
end

function convolution_integral_inclined_trapezoidal(a, p, x, s, l, m, k, K_sample)
    # --------------------------
    # 1. Generate geodesic data (radial fixed: r = p)
    # --------------------------
    KG = Kerr_Geodesics(a, p, 0.0, x)  # Assume this handles radial-fixed case
    
    # Extract frequencies and compute omega (polar-only mode)
    Frequencies = KG["Frequencies"]
    Γ = Frequencies["ϒt"]
    ϒθ = Frequencies["ϒθ"]
    ϒφ = Frequencies["ϒϕ"]
    omega = (m * ϒφ + k * ϒθ) / Γ  # No radial contribution (n omitted)
    
    # Radial solution parameters (r fixed, but retain for consistency)
    rsin = min(-50.0, 50.0 * log10(1 - a))
    rsout = max(500.0, 10pi / abs(omega))
    
    # --------------------------
    # 2. Solve radial function and spherical harmonics
    # --------------------------
    Y_soln = Y_solution(s, l, m, a, omega, 1; rsin = rsin, rsout = rsout)
    Binc = Y_soln.asymptotic.Binc
    SH = spin_weighted_spheroidal_harmonic(s, l, m, a * omega)  # Polar-dependent harmonics
    Y = Dict("params" => (s=s, l=l, m=m, a=a, omega=omega, lambda=SH.lambda),
            "asymptotic_amplitudes" => Y_soln.asymptotic.Binc,
            "Y" => Y_soln.solution.Y_hor(rstar_from_r(a, p)),
            "Yp" => Y_soln.solution.Yp_hor(rstar_from_r(a, p)),
            "X" => Y_soln.solution.X(rstar_from_r(a, p)))
    # --------------------------
    # 3. 1D sampling (polar only)
    # --------------------------
    KG_samp = kerr_geo_inclined_sample(KG, K_sample)  # Radial fixed sampling (r = p)              # Radial function (scalar for fixed r)
    SH_samp = swsh_sample(SH, KG_samp)                # Polar harmonics sampling
    
    # --------------------------
    # 4. 1D integrand and integral
    # --------------------------
    inte_samp = integrand_inclined_sample(KG_samp, Y, SH_samp, k)  # No n
    integral = trapezoidal_1d_integral(inte_samp) / Binc                        # 1D trapezoidal
    carter_samp = carter_ingredients_sample(KG_samp, a, m, omega)
    carter_factor = trapezoidal_1d_integral(carter_samp)
    
    return Dict(
        "Amplitude" => integral,
        "omega" => omega,
        "EnergyFlux" => abs2(integral)/(4.0pi*omega^2),
        "AngularMomentumFlux" => m * abs2(integral) / (4.0pi * omega^3),
        "CarterConstantFlux" => abs2(integral) * (carter_factor + k * ϒθ) / (2.0pi*omega^3),
        "Trajectory" => KG,
        "YSolution" => Y_soln,
        "SWSH" => SH
    )
end

function convolution_integral_inclined_levin(a, p, e, x, s, l, m, k, K_sample)

    KG = Kerr_Geodesics(a, p, e, x)

    Frequencies = KG["Frequencies"]
    Γ = Frequencies["ϒt"]
    ϒθ = Frequencies["ϒθ"]
    ϒφ = Frequencies["ϒϕ"]
    omega = (m * ϒφ + k * ϒθ) / Γ
    rsin = min(-50.0, 50.0 * log10(1 - a))
    rsout = max(500.0, 10pi / abs(omega))

    Y_soln = Y_solution(s, l, m, a, omega, IN; rsin = rsin, rsout = rsout)
    Binc = Y_soln.asymptotic.Binc
    SH = spin_weighted_spheroidal_harmonic(s, l, m, a * omega)

    KG_samp = kerr_geo_inclined_sample_cheby(KG, K_sample)
    Y = Dict("params" => (s=s, l=l, m=m, a=a, omega=omega, lambda=SH.lambda),
            "asymptotic_amplitudes" => Y_soln.asymptotic.Binc,
            "Y" => Y_soln.solution.Y_hor(rstar_from_r(a, p)),
            "Yp" => Y_soln.solution.Yp_hor(rstar_from_r(a, p)),
            "X" => Y_soln.solution.X(rstar_from_r(a, p)))
    SH_samp = swsh_sample(SH, KG_samp)

    Jp, Jm, dθphase, θphaseL, θphaseR, prefactor = integrand_inclined_sample_cheby(KG_samp, Y, SH_samp, k)
    
    integralp = levin_1d_integral(Jp, dθphase, θphaseL, θphaseR)
    integralm = levin_1d_integral(Jm, - dθphase, - θphaseL, - θphaseR)
    integral = (integralp + integralm) * prefactor / Binc
    KG_trap = kerr_geo_inclined_sample(KG, K_sample)
    carter_samp = carter_ingredients_sample(KG_trap, a, m, omega)
    carter_factor = trapezoidal_1d_integral(carter_samp)

    return Dict(
        "Amplitude" => integral,
        "omega" => omega,
        "EnergyFlux" => abs2(integral)/(4.0pi*omega^2),
        "AngularMomentumFlux" => m * abs2(integral) / (4.0pi * omega^3),
        "CarterConstantFlux" => abs2(integral) * (carter_factor + k * ϒθ) / (2.0pi*omega^3),
        "Trajectory" => KG,
        "YSolution" => Y_soln,
        "SWSH" => SH
    )
end

function convolution_integral_circular_equatorial(a, p, s, l, m)

    KG = Kerr_Geodesics(a, p, 0.0, 1.0)
    E = KG["Energy"]
    Lz = KG["AngularMomentum"]
    Γ = KG["Frequencies"]["ϒt"]
    ϒφ = KG["Frequencies"]["ϒϕ"]
    ω = m * ϒφ / Γ 

    r = p
    θ = π/2
    r2 = r^2
    Δ = r2 - 2r + a^2 
    rp = 1.0 + sqrt(1.0 - a^2)
    rm = 1.0 - sqrt(1.0 - a^2)

    rsin = min(-50.0, 50.0 * log10(1 - a))
    rsout = max(500.0, 10π / abs(ω))
    Y_soln = Y_solution(s, l, m, a, ω, IN; rsin=rsin, rsout=rsout)
    Y = Y_soln.solution.Y_hor(rstar_from_r(a, r))
    Yp = Y_soln.solution.Yp_hor(rstar_from_r(a, r))
    X = Y_soln.solution.X(rstar_from_r(a, r))
    Binc = Y_soln.asymptotic.Binc 

    c = a * ω
    SH = spin_weighted_spheroidal_harmonic(s, l, m, c)
    S0 = SH(θ, 0.0)
    S1 = SH(θ, 0.0; theta_derivative=1)
    λ = SH.lambda
    S2 = (m^2 + c^2 - s - 2c*m - λ) * S0

    ρ = -1.0 / r
    ρbar = -1.0 / r

    sinθ = sin(θ)
    inv_sinθ = 1/sinθ
    cotθ = cot(θ)
    L1 = -m + c
    L2 = -m + c
    L2p = -2.0

    L2S = S1 + L2 * S0
    L1Sp = S2 + L1 * S1
    L1L2S = L1Sp + L2p * S0 + L2 * S1 + L1 * L2 * S0
    L1pL2pS = -r * L1L2S + 3im * a * L1 * S0 + 2im * a * S1 -1im * a * L2 * S0

    rs = rstar_from_r(a, r)
    log_term = log((r - rp)/(r - rm))
    phase = exp(im*ω*rs - im*a*m*log_term/(2*sqrt(1 - a^2)))

    Wnn = ρbar * L1pL2pS * r2 * Y * phase / 2.0

    term_A = L2S + im*a*sinθ*(ρ - ρbar)*S0 
    term_B = 2Y + r*Yp
    term_C = L2S*(ρ + ρbar)*r*Y
    Wnmbar = (-r * ρbar * (term_A * term_B + term_C) / ρ) * phase
    
    term_D = X / (2*sqrt(r2 + a^2))
    term_E = (Y + 2r*Yp) * phase
    term_F = ρ*r*(2Y + r*Yp) * phase
    Wmbarmbar = S0 * ρbar * (term_D + term_E + term_F) / ρ

    numerator_N = E * (r2 + a^2) - a * Lz 
    N = numerator_N / Δ 
    term_M = im * sinθ * (a*E - Lz / (sinθ^2))
    Mbar = term_M 

    J = Wnn*N^2 + Wnmbar*N*Mbar + Wmbarmbar*Mbar^2
    integral = 4im * π * ω * J / (Binc * Γ)

    return Dict(
        "Amplitude" => integral,
        "omega" => ω,
        "EnergyFlux" => abs2(integral)/(4.0pi*ω^2),
        "AngularMomentumFlux" => m * abs2(integral) / (4.0pi * ω^3),
        "CarterConstantFlux" => 0.0,
        "Trajectory" => KG,
        "YSolution" => Y_soln,
        "SWSH" => SH
    )
end

function convolution_integral_trapezoidal(a, p, e, x, s, l, m, n, k; N = 256, K = 256)
    KG = Kerr_Geodesics(a, p, e, x)
    if typeof(KG) == Vector{String}
        return KG
    end

    if isapprox(e, 0.0; atol=1e-12) && isapprox(abs(x), 1.0; atol=1e-12)
        if n != 0 || k != 0
            return Dict(
                    "Amplitude" => 0.0 + 0.0im,
                    "omega" => nothing,
                    "EnergyFlux" => 0.0,
                    "AngularMomentumFlux" => 0.0,
                    "CarterConstantFlux" => 0.0,
                    "Trajectory" => KG,
                    "YSolution" => nothing,
                    "SWSH" => nothing
                    )
        end
        return convolution_integral_circular_equatorial(a, p, s, l, m)
    elseif isapprox(e, 0.0; atol=1e-12) && !isapprox(abs(x), 1.0; atol=1e-12)
        if n != 0
            return Dict(
                    "Amplitude" => 0.0 + 0.0im,
                    "omega" => nothing,
                    "EnergyFlux" => 0.0,
                    "AngularMomentumFlux" => 0.0,
                    "CarterConstantFlux" => 0.0,
                    "Trajectory" => KG,
                    "YSolution" => nothing,
                    "SWSH" => nothing
                    )
        end
        return convolution_integral_inclined_trapezoidal(a, p, x, s, l, m, k, K)
    elseif !isapprox(e, 0.0; atol=1e-12) && isapprox(abs(x), 1.0; atol=1e-12)
        if k != 0
            return Dict(
                    "Amplitude" => 0.0 + 0.0im,
                    "omega" => nothing,
                    "EnergyFlux" => 0.0,
                    "AngularMomentumFlux" => 0.0,
                    "CarterConstantFlux" => 0.0,
                    "Trajectory" => KG,
                    "YSolution" => nothing,
                    "SWSH" => nothing
                    )
        end
        return convolution_integral_eccentric_trapezoidal(a, p, e, s, l, m, n, N)
    else
        return convolution_integral_generic_trapezoidal(a, p, e, x, s, l, m, n, k, N, K)
    end
end

function convolution_integral_levin(a, p, e, x, s, l, m, n, k; N = 256, K = 32)
    KG = Kerr_Geodesics(a, p, e, x)
    if typeof(KG) == Vector{String}
        return KG
    end

    if isapprox(e, 0.0; atol=1e-12) && isapprox(abs(x), 1.0; atol=1e-12)
        if n != 0 || k != 0
            return Dict(
                    "Amplitude" => 0.0 + 0.0im,
                    "omega" => nothing,
                    "EnergyFlux" => 0.0,
                    "AngularMomentumFlux" => 0.0,
                    "CarterConstantFlux" => 0.0,
                    "Trajectory" => KG,
                    "YSolution" => nothing,
                    "SWSH" => nothing
                    )
        end
        return convolution_integral_circular_equatorial(a, p, s, l, m)
    elseif isapprox(e, 0.0; atol=1e-12) && !isapprox(abs(x), 1.0; atol=1e-12)
        if n != 0
            return Dict(
                    "Amplitude" => 0.0 + 0.0im,
                    "omega" => nothing,
                    "EnergyFlux" => 0.0,
                    "AngularMomentumFlux" => 0.0,
                    "CarterConstantFlux" => 0.0,
                    "Trajectory" => KG,
                    "YSolution" => nothing,
                    "SWSH" => nothing
                    )
        end
        return convolution_integral_inclined_levin(a, p, e, x, s, l, m, k, K)
    elseif !isapprox(e, 0.0; atol=1e-12) && isapprox(abs(x), 1.0; atol=1e-12)
        if k != 0
            return Dict(
                    "Amplitude" => 0.0 + 0.0im,
                    "omega" => nothing,
                    "EnergyFlux" => 0.0,
                    "AngularMomentumFlux" => 0.0,
                    "CarterConstantFlux" => 0.0,
                    "Trajectory" => KG,
                    "YSolution" => nothing,
                    "SWSH" => nothing
                    )
        end
        return convolution_integral_eccentric_levin(a, p, e, x, s, l, m, n, N)
    else
        return convolution_integral_generic_levin(a, p, e, x, s, l, m, n, k, N, K)
    end
end

end