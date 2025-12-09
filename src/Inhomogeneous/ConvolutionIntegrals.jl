module ConvolutionIntegrals

using LinearAlgebra
using SpinWeightedSpheroidalHarmonics
using KerrGeodesics

using ..GeneralizedSasakiNakamura: BoundaryCondition, IN, UP
using ..Coordinates: rstar_from_r
using ..SolutionsY: Y_solution
using ..GridSampling

export convolution_integral_trapezoidal, convolution_integral_levin

function horizon_factor(ω, a, m)
    rp = 1 + sqrt(1 - a^2)
    κ = ω - a * m / (2 * rp)
    ϵ0 = sqrt(1 - a^2) / (4 * rp)
    return ω / (κ * (2 * rp)^3 * (κ^2 + 4 * ϵ0^2) * 64pi)
end

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
    KG_samp = kerr_geo_generic_sample(KG, N_sample, K_sample)
    carter_samp = carter_ingredients_sample(KG_samp, a, m, omega)
    carter_factor = trapezoidal_1d_integral(carter_samp)

    if s == 2
        SH_p2 = spin_weighted_spheroidal_harmonic(s, l, m, a * omega)
        Yup_soln = Y_solution(s, l, m, a, omega, UP; rsin = rsin, rsout = rsout)
        Yup_samp = y_sample_p2(Yup_soln, KG_samp)
        SH_p2_samp = swsh_sample(SH_p2, KG_samp)
        inte_in_samp_p2 = integrand_generic_sample_p2(KG_samp, Yup_samp, SH_p2_samp, n, k)
        integral_p2 = trapezoidal_2d_integral(inte_in_samp_p2)
        hf = horizon_factor(omega, a, m)
        return Dict(
            "Amplitude" => integral_p2,
            "omega" => omega,
            "EnergyFlux" => hf * abs2(integral_p2),
            "AngularMomentumFlux" => hf * m * abs2(integral_p2) / omega,
            "CarterConstantFlux" => 2 * hf * abs2(integral_p2) * (carter_factor + k * ϒθ) / omega,
            "Trajectory" => KG,
            "YSolution" => Yup_soln,
            "SWSH" => SH_p2
        )
    elseif s == -2
        SH_m2 = spin_weighted_spheroidal_harmonic(s, l, m, a * omega)
        Yin_soln = Y_solution(s, l, m, a, omega, IN; rsin = rsin, rsout = rsout)
        Yin_samp = y_sample_m2(Yin_soln, KG_samp)
        SH_m2_samp = swsh_sample(SH_m2, KG_samp)
        inte_in_samp_m2 = integrand_generic_sample_m2(KG_samp, Yin_samp, SH_m2_samp, n, k)
        integral_m2 = trapezoidal_2d_integral(inte_in_samp_m2)
        return Dict(
            "Amplitude" => integral_m2,
            "omega" => omega,
            "EnergyFlux" => abs2(integral_m2)/(4.0pi*omega^2),
            "AngularMomentumFlux" => m * abs2(integral_m2) / (4.0pi * omega^3),
            "CarterConstantFlux" => abs2(integral_m2) * (carter_factor + k * ϒθ) / (2.0pi*omega^3),
            "Trajectory" => KG,
            "YSolution" => Yin_soln,
            "SWSH" => SH_m2
        )
    else
        error("Spin weight s must be either 2 or -2.")
    end
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
    KG_samp = kerr_geo_generic_sample_cheby(KG, N_sample, K_sample)
    KG_trap = kerr_geo_generic_sample(KG, N_sample, K_sample)
    carter_samp = carter_ingredients_sample(KG_trap, a, m, omega)
    carter_factor = trapezoidal_1d_integral(carter_samp)

    if s == 2
        SH_p2 = spin_weighted_spheroidal_harmonic(s, l, m, a * omega)
        SH_p2_samp = swsh_sample(SH_p2, KG_samp)
        Yup_soln = Y_solution(s, l, m, a, omega, UP; rsin = rsin, rsout = rsout)
        Yup_samp = y_sample_p2(Yup_soln, KG_samp)
        Jpp_grid_up, Jpm_grid_up, Jmp_grid_up, Jmm_grid_up, drphase_up, dθphase_up, rphaseL_up, rphaseR_up, θphaseL_up, θphaseR_up, prefactor_up = integrand_generic_sample_cheby_p2(KG_samp, Yup_samp, SH_p2_samp, n, k)
        integralpp_up = levin_2d_integral(Jpp_grid_up, drphase_up, dθphase_up, rphaseL_up, rphaseR_up, θphaseL_up, θphaseR_up)
        integralpm_up = levin_2d_integral(Jpm_grid_up, drphase_up, - dθphase_up, rphaseL_up, rphaseR_up, - θphaseL_up, - θphaseR_up)
        integralmp_up = levin_2d_integral(Jmp_grid_up, - drphase_up, dθphase_up, - rphaseL_up, - rphaseR_up, θphaseL_up, θphaseR_up)
        integralmm_up = levin_2d_integral(Jmm_grid_up, - drphase_up, - dθphase_up, - rphaseL_up, - rphaseR_up, - θphaseL_up, - θphaseR_up)
        integral_up = (integralpp_up + integralpm_up + integralmp_up + integralmm_up) * prefactor_up
        
        hf = horizon_factor(omega, a, m)
        return Dict(
            "Amplitude" => integral_up,
            "omega" => omega,
            "EnergyFlux" => hf * abs2(integral_up),
            "AngularMomentumFlux" => hf * m * abs2(integral_up) / omega,
            "CarterConstantFlux" => 2 * hf * abs2(integral_up) * (carter_factor + k * ϒθ) / omega,
            "Trajectory" => KG,
            "YSolution" => Yup_soln,
            "SWSH" => SH_p2
        )
    elseif s == -2
        SH_m2 = spin_weighted_spheroidal_harmonic(s, l, m, a * omega)
        Yin_soln = Y_solution(s, l, m, a, omega, IN; rsin = rsin, rsout = rsout)
        Yin_samp = y_sample_m2(Yin_soln, KG_samp)
        SH_m2_samp = swsh_sample(SH_m2, KG_samp)
        Jpp_grid_in, Jpm_grid_in, Jmp_grid_in, Jmm_grid_in, drphase_in, dθphase_in, rphaseL_in, rphaseR_in, θphaseL_in, θphaseR_in, prefactor_in = integrand_generic_sample_cheby_m2(KG_samp, Yin_samp, SH_m2_samp, n, k)
        integralpp_in = levin_2d_integral(Jpp_grid_in, drphase_in, dθphase_in, rphaseL_in, rphaseR_in, θphaseL_in, θphaseR_in)
        integralpm_in = levin_2d_integral(Jpm_grid_in, drphase_in, - dθphase_in, rphaseL_in, rphaseR_in, - θphaseL_in, - θphaseR_in)
        integralmp_in = levin_2d_integral(Jmp_grid_in, - drphase_in, dθphase_in, - rphaseL_in, - rphaseR_in, θphaseL_in, θphaseR_in)
        integralmm_in = levin_2d_integral(Jmm_grid_in, - drphase_in, - dθphase_in, - rphaseL_in, - rphaseR_in, - θphaseL_in, - θphaseR_in)
        integral_in = (integralpp_in + integralpm_in + integralmp_in + integralmm_in) * prefactor_in
        return Dict(
            "Amplitude" => integral_in,
            "omega" => omega,
            "EnergyFlux" => abs2(integral_in)/(4.0pi*omega^2),
            "AngularMomentumFlux" => m * abs2(integral_in) / (4.0pi * omega^3),
            "CarterConstantFlux" => abs2(integral_in) * (carter_factor + k * ϒθ) / (2.0pi*omega^3),
            "Trajectory" => KG,
            "YSolution" => Yin_soln,
            "SWSH" => SH_m2
        )
    else
        error("Spin weight s must be either 2 or -2.")
    end
end

function convolution_integral_eccentric_trapezoidal(a, p, e, s, l, m, n, N_sample)
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
    KG_samp = kerr_geo_eccentric_sample(KG, N_sample)  # Equatorial sampling (θ fixed)
    
    if s == 2
        Yup_soln = Y_solution(s, l, m, a, omega, UP; rsin = rsin, rsout = rsout)
        SH_p2 = spin_weighted_spheroidal_harmonic(s, l, m, a * omega)
        Yup_samp = y_sample_p2(Yup_soln, KG_samp)
        S0_p2 = SH_p2(π/2, 0.0)      
        S1_p2 = SH_p2(π/2, 0.0; theta_derivative=1)
        S2_p2 = (m^2 + (a * omega)^2 - 2 - 2 * a * omega * m - SH_p2.lambda) * S0_p2
        inte_samp_p2 = integrand_eccentric_sample_p2(KG_samp, Yup_samp, (S0_p2, S1_p2, S2_p2), n)
        integral_p2 = trapezoidal_1d_integral(inte_samp_p2)
        hf = horizon_factor(omega, a, m)
        return Dict(
            "Amplitude" => integral_p2,
            "omega" => omega,
            "EnergyFlux" => hf * abs2(integral_p2),
            "AngularMomentumFlux" => hf * m * abs2(integral_p2) / omega,
            "CarterConstantFlux" => 0.0,
            "Trajectory" => KG,
            "YSolution" => Yup_soln,
            "SWSH" => SH_p2
        )
    elseif s == -2
        Yin_soln = Y_solution(s, l, m, a, omega, IN; rsin = rsin, rsout = rsout)
        SH_m2 = spin_weighted_spheroidal_harmonic(s, l, m, a * omega)  # θ-fixed harmonics
        Yin_samp = y_sample_m2(Yin_soln, KG_samp) 
        S0_m2 = SH_m2(π/2, 0.0)      
        S1_m2 = SH_m2(π/2, 0.0; theta_derivative=1)
        S2_m2 = (m^2 + (a * omega)^2 + 2 - 2 * a * omega * m - SH_m2.lambda) * S0_m2
        inte_samp_m2 = integrand_eccentric_sample_m2(KG_samp, Yin_samp, (S0_m2, S1_m2, S2_m2), n)
        integral_m2 = trapezoidal_1d_integral(inte_samp_m2)
        return Dict(
            "Amplitudef" => integral_m2,
            "omega" => omega,
            "EnergyFlux" => abs2(integral_m2)/(4.0pi*omega^2),
            "AngularMomentumFlux" => m * abs2(integral_m2) / (4.0pi * omega^3),
            "CarterConstantFlux" => 0.0,
            "Trajectory" => KG,
            "YSolution" => Yin_soln,
            "SWSH" => SH_m2
        )
    else
        error("Spin weight s must be either 2 or -2.")
    end
end

function convolution_integral_eccentric_levin(a, p, e, s, l, m, n, N_sample)

    KG = Kerr_Geodesics(a, p, e, 1.0)

    Frequencies = KG["Frequencies"]
    Γ = Frequencies["ϒt"]
    ϒr = Frequencies["ϒr"]
    ϒφ = Frequencies["ϒϕ"]
    omega = (m * ϒφ + n * ϒr) / Γ
    rsin = min(-50.0, 50.0 * log10(1 - a))
    rsout = max(500.0, 10pi / abs(omega))
    KG_samp = kerr_geo_eccentric_sample_cheby(KG, N_sample)

    if s == 2
        Yup_soln = Y_solution(s, l, m, a, omega, UP; rsin = rsin, rsout = rsout)
        SH_p2 = spin_weighted_spheroidal_harmonic(s, l, m, a * omega)
        Yup_samp = y_sample_p2(Yup_soln, KG_samp)
        S0_p2 = SH_p2(π/2, 0.0)      
        S1_p2 = SH_p2(π/2, 0.0; theta_derivative=1)
        S2_p2 = (m^2 + (a * omega)^2 - 2 - 2 * a * omega * m - SH_p2.lambda) * S0_p2
        Jp_up, Jm_up, drphase_up, rphaseL_up, rphaseR_up, prefactor_up = integrand_eccentric_sample_cheby_p2(KG_samp, Yup_samp, (S0_p2, S1_p2, S2_p2), n)
        integralp_p2 = levin_1d_integral(Jp_up, drphase_up, rphaseL_up, rphaseR_up)
        integralm_p2 = levin_1d_integral(Jm_up, - drphase_up, - rphaseL_up, - rphaseR_up)
        integral_p2 = (integralp_p2 + integralm_p2) * prefactor_up
        hf = horizon_factor(omega, a, m)
        return Dict(
            "Amplitude" => integral_p2,
            "omega" => omega,
            "EnergyFlux_hor" => hf * abs2(integral_p2),
            "AngularMomentumFlux_hor" => hf * m * abs2(integral_p2) / omega,
            "CarterConstantFlux_hor" => 0.0,
            "Trajectory" => KG,
            "YSolution" => Yup_soln,
            "SWSH" => SH_p2
        )
    elseif s == -2
        Yin_soln = Y_solution(s, l, m, a, omega, IN; rsin = rsin, rsout = rsout)
        SH_m2 = spin_weighted_spheroidal_harmonic(s, l, m, a * omega)
        Yin_samp = y_sample_m2(Yin_soln, KG_samp)   
        S0_m2 = SH_m2(π/2, 0.0)      
        S1_m2 = SH_m2(π/2, 0.0; theta_derivative=1)
        S2_m2 = (m^2 + (a * omega)^2 + 2 - 2 * a * omega * m - SH_m2.lambda) * S0_m2
        Jp_in, Jm_in, drphase_in, rphaseL_in, rphaseR_in, prefactor_in = integrand_eccentric_sample_cheby_m2(KG_samp, Yin_samp, (S0_m2, S1_m2, S2_m2), n)
        integralp_m2 = levin_1d_integral(Jp_in, drphase_in, rphaseL_in, rphaseR_in)
        integralm_m2 = levin_1d_integral(Jm_in, - drphase_in, - rphaseL_in, - rphaseR_in)
        integral_m2 = (integralp_m2 + integralm_m2) * prefactor_in
        return Dict(
            "Amplitude" => integral_m2,
            "omega" => omega,
            "EnergyFlux" => abs2(integral_m2)/(4.0pi*omega^2),
            "AngularMomentumFlux" => m * abs2(integral_m2) / (4.0pi * omega^3),
            "CarterConstantFlux" => 0.0,
            "Trajectory" => KG,
            "YSolution" => Yin_soln,
            "SWSH" => SH_m2
        )
    else
        error("Spin weight s must be either 2 or -2.")
    end
end

function convolution_integral_inclined_trapezoidal(a, p, x, s, l, m, k, K_sample)
    KG = Kerr_Geodesics(a, p, 0.0, x)  # Assume this handles radial-fixed case
    
    # Extract frequencies and compute omega (polar-only mode)
    Frequencies = KG["Frequencies"]
    Γ = Frequencies["ϒt"]
    ϒθ = Frequencies["ϒθ"]
    ϒφ = Frequencies["ϒϕ"]
    omega = (m * ϒφ + k * ϒθ) / Γ  # No radial contribution (n omitted)
    
    # Radial solution parameters (r fixed, but retain for consistency)
    rsin = min(-30.0, 30.0 * log10(1 - a))
    rsout = max(200.0, 10pi / abs(omega))
    KG_samp = kerr_geo_inclined_sample(KG, K_sample)
    carter_samp = carter_ingredients_sample(KG_samp, a, m, omega)
    carter_factor = trapezoidal_1d_integral(carter_samp)

    if s == 2
        Yup_soln = Y_solution(s, l, m, a, omega, UP; rsin = rsin, rsout = rsout)
        SH_p2 = spin_weighted_spheroidal_harmonic(s, l, m, a * omega)
        Yup = Dict("params" => (s=2, l=l, m=m, a=a, omega=omega, lambda=SH_p2.lambda),
            "Cinc" => Yup_soln.asymptotic.Cinc,
            "Y" => Yup_soln.solution.Y_inf(rstar_from_r(a, p)),
            "Yp" => Yup_soln.solution.Yp_inf(rstar_from_r(a, p)),
            "X" => Yup_soln.solution.X(rstar_from_r(a, p)))
        SH_p2_samp = swsh_sample(SH_p2, KG_samp)
        inte_samp_up = integrand_inclined_sample_p2(KG_samp, Yup, SH_p2_samp, k)
        integral_up = trapezoidal_1d_integral(inte_samp_up)
        hf = horizon_factor(omega, a, m)
        return Dict(
            "Amplitude" => integral_up,
            "omega" => omega,
            "EnergyFlux" => hf * abs2(integral_up),
            "AngularMomentumFlux" => hf * m * abs2(integral_up) / omega,
            "CarterConstantFlux" => 2 * hf * abs2(integral_up) * (carter_factor + k * ϒθ) / omega,
            "Trajectory" => KG,
            "YSolution" => Yup_soln,
            "SWSH" => SH_p2
        )
    elseif s == -2
        Yin_soln = Y_solution(s, l, m, a, omega, IN; rsin = rsin, rsout = rsout)
        SH_m2 = spin_weighted_spheroidal_harmonic(s, l, m, a * omega) 
        Yin = Dict("params" => (s=-2, l=l, m=m, a=a, omega=omega, lambda=SH_m2.lambda),
            "Binc" => Yin_soln.asymptotic.Binc,
            "Y" => Yin_soln.solution.Y_hor(rstar_from_r(a, p)),
            "Yp" => Yin_soln.solution.Yp_hor(rstar_from_r(a, p)),
            "X" => Yin_soln.solution.X(rstar_from_r(a, p)))  
        SH_m2_samp = swsh_sample(SH_m2, KG_samp)
        inte_samp_in = integrand_inclined_sample_m2(KG_samp, Yin, SH_m2_samp, k)
        integral_in = trapezoidal_1d_integral(inte_samp_in)
        return Dict(
            "Amplitude" => integral_in,
            "omega" => omega,
            "EnergyFlux" => abs2(integral_in)/(4.0pi*omega^2),
            "AngularMomentumFlux" => m * abs2(integral_in) / (4.0pi * omega^3),
            "CarterConstantFlux" => abs2(integral_in) * (carter_factor + k * ϒθ) / (2.0pi*omega^3),
            "Trajectory" => KG,
            "YSolution" => Yin_soln,
            "SWSH" => SH_m2
        )
    else
        error("Spin weight s must be either 2 or -2.")
    end
end

function convolution_integral_inclined_levin(a, p, x, s, l, m, k, K_sample)

    KG = Kerr_Geodesics(a, p, 0.0, x)

    Frequencies = KG["Frequencies"]
    Γ = Frequencies["ϒt"]
    ϒθ = Frequencies["ϒθ"]
    ϒφ = Frequencies["ϒϕ"]
    omega = (m * ϒφ + k * ϒθ) / Γ
    rsin = min(-30.0, 30.0 * log10(1 - a))
    rsout = max(200.0, 10pi / abs(omega))
    KG_samp = kerr_geo_inclined_sample_cheby(KG, K_sample)
    KG_trap = kerr_geo_inclined_sample(KG, K_sample)
    carter_samp = carter_ingredients_sample(KG_trap, a, m, omega)
    carter_factor = trapezoidal_1d_integral(carter_samp)

    if s == 2
        Yup_soln = Y_solution(s, l, m, a, omega, UP; rsin = rsin, rsout = rsout)
        SH_p2 = spin_weighted_spheroidal_harmonic(s, l, m, a * omega)
        Yup = Dict("params" => (s=2, l=l, m=m, a=a, omega=omega, lambda=SH_p2.lambda),
            "Cinc" => Yup_soln.asymptotic.Cinc,
            "Y" => Yup_soln.solution.Y_inf(rstar_from_r(a, p)),
            "Yp" => Yup_soln.solution.Yp_inf(rstar_from_r(a, p)),
            "X" => Yup_soln.solution.X(rstar_from_r(a, p)))
        SH_p2_samp = swsh_sample(SH_p2, KG_samp)
        Jp_up, Jm_up, dθphase_up, θphaseL_up, θphaseR_up, prefactor_up = integrand_inclined_sample_cheby_p2(KG_samp, Yup, SH_p2_samp, k)
        integralp_up = levin_1d_integral(Jp_up, dθphase_up, θphaseL_up, θphaseR_up)
        integralm_up = levin_1d_integral(Jm_up, - dθphase_up, - θphaseL_up, - θphaseR_up)
        integral_up = (integralp_up + integralm_up) * prefactor_up
        hf = horizon_factor(omega, a, m)
        return Dict(
            "Amplitude" => integral_up,
            "omega" => omega,
            "EnergyFlux" => hf * abs2(integral_up),
            "AngularMomentumFlux" => hf * m * abs2(integral_up) / omega,
            "CarterConstantFlux" => 2 * hf * abs2(integral_up) * (carter_factor + k * ϒθ) / omega,
            "Trajectory" => KG,
            "YSolution" => Yup_soln,
            "SWSH" => SH_p2
        )
    elseif s == -2
        Yin_soln = Y_solution(s, l, m, a, omega, IN; rsin = rsin, rsout = rsout)
        SH_m2 = spin_weighted_spheroidal_harmonic(s, l, m, a * omega)
        Yin = Dict("params" => (s=-2, l=l, m=m, a=a, omega=omega, lambda=SH_m2.lambda),
            "Binc" => Yin_soln.asymptotic.Binc,
            "Y" => Yin_soln.solution.Y_hor(rstar_from_r(a, p)),
            "Yp" => Yin_soln.solution.Yp_hor(rstar_from_r(a, p)),
            "X" => Yin_soln.solution.X(rstar_from_r(a, p)))
        SH_m2_samp = swsh_sample(SH_m2, KG_samp)
        Jp_in, Jm_in, dθphase_in, θphaseL_in, θphaseR_in, prefactor_in = integrand_inclined_sample_cheby_m2(KG_samp, Yin, SH_m2_samp, k)
        integralp_in = levin_1d_integral(Jp_in, dθphase_in, θphaseL_in, θphaseR_in)
        integralm_in = levin_1d_integral(Jm_in, - dθphase_in, - θphaseL_in, - θphaseR_in)
        integral_in = (integralp_in + integralm_in) * prefactor_in
        return Dict(
            "Amplitude" => integral_in,
            "omega" => omega,
            "EnergyFlux" => abs2(integral_in)/(4.0pi*omega^2),
            "AngularMomentumFlux" => m * abs2(integral_in) / (4.0pi * omega^3),
            "CarterConstantFlux" => abs2(integral_in) * (carter_factor + k * ϒθ) / (2.0pi*omega^3),
            "Trajectory" => KG,
            "YSolution" => Yin_soln,
            "SWSH" => SH_m2
        )
    else
        error("Spin weight s must be either 2 or -2.")
    end
end

function convolution_integral_circular_equatorial_m2(a, p, l, m)

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
    Y_soln = Y_solution(-2, l, m, a, ω, IN; rsin=rsin, rsout=rsout)
    Y = Y_soln.solution.Y_hor(rstar_from_r(a, r))
    Yp = Y_soln.solution.Yp_hor(rstar_from_r(a, r))
    X = Y_soln.solution.X(rstar_from_r(a, r))
    Binc = Y_soln.asymptotic.Binc 

    c = a * ω
    SH = spin_weighted_spheroidal_harmonic(-2, l, m, c)
    S0 = SH(θ, 0.0)
    S1 = SH(θ, 0.0; theta_derivative=1)
    λ = SH.lambda
    S2 = (m^2 + c^2 + 2 - 2c*m - λ) * S0

    ρ = -1.0 / r
    ρbar = -1.0 / r

    sinθ = sin(θ)
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

function convolution_integral_circular_equatorial_p2(a, p, l, m)

    KG = Kerr_Geodesics(a, p, 0.0, 1.0)
    E = KG["Energy"]
    Lz = KG["AngularMomentum"]
    Γ = KG["Frequencies"]["ϒt"]
    ϒφ = KG["Frequencies"]["ϒϕ"]
    ω = m * ϒφ / Γ 

    r = p
    θ = π/2
    r2 = r^2
    Δ = r2 - 2*r + a^2 
    rp = 1.0 + sqrt(1.0 - a^2)
    rm = 1.0 - sqrt(1.0 - a^2)
    rsin = min(-50.0, 50.0 * log10(1 - a))
    rsout = max(500.0, 10π / abs(ω))

    Y_soln = Y_solution(2, l, m, a, ω, UP; rsin=rsin, rsout=rsout)
    Y = Y_soln.solution.Y_inf(rstar_from_r(a, r))
    Yp = Y_soln.solution.Yp_inf(rstar_from_r(a, r))
    X = Y_soln.solution.X(rstar_from_r(a, r))
    Cinc = Y_soln.asymptotic.Cinc

    c = a * ω
    SH = spin_weighted_spheroidal_harmonic(2, l, m, c)
    S0 = SH(θ, 0.0)
    S1 = SH(θ, 0.0; theta_derivative=1)
    λ = SH.lambda
    S2 = (m^2 + c^2 - 2 - 2c*m - λ) * S0

    ρ = -1.0 / r
    ρbar = -1.0 / r

    sinθ = sin(θ)
    L1 = m - c
    L2 = m - c
    L2p = -2.0

    L2S = S1 + L2 * S0
    L1Sp = S2 + L1 * S1
    L1L2S = L1Sp + L2p * S0 + L2 * S1 + L1 * L2 * S0
    L1pL2pS = -r * L1L2S + 3im * a * L1 * S0 + 2im * a * S1 -1im * a * L2 * S0

    rs = rstar_from_r(a, r)
    log_term = log((r - rp)/(r - rm))
    phase = exp(-(im*ω*rs - im*a*m*log_term/(2*sqrt(1 - a^2))))

    Wll = 2.0 * ρbar * L1pL2pS * r2 * Y * phase
    term_A = L2S + im*a*sinθ*(ρ - ρbar)*S0 
    term_B = 2.0 * Y + r*Yp
    term_C = L2S*(ρ + ρbar)*r*Y
    Wlm = 4.0 * (r * ρbar * (term_A * term_B + term_C) / ρ) * phase
    Wmm = 2.0 * S0 * ρbar * (X / sqrt(r^2 + a^2) + 2 * phase * (Y + 2 * r * Y * ρ + r * Yp * (2 + r * ρ))) / ρ

    L = (E * (r2 + a^2) - a * Lz) / Δ 
    M = - im * sinθ * (a*E - Lz / (sinθ^2))

    J = Wll*L^2 + Wlm*L*M + Wmm*M^2

    κ = ω - a * m / (2 * rp)
    c0 = 24 + 12im * ω + λ * (10 + λ) - 12 * a * ω * (a * ω - m)
    c1 = -32im * a * m - 8im * a * m * λ + 8im * a^2 * ω * (1 + λ)
    c2 = 12 * a^2 - 24im * a * m - 24 * a^2 * m^2 + 24im * a^2 * ω + 48 * a^3 * m * ω - 24 * a^4 * ω^2
    c3 = - 24im * a^3 * (a * ω - m) - 24 * a^2
    c4 = 12 * a^4
    η = c0 + c1 / rp + c2 / rp^2 + c3 / rp^3 + c4 / rp^4
    factor = 2*sqrt(2)*rp^(3/2)*(((4*ω)*(im-4*ω)-a*m*(im-8*ω)-a^2*(m^2+2im*ω-4*ω^2))*rp^2
            +a^2*(im-4*ω)*(a*m-2*ω)*rp)/(2*rp^3*(24+10*λ+λ^2+12im*ω)-rp^2*(8im*a*m*(11
            +2*λ+6im*ω)+a^2*(24+24*m^2+10*λ+λ^2-28im*ω-16im*λ*ω+48*ω^2))+8im*a^3*rp*(m
            *(7+λ-6im*ω)-a*ω*(4+λ))+12*a^5*ω*(a*ω-3*m))
    ϵ0 = sqrt(1 - a^2) / (4 * rp)
    energy_factor = ω / (κ * (2 * rp)^3 * (κ^2 + 4 * ϵ0^2) * 64pi)

    integral = - im * π * factor * J * η / (κ * Cinc * Γ)

    return Dict(
        "Amplitude" => integral,
        "omega" => ω,
        "EnergyFlux" => abs2(integral) * energy_factor,
        "AngularMomentumFlux" => abs2(integral) * energy_factor * m / ω,
        "CarterConstantFlux" => 0.0,
        "Trajectory" => KG,
        "YSolution" => Y_soln,
        "SWSH" => SH
    )
end

function convolution_integral_circular_equatorial(a, p, s, l, m)
    if s == 2
        return convolution_integral_circular_equatorial_p2(a, p, l, m)
    elseif s == -2
        return convolution_integral_circular_equatorial_m2(a, p, l, m)
    else
        error("Spin weight s must be either 2 or -2.")
    end
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
        return convolution_integral_inclined_levin(a, p, x, s, l, m, k, K)
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
        return convolution_integral_eccentric_levin(a, p, e, s, l, m, n, N)
    else
        return convolution_integral_generic_levin(a, p, e, x, s, l, m, n, k, N, K)
    end
end

end
