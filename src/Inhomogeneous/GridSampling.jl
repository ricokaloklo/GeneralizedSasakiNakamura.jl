module GridSampling

using ..Coordinates: r_from_rstar, rstar_from_r
using ..GeneralizedSasakiNakamura: IN, UP, Mode
using ..ConversionFactors
using ..ISEM: Y_radial
using ..AsymptoticExpansionCoefficientsY: Ypp_outgoing_inf_m2, Ypp_ingoing_inf_m2, Ypp_outgoing_hor_p2, Ypp_ingoing_hor_p2
using ..InitialConditionsY: Y_outgoing_inf_m2, Yp_outgoing_inf_m2, Y_ingoing_inf_m2, Yp_ingoing_inf_m2, Y_outgoing_hor_p2, Yp_outgoing_hor_p2, Y_ingoing_hor_p2, Yp_ingoing_hor_p2

export kerr_geo_generic_sample, kerr_geo_eccentric_sample, kerr_geo_inclined_sample
export kerr_geo_inclined_sample_dense, subsample_inclined_sample, refine_inclined_swsh_sample
export kerr_geo_eccentric_sample_dense, subsample_eccentric_sample, refine_eccentric_y_sample
export swsh_sample, y_sample_m2, y_sample_p2
export y_sample_m2_isem, y_sample_p2_isem
export integrand_generic_sample_m2, integrand_generic_sample_p2, integrand_eccentric_sample_m2, integrand_eccentric_sample_p2, integrand_inclined_sample_m2, integrand_inclined_sample_p2
export kerr_geo_generic_sample_cheby, kerr_geo_eccentric_sample_cheby, kerr_geo_inclined_sample_cheby
export integrand_generic_sample_cheby_m2, integrand_generic_sample_cheby_p2, integrand_eccentric_sample_cheby_m2, integrand_eccentric_sample_cheby_p2, integrand_inclined_sample_cheby_m2, integrand_inclined_sample_cheby_p2
export carter_ingredients_sample

function kerr_geo_generic_sample(KG::Dict, N_sample::Int64, K_sample::Int64)
    # === Check input ===
    ispow2(N_sample) || throw(ArgumentError("N_sample must be a power of 2"))
    ispow2(K_sample) || throw(ArgumentError("K_sample must be a power of 2"))
    N_half = div(N_sample, 2)
    K_half = div(K_sample, 2)

    # === Extract data from KG ===
    a = KG["a"]
    _, r, θ, _ = KG["Trajectory"]
    _, ur, uθ, _ = KG["FourVelocity"]
    Δtr, Δtθ, Δφr, Δφθ = KG["CrossFunction"]
    qt0, qr0, qθ0, qφ0 = KG["InitialPhases"]
    Frequencies = KG["Frequencies"]
    ϒr = Frequencies["ϒr"]
    ϒθ = Frequencies["ϒθ"]
    Γ = Frequencies["ϒt"]

    # === Map phase to physical quantities ===
    rq(qr)  = r((qr - qr0) / ϒr)
    θq(qθ)  = θ((qθ - qθ0) / ϒθ)
    urq(qr) = (r((qr - qr0) / ϒr)^2 + a^2 * cos(θ((qr - qr0) / ϒr))^2) * ur((qr - qr0) / ϒr)
    uθq(qθ) = (r((qθ - qθ0) / ϒθ)^2 + a^2 * cos(θ((qθ - qθ0) / ϒθ))^2) * uθ((qθ - qθ0) / ϒθ)

    # === Sampling grids ===
    qr_full = collect(range(0, π; length=N_sample))   # full radial phase
    qθ_full = collect(range(0, π; length=K_sample))   # full polar phase
    qr_half = qr_full[1:N_half]                        # half grid for symmetric functions
    qθ_half = qθ_full[1:K_half]

    # === Helper to extend half-grid using symmetry ===
    # antisym=true -> odd symmetry (negate reverse)
    # antisym=false -> even symmetry (just reverse)
    extend_half(half::Vector{T}, antisym::Bool) where T =
        antisym ? vcat(half, -reverse(half)) : vcat(half, reverse(half))

    # === No-symmetry quantities ===
    rq_vals  = [rq(q) for q in qr_full]
    rs_vals = rstar_from_r.(a, rq_vals)
    θq_vals  = [θq(q) for q in qθ_full]
    Δtr_vals = [Δtr(q) for q in qr_full]
    Δφr_vals = [Δφr(q) for q in qr_full]  # no symmetry

    # === Symmetric quantities (compute half then extend) ===
    Δtθ_half = [Δtθ(q) for q in qθ_half]   # antisymmetric
    Δφθ_half = [Δφθ(q) for q in qθ_half]   # antisymmetric

    Δtθ_vals = extend_half(Δtθ_half, true)
    Δφθ_vals = extend_half(Δφθ_half, true)

    # === urq / uθq values ===
    urq_fwd = [urq(q) for q in qr_full]           # full forward
    urq_rev = -urq_fwd                             # odd symmetry

    uθq_half = [uθq(q) for q in qθ_half]
    uθq_fwd  = extend_half(uθq_half, false)       # even symmetry
    uθq_rev  = -reverse(uθq_fwd)                  # reversed then negated

    return Dict(
        "r" => rq_vals,
        "rs" => rs_vals,
        "θ" => θq_vals,
        "Δtr" => Δtr_vals,
        "Δtθ" => Δtθ_vals,
        "Δφr" => Δφr_vals,
        "Δφθ" => Δφθ_vals,
        "ur_fwd" => urq_fwd,
        "ur_rev" => urq_rev,
        "uθ_fwd" => uθq_fwd,
        "uθ_rev" => uθq_rev,
        "N_sample" => N_sample,
        "K_sample" => K_sample,
        "E" => KG["Energy"],
        "Lz" => KG["AngularMomentum"],
        "Γ" => Γ,
        "CrossFunction" => KG["CrossFunction"],
        "initialPhases" => KG["InitialPhases"]
    )
end

function kerr_geo_generic_sample_cheby(KG::Dict, N_sample::Int64, K_sample::Int64)
    # === Input validation ===
    ispow2(N_sample) || throw(ArgumentError("N_sample must be a power of 2"))
    ispow2(K_sample) || throw(ArgumentError("K_sample must be a power of 2"))
    N_half = div(N_sample, 2)
    K_half = div(K_sample, 2)

    # === Extract data from Kerr geometry dictionary ===
    a = KG["a"]
    _, r, θ, _ = KG["Trajectory"]
    _, ur, uθ, _ = KG["FourVelocity"]
    Δtr, Δtθ, Δφr, Δφθ = KG["CrossFunction"]
    dtr, dtθ, dφr, dφθ = KG["DerivativesCrossFunction"]
    qt0, qr0, qθ0, qφ0 = KG["InitialPhases"]
    Frequencies = KG["Frequencies"]
    ϒr = Frequencies["ϒr"]
    ϒθ = Frequencies["ϒθ"]
    Γ = Frequencies["ϒt"]

    # === Map phase variables to physical quantities ===
    rq(qr)  = r((qr - qr0) / ϒr)  # Radial coordinate from radial phase
    θq(qθ)  = θ((qθ - qθ0) / ϒθ)  # Polar angle from polar phase
    # Radial four-velocity component (scaled)
    urq(qr) = (r((qr - qr0) / ϒr)^2 + a^2 * cos(θ((qr - qr0) / ϒr))^2) * ur((qr - qr0) / ϒr)
    # Polar four-velocity component (scaled)
    uθq(qθ) = (r((qθ - qθ0) / ϒθ)^2 + a^2 * cos(θ((qθ - qθ0) / ϒθ))^2) * uθ((qθ - qθ0) / ϒθ)

    # === Chebyshev node generation ===
    # Generates Chebyshev extremal points on [a, b] (clustered at boundaries)
    function chebyshev_nodes(n::Int, a::Float64, b::Float64)
        # Chebyshev extremal points on reference interval [-1, 1]: cos(kπ/(n-1)) for k=0,...,n-1
        x_ref = [cos(k*π/(n-1)) for k in 0:n-1]
        # Map to target interval [a, b]
        return (b - a) * (x_ref .+ 1) ./ 2 .+ a
    end

    # === Sampling grids (Chebyshev-distributed) ===
    # Full radial phase grid on [0, π] with Chebyshev nodes
    qr_full = chebyshev_nodes(N_sample, 0.0, 1.0pi)
    # Full polar phase grid on [0, π] with Chebyshev nodes
    qθ_full = chebyshev_nodes(K_sample, 0.0, 1.0pi)
    
    # Half-grids for symmetric function extension (first half of nodes)
    qr_half = qr_full[1:N_half]
    qθ_half = qθ_full[1:K_half]

    # === Helper function to extend half-grids using symmetry ===
    # antisym=true: odd symmetry (negate reversed half)
    # antisym=false: even symmetry (direct reversal)
    extend_half(half::Vector{T}, antisym::Bool) where T =
        antisym ? vcat(half, -reverse(half)) : vcat(half, reverse(half))

    # === Quantities without symmetry ===
    rq_vals  = [rq(q) for q in qr_full]         # Radial coordinate values
    rs_vals = rstar_from_r.(a, rq_vals)        # Tortoise coordinate (assumes rstar_from_r is defined)
    θq_vals  = [θq(q) for q in qθ_full]        # Polar angle values
    Δtr_vals = [Δtr(q) for q in qr_full]       # Radial cross function Δt
    Δφr_vals = [Δφr(q) for q in qr_full]       # Radial cross function Δφ (no symmetry)

    # === Symmetric quantities (compute half then extend) ===
    Δtθ_half = [Δtθ(q) for q in qθ_half]       # Polar cross function Δt (antisymmetric)
    Δφθ_half = [Δφθ(q) for q in qθ_half]       # Polar cross function Δφ (antisymmetric)

    Δtθ_vals = extend_half(Δtθ_half, true)     # Extend with odd symmetry
    Δφθ_vals = extend_half(Δφθ_half, true)     # Extend with odd symmetry

    dtr_vals = [dtr(q) for q in qr_full]
    dφr_vals = [dφr(q) for q in qr_full]
    dtθ_vals = [dtθ(q) for q in qθ_full]
    dφθ_vals = [dφθ(q) for q in qθ_full]

    # === urq / uθq values ===
    urq_fwd = [urq(q) for q in qr_full]        # Forward radial velocity
    urq_rev = -urq_fwd                         # Reverse radial velocity (odd symmetry)

    uθq_half = [uθq(q) for q in qθ_half]       # Half-grid polar velocity
    uθq_fwd  = extend_half(uθq_half, false)    # Forward polar velocity (even symmetry)
    uθq_rev  = -reverse(uθq_fwd)               # Reverse polar velocity

    return Dict(
        "r" => rq_vals,
        "rs" => rs_vals,
        "θ" => θq_vals,
        "Δtr" => Δtr_vals,
        "Δtθ" => Δtθ_vals,
        "Δφr" => Δφr_vals,
        "Δφθ" => Δφθ_vals,
        "dtr" => dtr_vals,
        "dtθ" => dtθ_vals,
        "dφr" => dφr_vals,
        "dφθ" => dφθ_vals,
        "ur_fwd" => urq_fwd,
        "ur_rev" => urq_rev,
        "uθ_fwd" => uθq_fwd,
        "uθ_rev" => uθq_rev,
        "N_sample" => N_sample,
        "K_sample" => K_sample,
        "E" => KG["Energy"],
        "Lz" => KG["AngularMomentum"],
        "Γ" => Γ,
        "CrossFunction" => KG["CrossFunction"],
        "initialPhases" => KG["InitialPhases"]
    )
end

function kerr_geo_eccentric_sample(KG::Dict, N_sample::Int64)
    # === Input validation ===
    ispow2(N_sample) || throw(ArgumentError("N_sample must be a power of 2"))
    N_half = div(N_sample, 2)

    # === Extract data from KG dictionary ===
    a = KG["a"]
    _, r, _, _ = KG["Trajectory"]  # Omit θ components (fixed at π/2)
    _, ur, _, _ = KG["FourVelocity"]  # Omit uθ components
    Δtr, _, Δφr, _ = KG["CrossFunction"]  # Omit θ-direction cross functions
    _, qr0, _, _ = KG["InitialPhases"]  # Omit θ initial phase
    Frequencies = KG["Frequencies"]
    ϒr = Frequencies["ϒr"]
    Γ = Frequencies["ϒt"]

    # === Map phase to physical quantities (radial only) ===
    rq(qr)  = r((qr - qr0) / ϒr)  # Radial coordinate from radial phase
    # 4-velocity radial component (simplified for θ=π/2 where cos(π/2)=0)
    urq(qr) = (rq(qr)^2 + a^2 * cos(π/2)^2) * ur((qr - qr0) / ϒr)

    # === Radial sampling grid ===
    qr_full = collect(range(0, π; length=N_sample))  # Full radial phase grid
    qr_half = qr_full[1:N_half]                       # Half grid for symmetric functions

    # === Quantities without symmetry ===
    rq_vals  = [rq(q) for q in qr_full]               # Radial coordinate values
    rs_vals = rstar_from_r.(a, rq_vals)               # Tortoise coordinate from radial coordinate
    Δtr_vals = [Δtr(q) for q in qr_full]              # Radial time difference values
    Δφr_vals = [Δφr(q) for q in qr_full]              # Radial azimuthal difference (no symmetry)

    # === Radial 4-velocity values (with symmetry) ===
    urq_fwd = [urq(q) for q in qr_full]               # Forward radial velocity
    urq_rev = -urq_fwd                                # Odd symmetry for reverse direction

    # === Return results dictionary ===
    return Dict(
        "r" => rq_vals,
        "rs" => rs_vals,
        "θ" => π/2,  # θ fixed at equatorial plane
        "Δtr" => Δtr_vals,
        "Δφr" => Δφr_vals,
        "ur_fwd" => urq_fwd,
        "ur_rev" => urq_rev,
        "N_sample" => N_sample,
        "K_sample" => 0,  # No θ-direction sampling
        "E" => KG["Energy"],
        "Lz" => KG["AngularMomentum"],
        "Γ" => Γ,
        "CrossFunction" => (Δtr, Δφr),  # Preserve function structure with zeros
        "InitialPhases" => KG["InitialPhases"] 
    )
end

function kerr_geo_eccentric_sample_dense(KG::Dict, N_sample::Int64)
    ispow2(N_sample) || throw(ArgumentError("N_sample must be a power of 2"))
    N_points = N_sample + 1

    a = KG["a"]
    _, r, _, _ = KG["Trajectory"]
    _, ur, _, _ = KG["FourVelocity"]
    Δtr, _, Δφr, _ = KG["CrossFunction"]
    _, qr0, _, _ = KG["InitialPhases"]
    Frequencies = KG["Frequencies"]
    ϒr = Frequencies["ϒr"]
    Γ = Frequencies["ϒt"]

    rq(qr) = r((qr - qr0) / ϒr)
    urq(qr) = (rq(qr)^2 + a^2 * cos(π / 2)^2) * ur((qr - qr0) / ϒr)
    qr_full = collect(range(0, π; length = N_points))
    rq_vals = [rq(q) for q in qr_full]
    rs_vals = rstar_from_r.(a, rq_vals)
    Δtr_vals = [Δtr(q) for q in qr_full]
    Δφr_vals = [Δφr(q) for q in qr_full]
    urq_fwd = [urq(q) for q in qr_full]
    urq_rev = -urq_fwd

    return Dict(
        "a" => a,
        "p" => KG["p"],
        "e" => KG["e"],
        "r" => rq_vals,
        "rs" => rs_vals,
        "θ" => π/2,
        "Δtr" => Δtr_vals,
        "Δφr" => Δφr_vals,
        "ur_fwd" => urq_fwd,
        "ur_rev" => urq_rev,
        "N_sample" => N_points,
        "K_sample" => 0,
        "E" => KG["Energy"],
        "Lz" => KG["AngularMomentum"],
        "Γ" => Γ,
        "Frequencies" => Frequencies,
        "Trajectory" => KG["Trajectory"],
        "FourVelocity" => KG["FourVelocity"],
        "CrossFunction" => (Δtr, Δφr),
        "InitialPhases" => KG["InitialPhases"],
        "initialPhases" => KG["InitialPhases"]
    )
end

function subsample_eccentric_sample(KG_sample::Dict, N_sample::Int64)
    ispow2(N_sample) || throw(ArgumentError("N_sample must be a power of 2"))
    full_points = length(KG_sample["r"])
    target_points = N_sample + 1
    full_intervals = full_points - 1
    target_intervals = target_points - 1
    full_intervals % target_intervals == 0 || throw(ArgumentError("N_sample does not divide the stored eccentric grid"))
    step = div(full_intervals, target_intervals)
    idx = 1:step:full_points
    sample = Dict{Any, Any}()
    for key in ("r", "rs", "Δtr", "Δφr", "ur_fwd", "ur_rev")
        sample[key] = KG_sample[key][idx]
    end
    sample["θ"] = KG_sample["θ"]
    sample["N_sample"] = target_points
    sample["K_sample"] = 0
    sample["a"] = KG_sample["a"]
    sample["p"] = KG_sample["p"]
    sample["e"] = KG_sample["e"]
    sample["E"] = KG_sample["E"]
    sample["Lz"] = KG_sample["Lz"]
    sample["Γ"] = KG_sample["Γ"]
    sample["Frequencies"] = KG_sample["Frequencies"]
    sample["Trajectory"] = KG_sample["Trajectory"]
    sample["FourVelocity"] = KG_sample["FourVelocity"]
    sample["CrossFunction"] = KG_sample["CrossFunction"]
    sample["InitialPhases"] = KG_sample["InitialPhases"]
    sample["initialPhases"] = KG_sample["InitialPhases"]
    return sample
end

function kerr_geo_eccentric_sample_cheby(KG::Dict, N_sample::Int64)
    # === Input validation ===
    ispow2(N_sample) || throw(ArgumentError("N_sample must be a power of 2"))
    N_half = div(N_sample, 2)

    # === Extract data from KG dictionary ===
    a = KG["a"]
    _, r, _, _ = KG["Trajectory"]  # Omit θ components (fixed at π/2)
    _, ur, _, _ = KG["FourVelocity"]  # Omit uθ components
    Δtr, _, Δφr, _ = KG["CrossFunction"]  # Omit θ-direction cross functions
    dtr, _, dφr, _ = KG["DerivativesCrossFunction"]  # Added cross function derivatives
    _, qr0, _, _ = KG["InitialPhases"]  # Omit θ initial phase
    Frequencies = KG["Frequencies"]
    ϒr = Frequencies["ϒr"]
    Γ = Frequencies["ϒt"]

    # === Map phase to physical quantities (radial only) ===
    rq(qr)  = r((qr - qr0) / ϒr)  # Radial coordinate from radial phase
    # 4-velocity radial component (simplified for θ=π/2 where cos(π/2)=0)
    urq(qr) = (rq(qr)^2 + a^2 * cos(π/2)^2) * ur((qr - qr0) / ϒr)

    # === Chebyshev node generation ===
    # Generates Chebyshev extremal points on [a, b] (clustered at boundaries)
    function chebyshev_nodes(n::Int, a::Float64, b::Float64)
        # Chebyshev extremal points on reference interval [-1, 1]: cos(kπ/(n-1)) for k=0,...,n-1
        x_ref = [cos(k*π/(n-1)) for k in 0:n-1]
        # Map to target interval [a, b]
        return (b - a) * (x_ref .+ 1) ./ 2 .+ a
    end

    # === Radial sampling grid (Chebyshev-distributed) ===
    qr_full = chebyshev_nodes(N_sample, 0.0, 1.0pi)  # Full radial phase grid with Chebyshev nodes
    qr_half = qr_full[1:N_half]                  # Half grid for symmetric functions

    # === Quantities without symmetry ===
    rq_vals  = [rq(q) for q in qr_full]          # Radial coordinate values
    rs_vals = rstar_from_r.(a, rq_vals)          # Tortoise coordinate from radial coordinate
    Δtr_vals = [Δtr(q) for q in qr_full]         # Radial time difference values
    Δφr_vals = [Δφr(q) for q in qr_full]         # Radial azimuthal difference (no symmetry)
    
    # === Cross function derivative values (added) ===
    dtr_vals = [dtr(q) for q in qr_full]         # Derivative of radial cross function Δt
    dφr_vals = [dφr(q) for q in qr_full]         # Derivative of radial cross function Δφ

    # === Radial 4-velocity values (with symmetry) ===
    urq_fwd = [urq(q) for q in qr_full]          # Forward radial velocity
    urq_rev = -urq_fwd                           # Odd symmetry for reverse direction

    # === Return results dictionary ===
    return Dict(
        "r" => rq_vals,
        "rs" => rs_vals,
        "θ" => π/2,  # θ fixed at equatorial plane
        "Δtr" => Δtr_vals,
        "Δφr" => Δφr_vals,
        "dtr" => dtr_vals,  # Added derivative output
        "dφr" => dφr_vals,  # Added derivative output
        "ur_fwd" => urq_fwd,
        "ur_rev" => urq_rev,
        "N_sample" => N_sample,
        "K_sample" => 0,  # No θ-direction sampling
        "E" => KG["Energy"],
        "Lz" => KG["AngularMomentum"],
        "Γ" => Γ,
        "CrossFunction" => (Δtr, Δφr),  # Preserve original function structure
        "DerivativesCrossFunction" => (dtr, dφr),  # Added derivatives to output dict
        "InitialPhases" => KG["InitialPhases"]
    )
end

function kerr_geo_inclined_sample(KG::Dict, K_sample::Int64)
    # === Input validation ===
    ispow2(K_sample) || throw(ArgumentError("K_sample must be a power of 2"))
    K_half = div(K_sample, 2)

    # === Extract data from KG dictionary ===
    a = KG["a"]
    _, _, θ, _ = KG["Trajectory"]  # Omit r components (fixed at p)
    _, _, uθ, _ = KG["FourVelocity"]  # Omit ur components
    _, Δtθ, _, Δφθ = KG["CrossFunction"]  # Omit r-direction cross functions
    _, _, qθ0, _ = KG["InitialPhases"]  # Omit r initial phase
    Frequencies = KG["Frequencies"]
    ϒθ = Frequencies["ϒθ"]
    Γ = Frequencies["ϒt"]
    p = KG["p"]

    # === Map phase to physical quantities (polar only) ===
    θq(qθ)  = θ((qθ - qθ0) / ϒθ)  # Polar angle from polar phase
    # 4-velocity polar component (simplified for fixed r=p)
    uθq(qθ) = (p^2 + a^2 * cos(θq(qθ))^2) * uθ((qθ - qθ0) / ϒθ)

    # === Polar sampling grid ===
    qθ_full = collect(range(0, π; length=K_sample))  # Full polar phase grid
    qθ_half = qθ_full[1:K_half]                       # Half grid for symmetric functions

    # === Helper to extend half-grid using symmetry ===
    # antisym=true -> odd symmetry (negate reverse)
    # antisym=false -> even symmetry (just reverse)
    extend_half(half::Vector{T}, antisym::Bool) where T =
        antisym ? vcat(half, -reverse(half)) : vcat(half, reverse(half))

    # === No-symmetry quantities ===
    θq_vals  = [θq(q) for q in qθ_full]               # Polar angle values
    rs = rstar_from_r(a, p)
    Δtθ_half = [Δtθ(q) for q in qθ_half]   # antisymmetric
    Δφθ_half = [Δφθ(q) for q in qθ_half]   # antisymmetric

    Δtθ_vals = extend_half(Δtθ_half, true)
    Δφθ_vals = extend_half(Δφθ_half, true)

    # === Polar 4-velocity values (with symmetry) ===
    uθq_half = [uθq(q) for q in qθ_half]
    uθq_fwd  = extend_half(uθq_half, false)           # Even symmetry for forward direction
    uθq_rev  = - reverse(uθq_fwd)                      # Reversed then negated for reverse direction

    # === Return results dictionary ===
    return Dict(
        "r" => p,  # Radial coordinate fixed at p
        "rs" => rs,
        "θ" => θq_vals,
        "Δtθ" => Δtθ_vals,
        "Δφθ" => Δφθ_vals,
        "uθ_fwd" => uθq_fwd,
        "uθ_rev" => uθq_rev,
        "K_sample" => K_sample,
        "E" => KG["Energy"],
        "Lz" => KG["AngularMomentum"],
        "Γ" => Γ,
        "CrossFunction" => (Δtθ, Δφθ),
        "InitialPhases" => KG["InitialPhases"]
    )
end

function kerr_geo_inclined_sample_dense(KG::Dict, x::Real, K_interval::Int64)
    ispow2(K_interval) || throw(ArgumentError("K_interval must be a power of 2"))
    K_points = K_interval + 1
    a = KG["a"]
    p = KG["p"]
    _, _, θ, _ = KG["Trajectory"]
    _, _, uθ, _ = KG["FourVelocity"]
    _, Δtθ, _, Δφθ = KG["CrossFunction"]
    _, _, qθ0, _ = KG["InitialPhases"]
    Frequencies = KG["Frequencies"]
    ϒθ = Frequencies["ϒθ"]
    Γ = Frequencies["ϒt"]
    θ_fun(qθ) = θ((qθ - qθ0) / ϒθ)
    uθ_fun(qθ) = (p^2 + a^2 * cos(θ_fun(qθ))^2) * uθ((qθ - qθ0) / ϒθ)
    qθ_full = collect(range(0, π; length = K_points))
    θ_vals = [θ_fun(q) for q in qθ_full]
    Δtθ_vals = [Δtθ(q) for q in qθ_full]
    Δφθ_vals = [Δφθ(q) for q in qθ_full]
    uθ_fwd = [uθ_fun(q) for q in qθ_full]
    return Dict(
        "a" => a,
        "p" => p,
        "x" => x,
        "r" => p,
        "rs" => rstar_from_r(a, p),
        "θ" => θ_vals,
        "Δtθ" => Δtθ_vals,
        "Δφθ" => Δφθ_vals,
        "uθ_fwd" => uθ_fwd,
        "uθ_rev" => -uθ_fwd,
        "K_interval" => K_interval,
        "K_sample" => K_points,
        "E" => KG["Energy"],
        "Lz" => KG["AngularMomentum"],
        "Γ" => Γ,
        "Frequencies" => Frequencies,
        "Trajectory" => KG["Trajectory"],
        "FourVelocity" => KG["FourVelocity"],
        "CrossFunction" => KG["CrossFunction"],
        "InitialPhases" => KG["InitialPhases"],
    )
end

function subsample_inclined_sample(KG_sample::Dict, K_interval::Int64)
    ispow2(K_interval) || throw(ArgumentError("K_interval must be a power of 2"))
    full_points = length(KG_sample["θ"])
    target_points = K_interval + 1
    full_intervals = full_points - 1
    target_intervals = target_points - 1
    full_intervals % target_intervals == 0 || throw(ArgumentError("K_interval does not divide the stored inclined grid"))
    step = div(full_intervals, target_intervals)
    idx = 1:step:full_points
    sample = Dict{Any, Any}()
    for key in ("θ", "Δtθ", "Δφθ", "uθ_fwd", "uθ_rev")
        sample[key] = view(KG_sample[key], idx)
    end
    for key in ("a", "p", "x", "r", "rs", "E", "Lz", "Γ", "Frequencies", "Trajectory", "FourVelocity", "CrossFunction", "InitialPhases")
        sample[key] = KG_sample[key]
    end
    sample["K_interval"] = K_interval
    sample["K_sample"] = target_points
    return sample
end

@inline function _swsh_values(SH, θ, s, m, c, λ)
    s0 = SH(θ, 0.0)
    s1 = SH(θ, 0.0; theta_derivative = 1)
    sθ = sin(θ)
    cθ = cos(θ)
    cotθ = cθ / sθ
    s2 = -cotθ * s1 + (((m + s * cθ)^2) / (sθ^2) + 2 * c * s * cθ + c^2 * (sθ^2) - s - 2 * c * m - λ) * s0
    return s0, s1, s2
end

function refine_inclined_swsh_sample(SH, prev_SHsamp::Dict, next_sample::Dict)
    prev_points = length(prev_SHsamp["S0"])
    next_points = length(next_sample["θ"])
    next_points >= prev_points || throw(ArgumentError("next_sample must not be coarser than prev_SHsamp"))
    next_points == prev_points && return prev_SHsamp
    refined = Dict{String, Vector{ComplexF64}}(
        "S0" => Vector{ComplexF64}(undef, next_points),
        "S1" => Vector{ComplexF64}(undef, next_points),
        "S2" => Vector{ComplexF64}(undef, next_points),
    )
    @inbounds for i in 1:next_points
        if isodd(i)
            j = (i + 1) >>> 1
            refined["S0"][i] = prev_SHsamp["S0"][j]
            refined["S1"][i] = prev_SHsamp["S1"][j]
            refined["S2"][i] = prev_SHsamp["S2"][j]
        else
            s0, s1, s2 = _swsh_values(SH, next_sample["θ"][i], SH.params.s, SH.params.m, SH.params.c, SH.lambda)
            refined["S0"][i] = s0
            refined["S1"][i] = s1
            refined["S2"][i] = s2
        end
    end
    return refined
end

function kerr_geo_inclined_sample_cheby(KG::Dict, K_sample::Int64)
    # === Input validation ===
    ispow2(K_sample) || throw(ArgumentError("K_sample must be a power of 2"))
    K_half = div(K_sample, 2)

    # === Extract data from KG dictionary ===
    a = KG["a"]
    _, _, θ, _ = KG["Trajectory"]  # Omit r components (fixed at p)
    _, _, uθ, _ = KG["FourVelocity"]  # Omit ur components
    _, Δtθ, _, Δφθ = KG["CrossFunction"]  # Omit r-direction cross functions
    _, dtθ, _, dφθ = KG["DerivativesCrossFunction"]  # Added cross function derivatives
    _, _, qθ0, _ = KG["InitialPhases"]  # Omit r initial phase
    Frequencies = KG["Frequencies"]
    ϒθ = Frequencies["ϒθ"]
    Γ = Frequencies["ϒt"]
    p = KG["p"]

    # === Map phase to physical quantities (polar only) ===
    θq(qθ)  = θ((qθ - qθ0) / ϒθ)  # Polar angle from polar phase
    # 4-velocity polar component (simplified for fixed r=p)
    uθq(qθ) = (p^2 + a^2 * cos(θq(qθ))^2) * uθ((qθ - qθ0) / ϒθ)

    # === Chebyshev node generation ===
    # Generates Chebyshev extremal points on [a, b] (clustered at boundaries)
    function chebyshev_nodes(n::Int, a::Float64, b::Float64)
        # Chebyshev extremal points on reference interval [-1, 1]: cos(kπ/(n-1)) for k=0,...,n-1
        x_ref = [cos(k*π/(n-1)) for k in 0:n-1]
        # Map to target interval [a, b]
        return (b - a) * (x_ref .+ 1) ./ 2 .+ a
    end

    # === Polar sampling grid (Chebyshev-distributed) ===
    qθ_full = chebyshev_nodes(K_sample, 0.0, 1.0pi)  # Full polar phase grid with Chebyshev nodes
    qθ_half = qθ_full[1:K_half]                  # Half grid for symmetric functions

    # === Helper to extend half-grid using symmetry ===
    # antisym=true -> odd symmetry (negate reverse)
    # antisym=false -> even symmetry (just reverse)
    extend_half(half::Vector{T}, antisym::Bool) where T =
        antisym ? vcat(half, -reverse(half)) : vcat(half, reverse(half))

    # === No-symmetry quantities ===
    θq_vals  = [θq(q) for q in qθ_full]               # Polar angle values
    rs_vals = fill(rstar_from_r(a, p), K_sample)      # Tortoise coordinate (fixed r=p)
    Δtθ_vals = [Δtθ(q) for q in qθ_full]              # Polar time difference values
    Δφθ_vals = [Δφθ(q) for q in qθ_full]              # Polar azimuthal difference (no symmetry)
    
    # === Cross function derivative values (added) ===
    dtθ_vals = [dtθ(q) for q in qθ_full]              # Derivative of polar cross function Δt
    dφθ_vals = [dφθ(q) for q in qθ_full]              # Derivative of polar cross function Δφ

    # === Polar 4-velocity values (with symmetry) ===
    uθq_half = [uθq(q) for q in qθ_half]
    uθq_fwd  = extend_half(uθq_half, false)           # Even symmetry for forward direction
    uθq_rev  = -reverse(uθq_fwd)                      # Reversed then negated for reverse direction

    # === Return results dictionary ===
    return Dict(
        "r" => p,  # Radial coordinate fixed at p
        "rs" => rs_vals,
        "θ" => θq_vals,
        "Δtθ" => Δtθ_vals,
        "Δφθ" => Δφθ_vals,
        "dtθ" => dtθ_vals,  # Added derivative output
        "dφθ" => dφθ_vals,  # Added derivative output
        "uθ_fwd" => uθq_fwd,
        "uθ_rev" => uθq_rev,
        "K_sample" => K_sample,
        "N_sample" => 0,  # No r-direction sampling
        "E" => KG["Energy"],
        "Lz" => KG["AngularMomentum"],
        "Γ" => Γ,
        "CrossFunction" => (Δtθ, Δφθ),  # Preserve original function structure
        "DerivativesCrossFunction" => (dtθ, dφθ),  # Added derivatives to output dict
        "InitialPhases" => KG["InitialPhases"]
    )
end

function swsh_sample(SH, KG_sample::Dict)
    θ_vals = KG_sample["θ"]
    s = SH.params.s
    m = SH.params.m
    c = SH.params.c
    λ = SH.lambda
    n = length(θ_vals)
    S0_vals = Vector{ComplexF64}(undef, n)
    S1_vals = Vector{ComplexF64}(undef, n)
    S2_vals = Vector{ComplexF64}(undef, n)
    @inbounds for i in eachindex(θ_vals)
        s0, s1, s2 = _swsh_values(SH, θ_vals[i], s, m, c, λ)
        S0_vals[i] = s0
        S1_vals[i] = s1
        S2_vals[i] = s2
    end
    return Dict(
        "S0" => S0_vals,
        "S1" => S1_vals,
        "S2" => S2_vals
    )
end

function y_sample_m2(Y_soln, KG_sample::Dict)
    a = Y_soln.mode.a
    rs_vals = KG_sample["rs"]
    sol = Y_soln.solution
    Y = sol.Y_hor
    Yp = sol.Yp_hor
    X = sol.X
    Y_vals = [Y(rs) for rs in rs_vals]
    Yp_vals = [Yp(rs) for rs in rs_vals]
    X_vals = [X(rs) for rs in rs_vals]
    return Dict(
        "Y" => Y_vals,
        "Yp" => Yp_vals,
        "X" => X_vals,
        "params" => Y_soln.mode,
        "Binc" => Y_soln.asymptotic.Binc
    )
end

function y_sample_p2(Y_soln, KG_sample::Dict)
    a = Y_soln.mode.a
    rs_vals = KG_sample["rs"]
    sol = Y_soln.solution
    Y = sol.Y_inf
    Yp = sol.Yp_inf
    X = sol.X
    Y_vals = [Y(rs) for rs in rs_vals]
    Yp_vals = [Yp(rs) for rs in rs_vals]
    X_vals = [X(rs) for rs in rs_vals]
    return Dict(
        "Y" => Y_vals,
        "Yp" => Yp_vals,
        "X" => X_vals,
        "params" => Y_soln.mode,
        "Cinc" => Y_soln.asymptotic.Cinc
    )
end

@inline function _isem_gsn_incidence_amplitude(Y_soln)
    s = Y_soln.mode.s
    if s == -2
        return Y_soln.incidence_amplitude / ConversionFactors.Binc(s, Y_soln.mode.m, Y_soln.mode.a, Y_soln.mode.omega, Y_soln.mode.lambda)
    elseif s == 2
        return Y_soln.incidence_amplitude / ConversionFactors.Cinc(s, Y_soln.mode.m, Y_soln.mode.a, Y_soln.mode.omega, Y_soln.mode.lambda)
    else
        error("ISEM Y interfaces currently support only s = -2 and s = 2.")
    end
end

function _y_sample_isem(Y_soln, KG_sample::Dict, amplitude_key::String; threaded::Bool = false)
    r_vals = KG_sample["r"]
    n = length(r_vals)
    Y_vals = Vector{ComplexF64}(undef, n)
    Yp_vals = Vector{ComplexF64}(undef, n)
    X_vals = Vector{ComplexF64}(undef, n)
    Y = Y_soln.Y_solution
    if threaded
        Threads.@threads :static for i in eachindex(r_vals)
            Yi, Ypi, Xi, _ = Y(r_vals[i])
            Y_vals[i] = Yi
            Yp_vals[i] = Ypi
            X_vals[i] = Xi
        end
    else
        @inbounds for i in eachindex(r_vals)
            Yi, Ypi, Xi, _ = Y(r_vals[i])
            Y_vals[i] = Yi
            Yp_vals[i] = Ypi
            X_vals[i] = Xi
        end
    end
    return Dict(
        "Y" => Y_vals,
        "Yp" => Yp_vals,
        "X" => X_vals,
        "params" => Y_soln.mode,
        amplitude_key => _isem_gsn_incidence_amplitude(Y_soln)
    )
end

function y_sample_m2_isem(Y_soln, KG_sample::Dict)
    return _y_sample_isem(Y_soln, KG_sample, "Binc")
end

function y_sample_m2_isem(Y_soln, KG_sample::Dict, N_sample::Int64)
    return y_sample_m2_isem(Y_soln, subsample_eccentric_sample(KG_sample, N_sample))
end

function y_sample_p2_isem(Y_soln, KG_sample::Dict)
    return _y_sample_isem(Y_soln, KG_sample, "Cinc")
end

function y_sample_p2_isem(Y_soln, KG_sample::Dict, N_sample::Int64)
    return y_sample_p2_isem(Y_soln, subsample_eccentric_sample(KG_sample, N_sample))
end

function y_sample_m2_isem_threaded(Y_soln, KG_sample::Dict)
    return _y_sample_isem(Y_soln, KG_sample, "Binc"; threaded = true)
end

function y_sample_p2_isem_threaded(Y_soln, KG_sample::Dict)
    return _y_sample_isem(Y_soln, KG_sample, "Cinc"; threaded = true)
end

function swsh_sample_threaded(SH, KG_sample::Dict)
    θ_vals = KG_sample["θ"]
    s = SH.params.s
    m = SH.params.m
    c = SH.params.c
    λ = SH.lambda
    n = length(θ_vals)
    S0_vals = Vector{ComplexF64}(undef, n)
    S1_vals = Vector{ComplexF64}(undef, n)
    S2_vals = Vector{ComplexF64}(undef, n)
    Threads.@threads :static for i in eachindex(θ_vals)
        s0, s1, s2 = _swsh_values(SH, θ_vals[i], s, m, c, λ)
        S0_vals[i] = s0
        S1_vals[i] = s1
        S2_vals[i] = s2
    end
    return Dict(
        "S0" => S0_vals,
        "S1" => S1_vals,
        "S2" => S2_vals
    )
end

function _copy_y_metadata!(out::Dict, Y_sample::Dict)
    out["params"] = Y_sample["params"]
    if haskey(Y_sample, "Binc")
        out["Binc"] = Y_sample["Binc"]
    end
    if haskey(Y_sample, "Cinc")
        out["Cinc"] = Y_sample["Cinc"]
    end
    return out
end

function _refine_y_sample(Y_soln, Y_sample::Dict, next_sample::Dict; threaded::Bool = false)
    prev_points = length(Y_sample["Y"])
    next_points = length(next_sample["r"])
    next_points >= prev_points || throw(ArgumentError("next_sample must not be coarser than prev Y sample"))
    next_points == prev_points && return Y_sample
    out = Dict{Any, Any}()
    out["Y"] = Vector{ComplexF64}(undef, next_points)
    out["Yp"] = Vector{ComplexF64}(undef, next_points)
    out["X"] = Vector{ComplexF64}(undef, next_points)
    Y = Y_soln.Y_solution
    if threaded
        Threads.@threads :static for i in 1:next_points
            if isodd(i)
                j = (i + 1) >>> 1
                out["Y"][i] = Y_sample["Y"][j]
                out["Yp"][i] = Y_sample["Yp"][j]
                out["X"][i] = Y_sample["X"][j]
            else
                Yi, Ypi, Xi, _ = Y(next_sample["r"][i])
                out["Y"][i] = Yi
                out["Yp"][i] = Ypi
                out["X"][i] = Xi
            end
        end
    else
        @inbounds for i in 1:next_points
            if isodd(i)
                j = (i + 1) >>> 1
                out["Y"][i] = Y_sample["Y"][j]
                out["Yp"][i] = Y_sample["Yp"][j]
                out["X"][i] = Y_sample["X"][j]
            else
                Yi, Ypi, Xi, _ = Y(next_sample["r"][i])
                out["Y"][i] = Yi
                out["Yp"][i] = Ypi
                out["X"][i] = Xi
            end
        end
    end
    _copy_y_metadata!(out, Y_sample)
    return out
end

function refine_eccentric_y_sample(Y_soln, Y_sample::Dict, KG_sample::Dict, N_sample::Int64)
    return _refine_y_sample(Y_soln, Y_sample, subsample_eccentric_sample(KG_sample, 2N_sample))
end

@inline function _generic_integrand_sample_m2_fast(KG_samp, Y_samp, SH_samp, n::Int64, k::Int64)
    r_vec = KG_samp["r"]::Vector{Float64}
    rs_vec = KG_samp["rs"]::Vector{Float64}
    θ_vec = KG_samp["θ"]::Vector{Float64}
    N = KG_samp["N_sample"]::Int
    K = KG_samp["K_sample"]::Int
    Γ = KG_samp["Γ"]::Float64
    qt0, qr0, qθ0, qφ0 = KG_samp["initialPhases"]::NTuple{4, Float64}
    Δtr_func, Δtθ_func, Δφr_func, Δφθ_func = KG_samp["CrossFunction"]::Vector{Function}
    urp = KG_samp["ur_fwd"]::Vector{Float64}
    urm = KG_samp["ur_rev"]::Vector{Float64}
    Δtr = KG_samp["Δtr"]::Vector{Float64}
    Δtθ = KG_samp["Δtθ"]::Vector{Float64}
    Δφr = KG_samp["Δφr"]::Vector{Float64}
    Δφθ = KG_samp["Δφθ"]::Vector{Float64}
    uθp = KG_samp["uθ_fwd"]::Vector{Float64}
    uθm = KG_samp["uθ_rev"]::Vector{Float64}
    E = KG_samp["E"]::Float64
    Lz = KG_samp["Lz"]::Float64
    S0 = SH_samp["S0"]::Vector{ComplexF64}
    S1 = SH_samp["S1"]::Vector{ComplexF64}
    S2 = SH_samp["S2"]::Vector{ComplexF64}
    Y = Y_samp["Y"]::Vector{ComplexF64}
    Yp = Y_samp["Yp"]::Vector{ComplexF64}
    X = Y_samp["X"]::Vector{ComplexF64}
    Binc = Y_samp["Binc"]::ComplexF64
    params = Y_samp["params"]::Mode
    m = params.m::Int
    a = params.a::Float64
    ω = params.omega::Float64

    rp = 1.0 + sqrt(1.0 - a^2)
    rm = 1.0 - sqrt(1.0 - a^2)
    inv4π2 = 1.0 / (4.0 * π^2)
    inv2κ = 1.0 / (2.0 * sqrt(1.0 - a^2))
    qr_step = N > 1 ? π / (N - 1) : 0.0
    qθ_step = K > 1 ? π / (K - 1) : 0.0

    sinθ = Vector{Float64}(undef, K)
    cosθ = Vector{Float64}(undef, K)
    cotθ = Vector{Float64}(undef, K)
    inv_sinθ = Vector{Float64}(undef, K)
    inv_sin2θ = Vector{Float64}(undef, K)
    L1 = Vector{Float64}(undef, K)
    L2 = Vector{Float64}(undef, K)
    L2p = Vector{Float64}(undef, K)
    L2S = Vector{ComplexF64}(undef, K)
    L1Sp = Vector{ComplexF64}(undef, K)
    L1L2S = Vector{ComplexF64}(undef, K)
    theta_const = Vector{ComplexF64}(undef, K)
    theta_L1 = Vector{ComplexF64}(undef, K)
    θphase = Vector{Float64}(undef, K)
    termM = Vector{ComplexF64}(undef, K)
    S0v = Vector{ComplexF64}(undef, K)
    S1v = Vector{ComplexF64}(undef, K)
    S2v = Vector{ComplexF64}(undef, K)

    @inbounds for j in 1:K
        θ = θ_vec[j]
        sθ = sin(θ)
        cθ = cos(θ)
        cotθj = cθ / sθ
        invsθ = 1.0 / sθ
        invs2θ = invsθ * invsθ
        L1j = (-m) * invsθ + a * ω * sθ + cotθj
        L2j = (-m) * invsθ + a * ω * sθ + 2.0 * cotθj
        L2pj = m * cθ * invs2θ + a * ω * cθ - 2.0 * invs2θ
        S0j = S0[j]
        S1j = S1[j]
        S2j = S2[j]
        L2Sj = S1j + L2j * S0j
        L1Spj = S2j + L1j * S1j
        L1L2Sj = L1Spj + L2pj * S0j + L2j * S1j + L1j * L2j * S0j
        sinθ[j] = sθ
        cosθ[j] = cθ
        cotθ[j] = cotθj
        inv_sinθ[j] = invsθ
        inv_sin2θ[j] = invs2θ
        L1[j] = L1j
        L2[j] = L2j
        L2p[j] = L2pj
        S0v[j] = S0j
        S1v[j] = S1j
        S2v[j] = S2j
        L2S[j] = L2Sj
        L1Sp[j] = L1Spj
        L1L2S[j] = L1L2Sj
        theta_const[j] = 3im * a * cθ * S0j + 2im * a * sθ * S1j - 1im * a * sθ * L2j * S0j
        theta_L1[j] = 3im * a * sθ * S0j
        θphase[j] = ω * Δtθ[j] - m * Δφθ[j] + k * (K > 1 ? (j - 1) * qθ_step : 0.0)
        termM[j] = im * sθ * (a * E - Lz * invs2θ)
    end

    rphase = Vector{Float64}(undef, N)
    phase = Vector{ComplexF64}(undef, N)
    Np = Vector{Float64}(undef, N)
    Nm = Vector{Float64}(undef, N)
    Yv = Vector{ComplexF64}(undef, N)
    Ypv = Vector{ComplexF64}(undef, N)
    Xv = Vector{ComplexF64}(undef, N)

    @inbounds for i in 1:N
        r = r_vec[i]
        r2 = r * r
        Δ = r2 - 2.0 * r + a^2
        qr_i = N > 1 ? (i - 1) * qr_step : 0.0
        rphase[i] = ω * Δtr[i] - m * Δφr[i] + n * qr_i
        phase[i] = exp(im * ω * rs_vec[i] - im * a * m * log((r - rp) / (r - rm)) * inv2κ)
        numer = E * (r2 + a^2) - a * Lz
        Np[i] = (numer + urp[i]) / Δ
        Nm[i] = (numer + urm[i]) / Δ
        Yv[i] = Y[i]
        Ypv[i] = Yp[i]
        Xv[i] = X[i]
    end

    Xi = Float64(m * (Δφr_func(qr0) + Δφθ_func(qθ0) - qφ0) - ω * (Δtr_func(qr0) + Δtθ_func(qθ0) - qt0) - k * qθ0 - n * qr0)
    prefactor = 4im * π * ω * exp(im * Xi) / (Γ * Binc)

    out = Matrix{ComplexF64}(undef, N, K)
    @inbounds for j in 1:K
        sθ = sinθ[j]
        cθ = cosθ[j]
        invsθ = inv_sinθ[j]
        invs2θ = inv_sin2θ[j]
        L1j = L1[j]
        L2j = L2[j]
        L2pj = L2p[j]
        S0j = S0v[j]
        S1j = S1v[j]
        S2j = S2v[j]
        L2Sj = L2S[j]
        L1Spj = L1Sp[j]
        L1L2Sj = L1L2S[j]
        theta_const_j = theta_const[j]
        theta_L1_j = theta_L1[j]
        phaseθ = θphase[j]
        Mpj = termM[j] + uθp[j]
        Mmj = termM[j] + uθm[j]
        expθp = exp(im * phaseθ)
        expθm = exp(-im * phaseθ)

        for i in 1:N
            r = r_vec[i]
            r2 = r * r
            ρ = -1.0 / (r - im * a * cθ)
            ρbar = -1.0 / (r + im * a * cθ)
            invρ = 1.0 / ρ
            Yi = Yv[i]
            Ypi = Ypv[i]
            Xi = Xv[i]
            phase_r = phase[i]
            rphase_i = rphase[i]
            Npi = Np[i]
            Nmi = Nm[i]

            L1pL2pS = invρ * L1L2Sj + theta_const_j + theta_L1_j * L1j

            Wnn = ρbar * L1pL2pS * r2 * Yi * phase_r / 2.0
            termA = L2Sj + im * a * sθ * (ρ - ρbar) * S0j
            termB = 2.0 * Yi + r * Ypi
            termC = L2Sj * (ρ + ρbar) * r * Yi
            Wnmbar = (-r * ρbar * (termA * termB + termC) / ρ) * phase_r
            termD = Xi / (2.0 * sqrt(r2 + a^2))
            termE = (Yi + 2.0 * r * Ypi) * phase_r
            termF = ρ * r * (2.0 * Yi + r * Ypi) * phase_r
            Wmbarmbar = S0j * ρbar * (termD + termE + termF) / ρ

            Np2 = Npi * Npi
            Nm2 = Nmi * Nmi
            Mpp2 = Mpj * Mpj
            Mm2 = Mmj * Mmj
            exp_rp = exp(im * rphase_i)
            exp_rm = exp(-im * rphase_i)
            exp_pp = exp_rp * expθp
            exp_pm = exp_rp * expθm
            exp_mp = exp_rm * expθp
            exp_mm = exp_rm * expθm
            Jpp = (Wnn * Np2 + Wnmbar * Npi * Mpj + Wmbarmbar * Mpp2) * exp_pp * inv4π2
            Jpm = (Wnn * Np2 + Wnmbar * Npi * Mmj + Wmbarmbar * Mm2) * exp_pm * inv4π2
            Jmp = (Wnn * Nm2 + Wnmbar * Nmi * Mpj + Wmbarmbar * Mpp2) * exp_mp * inv4π2
            Jmm = (Wnn * Nm2 + Wnmbar * Nmi * Mmj + Wmbarmbar * Mm2) * exp_mm * inv4π2

            out[i, j] = prefactor * (Jpp + Jpm + Jmp + Jmm)
        end
    end

    return out
end

@inline function _generic_integrand_sample_p2_fast(KG_samp, Y_samp, SH_samp, n::Int64, k::Int64)
    r_vec = KG_samp["r"]::Vector{Float64}
    rs_vec = KG_samp["rs"]::Vector{Float64}
    θ_vec = KG_samp["θ"]::Vector{Float64}
    N = KG_samp["N_sample"]::Int
    K = KG_samp["K_sample"]::Int
    Γ = KG_samp["Γ"]::Float64
    qt0, qr0, qθ0, qφ0 = KG_samp["initialPhases"]::NTuple{4, Float64}
    Δtr_func, Δtθ_func, Δφr_func, Δφθ_func = KG_samp["CrossFunction"]::Vector{Function}
    urp = KG_samp["ur_fwd"]::Vector{Float64}
    urm = KG_samp["ur_rev"]::Vector{Float64}
    Δtr = KG_samp["Δtr"]::Vector{Float64}
    Δtθ = KG_samp["Δtθ"]::Vector{Float64}
    Δφr = KG_samp["Δφr"]::Vector{Float64}
    Δφθ = KG_samp["Δφθ"]::Vector{Float64}
    uθp = KG_samp["uθ_fwd"]::Vector{Float64}
    uθm = KG_samp["uθ_rev"]::Vector{Float64}
    E = KG_samp["E"]::Float64
    Lz = KG_samp["Lz"]::Float64
    S0 = SH_samp["S0"]::Vector{ComplexF64}
    S1 = SH_samp["S1"]::Vector{ComplexF64}
    S2 = SH_samp["S2"]::Vector{ComplexF64}
    Y = Y_samp["Y"]::Vector{ComplexF64}
    Yp = Y_samp["Yp"]::Vector{ComplexF64}
    X = Y_samp["X"]::Vector{ComplexF64}
    Cinc = Y_samp["Cinc"]::ComplexF64
    params = Y_samp["params"]::Mode
    m = params.m::Int
    a = params.a::Float64
    ω = params.omega::Float64
    λ = params.lambda::Float64

    rp = 1.0 + sqrt(1.0 - a^2)
    rm = 1.0 - sqrt(1.0 - a^2)
    inv4π2 = 1.0 / (4.0 * π^2)
    inv2κ = 1.0 / (2.0 * sqrt(1.0 - a^2))
    qr_step = N > 1 ? π / (N - 1) : 0.0
    qθ_step = K > 1 ? π / (K - 1) : 0.0

    sinθ = Vector{Float64}(undef, K)
    cosθ = Vector{Float64}(undef, K)
    cotθ = Vector{Float64}(undef, K)
    inv_sinθ = Vector{Float64}(undef, K)
    inv_sin2θ = Vector{Float64}(undef, K)
    L1 = Vector{Float64}(undef, K)
    L2 = Vector{Float64}(undef, K)
    L2p = Vector{Float64}(undef, K)
    L2S = Vector{ComplexF64}(undef, K)
    L1Sp = Vector{ComplexF64}(undef, K)
    L1L2S = Vector{ComplexF64}(undef, K)
    theta_const = Vector{ComplexF64}(undef, K)
    theta_L1 = Vector{ComplexF64}(undef, K)
    θphase = Vector{Float64}(undef, K)
    termM = Vector{ComplexF64}(undef, K)
    S0v = Vector{ComplexF64}(undef, K)
    S1v = Vector{ComplexF64}(undef, K)
    S2v = Vector{ComplexF64}(undef, K)

    @inbounds for j in 1:K
        θ = θ_vec[j]
        sθ = sin(θ)
        cθ = cos(θ)
        cotθj = cθ / sθ
        invsθ = 1.0 / sθ
        invs2θ = invsθ * invsθ
        L1j = m * invsθ - a * ω * sθ + cotθj
        L2j = m * invsθ - a * ω * sθ + 2.0 * cotθj
        L2pj = -m * cθ * invs2θ - a * ω * cθ - 2.0 * invs2θ
        S0j = S0[j]
        S1j = S1[j]
        S2j = S2[j]
        L2Sj = S1j + L2j * S0j
        L1Spj = S2j + L1j * S1j
        L1L2Sj = L1Spj + L2pj * S0j + L2j * S1j + L1j * L2j * S0j
        sinθ[j] = sθ
        cosθ[j] = cθ
        cotθ[j] = cotθj
        inv_sinθ[j] = invsθ
        inv_sin2θ[j] = invs2θ
        L1[j] = L1j
        L2[j] = L2j
        L2p[j] = L2pj
        S0v[j] = S0j
        S1v[j] = S1j
        S2v[j] = S2j
        L2S[j] = L2Sj
        L1Sp[j] = L1Spj
        L1L2S[j] = L1L2Sj
        theta_const[j] = 3im * a * cθ * S0j + 2im * a * sθ * S1j - 1im * a * sθ * L2j * S0j
        theta_L1[j] = 3im * a * sθ * S0j
        θphase[j] = ω * Δtθ[j] - m * Δφθ[j] + k * (K > 1 ? (j - 1) * qθ_step : 0.0)
        termM[j] = -im * sθ * (a * E - Lz * invs2θ)
    end

    rphase = Vector{Float64}(undef, N)
    phase = Vector{ComplexF64}(undef, N)
    Lp = Vector{Float64}(undef, N)
    Lm = Vector{Float64}(undef, N)
    Yv = Vector{ComplexF64}(undef, N)
    Ypv = Vector{ComplexF64}(undef, N)
    Xv = Vector{ComplexF64}(undef, N)

    @inbounds for i in 1:N
        r = r_vec[i]
        r2 = r * r
        Δ = r2 - 2.0 * r + a^2
        qr_i = N > 1 ? (i - 1) * qr_step : 0.0
        rphase[i] = ω * Δtr[i] - m * Δφr[i] + n * qr_i
        phase[i] = exp(-im * ω * rs_vec[i] + im * a * m * log((r - rp) / (r - rm)) * inv2κ)
        numer = E * (r2 + a^2) - a * Lz
        Lp[i] = (numer - urp[i]) / Δ
        Lm[i] = (numer - urm[i]) / Δ
        Yv[i] = Y[i]
        Ypv[i] = Yp[i]
        Xv[i] = X[i]
    end

    κ = ω - a * m / (2 * rp)
    c0 = 24 + 12im * ω + λ * (10 + λ) - 12 * a * ω * (a * ω - m)
    c1 = -32im * a * m - 8im * a * m * λ + 8im * a^2 * ω * (1 + λ)
    c2 = 12 * a^2 - 24im * a * m - 24 * a^2 * m^2 + 24im * a^2 * ω + 48 * a^3 * m * ω - 24 * a^4 * ω^2
    c3 = -24im * a^3 * (a * ω - m) - 24 * a^2
    c4 = 12 * a^4
    η = c0 + c1 / rp + c2 / rp^2 + c3 / rp^3 + c4 / rp^4
    factor = 2 * sqrt(2) * rp^(3 / 2) * (((4 * ω) * (im - 4 * ω) - a * m * (im - 8 * ω) - a^2 * (m^2 + 2im * ω - 4 * ω^2)) * rp^2 +
        a^2 * (im - 4 * ω) * (a * m - 2 * ω) * rp) / (2 * rp^3 * (24 + 10 * λ + λ^2 + 12im * ω) -
        rp^2 * (8im * a * m * (11 + 2 * λ + 6im * ω) + a^2 * (24 + 24 * m^2 + 10 * λ + λ^2 - 28im * ω - 16im * λ * ω + 48 * ω^2)) +
        8im * a^3 * rp * (m * (7 + λ - 6im * ω) - a * ω * (4 + λ)) + 12 * a^5 * ω * (a * ω - 3 * m))
    Xi = Float64(m * (Δφr_func(qr0) + Δφθ_func(qθ0) - qφ0) - ω * (Δtr_func(qr0) + Δtθ_func(qθ0) - qt0) - k * qθ0 - n * qr0)
    prefactor = -im * π * factor * η * exp(im * Xi) / (Γ * κ * Cinc)

    out = Matrix{ComplexF64}(undef, N, K)
    @inbounds for j in 1:K
        sθ = sinθ[j]
        cθ = cosθ[j]
        invsθ = inv_sinθ[j]
        invs2θ = inv_sin2θ[j]
        L1j = L1[j]
        L2j = L2[j]
        L2pj = L2p[j]
        S0j = S0v[j]
        S1j = S1v[j]
        S2j = S2v[j]
        L2Sj = L2S[j]
        L1Spj = L1Sp[j]
        L1L2Sj = L1L2S[j]
        theta_const_j = theta_const[j]
        theta_L1_j = theta_L1[j]
        phaseθ = θphase[j]
        Mpj = termM[j] + uθp[j]
        Mmj = termM[j] + uθm[j]
        expθp = exp(im * phaseθ)
        expθm = exp(-im * phaseθ)

        for i in 1:N
            r = r_vec[i]
            r2 = r * r
            ρ = -1.0 / (r - im * a * cθ)
            ρbar = -1.0 / (r + im * a * cθ)
            invρ = 1.0 / ρ
            Yi = Yv[i]
            Ypi = Ypv[i]
            Xi = Xv[i]
            phase_r = phase[i]
            rphase_i = rphase[i]
            Lpi = Lp[i]
            Lmi = Lm[i]

            L1pL2pS = invρ * L1L2Sj + theta_const_j + theta_L1_j * L1j

            Wll = ρbar * L1pL2pS * r2 * Yi * phase_r * 2.0
            termA = L2Sj + im * a * sθ * (ρ - ρbar) * S0j
            termB = 2.0 * Yi + r * Ypi
            termC = L2Sj * (ρ + ρbar) * r * Yi
            Wlm = 4.0 * (r * ρbar * (termA * termB + termC) / ρ) * phase_r
            termD = Xi / (2.0 * sqrt(r2 + a^2))
            termE = (Yi + 2.0 * r * Ypi) * phase_r
            termF = ρ * r * (2.0 * Yi + r * Ypi) * phase_r
            Wmm = 4.0 * S0j * ρbar * (termD + termE + termF) / ρ

            Lp2 = Lpi * Lpi
            Lm2 = Lmi * Lmi
            Mp2 = Mpj * Mpj
            Mm2 = Mmj * Mmj
            exp_rp = exp(im * rphase_i)
            exp_rm = exp(-im * rphase_i)
            exp_pp = exp_rp * expθp
            exp_pm = exp_rp * expθm
            exp_mp = exp_rm * expθp
            exp_mm = exp_rm * expθm
            Jpp = (Wll * Lp2 + Wlm * Lpi * Mpj + Wmm * Mp2) * exp_pp * inv4π2
            Jpm = (Wll * Lp2 + Wlm * Lpi * Mmj + Wmm * Mm2) * exp_pm * inv4π2
            Jmp = (Wll * Lm2 + Wlm * Lmi * Mpj + Wmm * Mp2) * exp_mp * inv4π2
            Jmm = (Wll * Lm2 + Wlm * Lmi * Mmj + Wmm * Mm2) * exp_mm * inv4π2

            out[i, j] = prefactor * (Jpp + Jpm + Jmp + Jmm)
        end
    end

    return out
end

function integrand_generic_sample_m2(KG_samp, Y_samp, SH_samp, n::Int64, k::Int64)
    # --------------------------
    # 1. Extract Parameters (1D Vectors)
    # --------------------------
    # Radial (N) and Polar (K) sampling vectors
    r_vec = KG_samp["r"]          # Length N (radial samples)
    rs_vec = KG_samp["rs"]        # Precomputed tortoise coordinates (length N)
    θ_vec = KG_samp["θ"]          # Length K (polar samples)
    N = KG_samp["N_sample"]       # Grid size N (radial)
    K = KG_samp["K_sample"]       # Grid size K (polar)
    Γ = KG_samp["Γ"]
    qt0, qr0, qθ0, qφ0 = KG_samp["initialPhases"]
    Δtr_func, Δtθ_func, Δφr_func, Δφθ_func = KG_samp["CrossFunction"]

    # Radial velocity (length N)
    urp = KG_samp["ur_fwd"]       # Forward radial velocity
    urm = KG_samp["ur_rev"]       # Reverse radial velocity

    # Polar velocity (length K)
    uθp = KG_samp["uθ_fwd"]       # Forward polar velocity
    uθm = KG_samp["uθ_rev"]       # Reverse polar velocity

    # Cross functions (length N for Δtr/Δφr, length K for Δtθ/Δφθ)
    Δtr = KG_samp["Δtr"]          # Radial cross function (N)
    Δtθ = KG_samp["Δtθ"]          # Polar cross function (K)
    Δφr = KG_samp["Δφr"]          # Radial cross function (N)
    Δφθ = KG_samp["Δφθ"]          # Polar cross function (K)

    # Orbital constants
    E = KG_samp["E"]              # Energy (scalar)
    Lz = KG_samp["Lz"]            # z-angular momentum (scalar)

    # Spherical harmonics (length K, polar-dependent)
    S0 = SH_samp["S0"]            # Base value (K)
    S1 = SH_samp["S1"]            # 1st θ-derivative (K)
    S2 = SH_samp["S2"]            # 2nd θ-derivative (K)

    # Radial functions (length N, radial-dependent)
    Y = Y_samp["Y"]               # Radial function (N)
    Yp = Y_samp["Yp"]             # 1st r-derivative (N)
    X = Y_samp["X"]               # Auxiliary radial function (N)
    Binc = Y_samp["Binc"]

    # Constant parameters
    params = Y_samp["params"]
    s = params.s                  # Spin parameter (scalar)
    l = params.l                  # Angular quantum number (scalar)
    m = params.m                  # Magnetic quantum number (scalar)
    a = params.a                  # Black hole spin (scalar)
    ω = params.omega              # Frequency (scalar)
    λ = params.lambda             # Eigenvalue (scalar)

    # Black hole horizons (scalar)
    rp = 1.0 + sqrt(1.0 - a^2)    # Outer horizon
    rm = 1.0 - sqrt(1.0 - a^2)    # Inner horizon

    # --------------------------
    # 2. Generate 2D Grid Matrices (Vectorization Core)
    # --------------------------
    # Create N×K grids for r, rs, θ (broadcast 1D vectors to 2D)
    r_grid = reshape(r_vec, :, 1)  # N×1: broadcast across columns
    θ_grid = reshape(θ_vec, 1, :)   # 1×K: broadcast across rows

    # Precompute θ-dependent grid matrices (K→N×K broadcast)
    sinθ_grid = sin.(θ_grid)
    cosθ_grid = cos.(θ_grid)
    cotθ_grid = cot.(θ_grid)
    inv_sinθ_grid = 1.0 ./ sinθ_grid
    inv_sin2θ_grid = 1.0 ./ (sinθ_grid .^ 2)  # Precompute 1/sin²θ

    # Precompute r-dependent grid matrices (N→N×K broadcast)
    r2_grid = r_grid .^ 2
    Δ_grid = r2_grid .- 2.0 .* r_grid .+ a^2  # Kerr Δ(r) grid

    # --------------------------
    # 3. Compute Complex Parameters (Vectorized)
    # --------------------------
    # ρ and ρbar (N×K, complex-valued)
    ρ_grid = -1.0 ./ (r_grid .- im .* a .* cosθ_grid)
    ρbar_grid = -1.0 ./ (r_grid .+ im .* a .* cosθ_grid)
    inv_ρ_grid = 1.0 ./ ρ_grid  # Precompute 1/ρ

    # --------------------------
    # 4. L-series Operators (Vectorized)
    # --------------------------
    L1_grid = (-m) .* inv_sinθ_grid .+ a .* ω .* sinθ_grid .+ cotθ_grid
    L2_grid = (-m) .* inv_sinθ_grid .+ a .* ω .* sinθ_grid .+ 2.0 .* cotθ_grid
    L2p_grid = m .* cosθ_grid .* inv_sin2θ_grid .+ a .* ω .* cosθ_grid .- 2.0 .* inv_sin2θ_grid

    # --------------------------
    # 5. Spherical Harmonic Combinations (Vectorized)
    # --------------------------
    S0_grid = reshape(S0, 1, :)   # 1×K: broadcast across rows
    S1_grid = reshape(S1, 1, :)   # 1×K: broadcast across rows
    S2_grid = reshape(S2, 1, :)   # 1×K: broadcast across rows

    L2S_grid = S1_grid .+ L2_grid .* S0_grid
    L1Sp_grid = S2_grid .+ L1_grid .* S1_grid
    L1L2S_grid = L1Sp_grid .+ L2p_grid .* S0_grid .+ L2_grid .* S1_grid .+ L1_grid .* L2_grid .* S0_grid

    # L1pL2pS (vectorized complex computation)
    term1 = inv_ρ_grid .* L1L2S_grid
    term2 = 3im .* a .* sinθ_grid .* L1_grid .* S0_grid
    term3 = 3im .* a .* cosθ_grid .* S0_grid
    term4 = 2im .* a .* sinθ_grid .* S1_grid
    term5 = -1im .* a .* sinθ_grid .* L2_grid .* S0_grid
    L1pL2pS_grid = term1 .+ term2 .+ term3 .+ term4 .+ term5

    # --------------------------
    # 6. Phase Terms (Vectorized, using precomputed rs_vec)
    # --------------------------
    # Radial phase (phase[i,j]): use precomputed rs_vec instead of rstar_from_r
    log_term_vec = log.((r_vec .- rp) ./ (r_vec .- rm))  # 1D (N)
    phase_vec = exp.(im .* ω .* rs_vec .- im .* a .* m .* log_term_vec ./ (2.0 .* sqrt(1.0 .- a^2)))  # 1D
    phase_grid = reshape(phase_vec, :, 1)  # N×1: broadcast across columns

    # rphase/θphase (1D, then broadcast to N×K for J terms)
    qr = range(0.0, π, length=N)
    qθ = range(0.0, π, length=K)
    rphase_vec = ω .* Δtr .- m .* Δφr .+ n .* qr  # 1D (N)
    θphase_vec = ω .* Δtθ .- m .* Δφθ .+ k .* qθ  # 1D (K)
    rphase_grid = reshape(rphase_vec, :, 1)       # N×1
    θphase_grid = reshape(θphase_vec, 1, :)       # 1×K

    # --------------------------
    # 7. W-series Quantities (Vectorized)
    # --------------------------
    Y_grid = reshape(Y, :, 1)
    Yp_grid = reshape(Yp, :, 1)
    X_grid = reshape(X, :, 1)

    # Wnn: use precomputed phase_grid
    Wnn_grid = ρbar_grid .* L1pL2pS_grid .* r2_grid .* Y_grid .* phase_grid ./ 2.0

    # Wnmbar
    term_A = L2S_grid .+ im .* a .* sinθ_grid .* (ρ_grid .- ρbar_grid) .* S0_grid
    term_B = 2.0 .* Y_grid .+ r_grid .* Yp_grid
    term_C = L2S_grid .* (ρ_grid .+ ρbar_grid) .* r_grid .* Y_grid
    Wnmbar_grid = (-r_grid .* ρbar_grid .* (term_A .* term_B .+ term_C) ./ ρ_grid) .* phase_grid

    # Wmbarmbar
    term_D = X_grid ./ (2.0 .* sqrt.(r2_grid .+ a^2))
    term_E = (Y_grid .+ 2.0 .* r_grid .* Yp_grid) .* phase_grid
    term_F = ρ_grid .* r_grid .* (2.0 .* Y_grid .+ r_grid .* Yp_grid) .* phase_grid
    Wmbarmbar_grid = S0_grid .* ρbar_grid .* (term_D .+ term_E .+ term_F) ./ ρ_grid

    # --------------------------
    # 8. N/Mbar Terms (Vectorized)
    # --------------------------
    urp_grid = reshape(urp, :, 1)
    urm_grid = reshape(urm, :, 1)
    uθp_grid = reshape(uθp, 1, :)
    uθm_grid = reshape(uθm, 1, :)

    # Np/Nm (using precomputed Δ_grid)
    numerator_N = E .* (r2_grid .+ a^2) .- a .* Lz
    Np_grid = (numerator_N .+ urp_grid) ./ Δ_grid
    Nm_grid = (numerator_N .+ urm_grid) ./ Δ_grid

    # Mbarp/Mbarm
    term_M = im .* sinθ_grid .* (a .* E .- Lz .* inv_sin2θ_grid)
    Mbarp_grid = term_M .+ uθp_grid
    Mbarm_grid = term_M .+ uθm_grid

    # --------------------------
    # 9. J-series Terms (Vectorized)
    # --------------------------
    Np_sq = Np_grid .^ 2
    Nm_sq = Nm_grid .^ 2
    Mbarp_sq = Mbarp_grid .^ 2
    Mbarm_sq = Mbarm_grid .^ 2

    exp_pp = exp.(im .* (rphase_grid .+ θphase_grid))
    exp_pm = exp.(im .* (rphase_grid .- θphase_grid))
    exp_mp = exp.(im .* (-rphase_grid .+ θphase_grid))
    exp_mm = exp.(im .* (-rphase_grid .- θphase_grid))

    Jpp_grid = (Wnn_grid .* Np_sq .+ Wnmbar_grid .* Np_grid .* Mbarp_grid .+ Wmbarmbar_grid .* Mbarp_sq) .* exp_pp
    Jpm_grid = (Wnn_grid .* Np_sq .+ Wnmbar_grid .* Np_grid .* Mbarm_grid .+ Wmbarmbar_grid .* Mbarm_sq) .* exp_pm
    Jmp_grid = (Wnn_grid .* Nm_sq .+ Wnmbar_grid .* Nm_grid .* Mbarp_grid .+ Wmbarmbar_grid .* Mbarp_sq) .* exp_mp
    Jmm_grid = (Wnn_grid .* Nm_sq .+ Wnmbar_grid .* Nm_grid .* Mbarm_grid .+ Wmbarmbar_grid .* Mbarm_sq) .* exp_mm

    # --------------------------
    # 10. Sum J Terms (Final Result)
    # --------------------------
    Jlmnk = (Jpp_grid .+ Jpm_grid .+ Jmp_grid .+ Jmm_grid) ./ (4 .* π .^2)

    Xi = m * (Δφr_func(qr0) + Δφθ_func(qθ0) - qφ0) - ω * (Δtr_func(qr0) + Δtθ_func(qθ0) - qt0) - k * qθ0 - n * qr0

    return 4im .* π .* ω .* Jlmnk .* exp(im * Xi) ./ (Γ .* Binc)
end

function integrand_generic_sample_p2(KG_samp, Y_samp, SH_samp, n::Int64, k::Int64)
    # --------------------------
    # 1. Extract Parameters (1D Vectors)
    # --------------------------
    # Radial (N) and Polar (K) sampling vectors
    r_vec = KG_samp["r"]          # Length N (radial samples)
    rs_vec = KG_samp["rs"]        # Precomputed tortoise coordinates (length N)
    θ_vec = KG_samp["θ"]          # Length K (polar samples)
    N = KG_samp["N_sample"]       # Grid size N (radial)
    K = KG_samp["K_sample"]       # Grid size K (polar)
    Γ = KG_samp["Γ"]
    qt0, qr0, qθ0, qφ0 = KG_samp["initialPhases"]
    Δtr_func, Δtθ_func, Δφr_func, Δφθ_func = KG_samp["CrossFunction"]

    # Radial velocity (length N)
    urp = KG_samp["ur_fwd"]       # Forward radial velocity
    urm = KG_samp["ur_rev"]       # Reverse radial velocity

    # Polar velocity (length K)
    uθp = KG_samp["uθ_fwd"]       # Forward polar velocity
    uθm = KG_samp["uθ_rev"]       # Reverse polar velocity

    # Cross functions (length N for Δtr/Δφr, length K for Δtθ/Δφθ)
    Δtr = KG_samp["Δtr"]          # Radial cross function (N)
    Δtθ = KG_samp["Δtθ"]          # Polar cross function (K)
    Δφr = KG_samp["Δφr"]          # Radial cross function (N)
    Δφθ = KG_samp["Δφθ"]          # Polar cross function (K)

    # Orbital constants
    E = KG_samp["E"]              # Energy (scalar)
    Lz = KG_samp["Lz"]            # z-angular momentum (scalar)

    # Spherical harmonics (length K, polar-dependent)
    S0 = SH_samp["S0"]            # Base value (K)
    S1 = SH_samp["S1"]            # 1st θ-derivative (K)
    S2 = SH_samp["S2"]            # 2nd θ-derivative (K)

    # Radial functions (length N, radial-dependent)
    Y = Y_samp["Y"]               # Radial function (N)
    Yp = Y_samp["Yp"]             # 1st r-derivative (N)
    X = Y_samp["X"]               # Auxiliary radial function (N)
    Cinc = Y_samp["Cinc"]

    # Constant parameters
    params = Y_samp["params"]
    s = params.s                  # Spin parameter (scalar)
    l = params.l                  # Angular quantum number (scalar)
    m = params.m                  # Magnetic quantum number (scalar)
    a = params.a                  # Black hole spin (scalar)
    ω = params.omega              # Frequency (scalar)
    λ = params.lambda             # Eigenvalue (scalar)

    # Black hole horizons (scalar)
    rp = 1.0 + sqrt(1.0 - a^2)    # Outer horizon
    rm = 1.0 - sqrt(1.0 - a^2)    # Inner horizon

    # --------------------------
    # 2. Generate 2D Grid Matrices (Vectorization Core)
    # --------------------------
    # Create N×K grids for r, rs, θ (broadcast 1D vectors to 2D)
    r_grid = reshape(r_vec, :, 1)  # N×1: broadcast across columns
    θ_grid = reshape(θ_vec, 1, :)  # 1×K: broadcast across rows

    # Precompute θ-dependent grid matrices (K→N×K broadcast)
    sinθ_grid = sin.(θ_grid)
    cosθ_grid = cos.(θ_grid)
    cotθ_grid = cot.(θ_grid)
    inv_sinθ_grid = 1.0 ./ sinθ_grid
    inv_sin2θ_grid = 1.0 ./ (sinθ_grid .^ 2)  # Precompute 1/sin²θ

    # Precompute r-dependent grid matrices (N→N×K broadcast)
    r2_grid = r_grid .^ 2
    Δ_grid = r2_grid .- 2.0 .* r_grid .+ a^2  # Kerr Δ(r) grid

    # --------------------------
    # 3. Compute Complex Parameters (Vectorized)
    # --------------------------
    # ρ and ρbar (N×K, complex-valued)
    ρ_grid = -1.0 ./ (r_grid .- im .* a .* cosθ_grid)
    ρbar_grid = -1.0 ./ (r_grid .+ im .* a .* cosθ_grid)
    inv_ρ_grid = 1.0 ./ ρ_grid  # Precompute 1/ρ

    # --------------------------
    # 4. L-series Operators (Vectorized)
    # --------------------------
    L1_grid = m .* inv_sinθ_grid .- a .* ω .* sinθ_grid .+ cotθ_grid
    L2_grid = m .* inv_sinθ_grid .- a .* ω .* sinθ_grid .+ 2.0 .* cotθ_grid
    L2p_grid = - m .* cosθ_grid .* inv_sin2θ_grid .- a .* ω .* cosθ_grid .- 2.0 .* inv_sin2θ_grid

    # --------------------------
    # 5. Spherical Harmonic Combinations (Vectorized)
    # --------------------------
    S0_grid = reshape(S0, 1, :)   # 1×K: broadcast across rows
    S1_grid = reshape(S1, 1, :)   # 1×K: broadcast across rows
    S2_grid = reshape(S2, 1, :)   # 1×K: broadcast across rows

    L2S_grid = S1_grid .+ L2_grid .* S0_grid
    L1Sp_grid = S2_grid .+ L1_grid .* S1_grid
    L1L2S_grid = L1Sp_grid .+ L2p_grid .* S0_grid .+ L2_grid .* S1_grid .+ L1_grid .* L2_grid .* S0_grid

    # L1pL2pS (vectorized complex computation)
    term1 = inv_ρ_grid .* L1L2S_grid
    term2 = 3im .* a .* sinθ_grid .* L1_grid .* S0_grid
    term3 = 3im .* a .* cosθ_grid .* S0_grid
    term4 = 2im .* a .* sinθ_grid .* S1_grid
    term5 = -1im .* a .* sinθ_grid .* L2_grid .* S0_grid
    L1pL2pS_grid = term1 .+ term2 .+ term3 .+ term4 .+ term5

    # --------------------------
    # 6. Phase Terms (Vectorized, using precomputed rs_vec)
    # --------------------------
    # Radial phase (phase[i,j]): use precomputed rs_vec instead of rstar_from_r
    log_term_vec = log.((r_vec .- rp) ./ (r_vec .- rm))  # 1D (N)
    phase_vec = exp.(.- im .* ω .* rs_vec .+ im .* a .* m .* log_term_vec ./ (2.0 .* sqrt(1.0 .- a^2)))  # 1D
    phase_grid = reshape(phase_vec, :, 1)  # N×1: broadcast across columns

    # rphase/θphase (1D, then broadcast to N×K for J terms)
    qr = range(0.0, π, length=N)
    qθ = range(0.0, π, length=K)
    rphase_vec = ω .* Δtr .- m .* Δφr .+ n .* qr  # 1D (N)
    θphase_vec = ω .* Δtθ .- m .* Δφθ .+ k .* qθ  # 1D (K)
    rphase_grid = reshape(rphase_vec, :, 1)       # N×1
    θphase_grid = reshape(θphase_vec, 1, :)       # 1×K

    # --------------------------
    # 7. W-series Quantities (Vectorized)
    # --------------------------
    Y_grid = reshape(Y, :, 1)
    Yp_grid = reshape(Yp, :, 1)
    X_grid = reshape(X, :, 1)

    # Wll: use precomputed phase_grid
    Wll_grid = ρbar_grid .* L1pL2pS_grid .* r2_grid .* Y_grid .* phase_grid .* 2.0

    # Wlm
    term_A = L2S_grid .+ im .* a .* sinθ_grid .* (ρ_grid .- ρbar_grid) .* S0_grid
    term_B = 2.0 .* Y_grid .+ r_grid .* Yp_grid
    term_C = L2S_grid .* (ρ_grid .+ ρbar_grid) .* r_grid .* Y_grid
    Wlm_grid = 4.0 .* (r_grid .* ρbar_grid .* (term_A .* term_B .+ term_C) ./ ρ_grid) .* phase_grid

    # Wmm
    term_D = X_grid ./ (2.0 .* sqrt.(r2_grid .+ a^2))
    term_E = (Y_grid .+ 2.0 .* r_grid .* Yp_grid) .* phase_grid
    term_F = ρ_grid .* r_grid .* (2.0 .* Y_grid .+ r_grid .* Yp_grid) .* phase_grid
    Wmm_grid = 4.0 .* S0_grid .* ρbar_grid .* (term_D .+ term_E .+ term_F) ./ ρ_grid

    # --------------------------
    # 8. L/M Terms (Vectorized)
    # --------------------------
    urp_grid = reshape(urp, :, 1)
    urm_grid = reshape(urm, :, 1)
    uθp_grid = reshape(uθp, 1, :)
    uθm_grid = reshape(uθm, 1, :)

    # Lp/Lm (using precomputed Δ_grid)
    numerator_L = E .* (r2_grid .+ a^2) .- a .* Lz
    Lp_grid = (numerator_L .- urp_grid) ./ Δ_grid
    Lm_grid = (numerator_L .- urm_grid) ./ Δ_grid

    # Mp/Mm
    term_M = .- im .* sinθ_grid .* (a .* E .- Lz .* inv_sin2θ_grid)
    Mp_grid = term_M .+ uθp_grid
    Mm_grid = term_M .+ uθm_grid

    # --------------------------
    # 9. J-series Terms (Vectorized)
    # --------------------------
    Lp_sq = Lp_grid .^ 2
    Lm_sq = Lm_grid .^ 2
    Mp_sq = Mp_grid .^ 2
    Mm_sq = Mm_grid .^ 2

    exp_pp = exp.(im .* (rphase_grid .+ θphase_grid))
    exp_pm = exp.(im .* (rphase_grid .- θphase_grid))
    exp_mp = exp.(im .* (-rphase_grid .+ θphase_grid))
    exp_mm = exp.(im .* (-rphase_grid .- θphase_grid))

    Jpp_grid = (Wll_grid .* Lp_sq .+ Wlm_grid .* Lp_grid .* Mp_grid .+ Wmm_grid .* Mp_sq) .* exp_pp
    Jpm_grid = (Wll_grid .* Lp_sq .+ Wlm_grid .* Lp_grid .* Mm_grid .+ Wmm_grid .* Mm_sq) .* exp_pm
    Jmp_grid = (Wll_grid .* Lm_sq .+ Wlm_grid .* Lm_grid .* Mp_grid .+ Wmm_grid .* Mp_sq) .* exp_mp
    Jmm_grid = (Wll_grid .* Lm_sq .+ Wlm_grid .* Lm_grid .* Mm_grid .+ Wmm_grid .* Mm_sq) .* exp_mm

    # --------------------------
    # 10. Sum J Terms (Final Result)
    # --------------------------
    Jlmnk = (Jpp_grid .+ Jpm_grid .+ Jmp_grid .+ Jmm_grid) ./ (4 .* π .^2)

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

    Xi = m * (Δφr_func(qr0) + Δφθ_func(qθ0) - qφ0) - ω * (Δtr_func(qr0) + Δtθ_func(qθ0) - qt0) - k * qθ0 - n * qr0

    return .- im .* π .* factor .* Jlmnk .* exp(im * Xi) .* η ./ (Γ .* κ .* Cinc)
end

function integrand_generic_sample_cheby_m2(KG_samp, Y_samp, SH_samp, n::Int64, k::Int64)
    # --------------------------
    # 1. Extract Parameters (1D Vectors)
    # --------------------------
    # Radial (N) and Polar (K) sampling vectors
    r_vec = KG_samp["r"]          # Length N (radial samples)
    rs_vec = KG_samp["rs"]        # Precomputed tortoise coordinates (length N)
    θ_vec = KG_samp["θ"]          # Length K (polar samples)
    N = KG_samp["N_sample"]       # Grid size N (radial)
    K = KG_samp["K_sample"]       # Grid size K (polar)
    Γ = KG_samp["Γ"]
    qt0, qr0, qθ0, qφ0 = KG_samp["initialPhases"]
    Δtr_func, Δtθ_func, Δφr_func, Δφθ_func = KG_samp["CrossFunction"]
    dtr = KG_samp["dtr"]
    dφr = KG_samp["dφr"]
    dtθ = KG_samp["dtθ"]
    dφθ = KG_samp["dφθ"]

    # Radial velocity (length N)
    urp = KG_samp["ur_fwd"]       # Forward radial velocity
    urm = KG_samp["ur_rev"]       # Reverse radial velocity

    # Polar velocity (length K)
    uθp = KG_samp["uθ_fwd"]       # Forward polar velocity
    uθm = KG_samp["uθ_rev"]       # Reverse polar velocity

    # Cross functions (length N for Δtr/Δφr, length K for Δtθ/Δφθ)
    Δtr = KG_samp["Δtr"]          # Radial cross function (N)
    Δtθ = KG_samp["Δtθ"]          # Polar cross function (K)
    Δφr = KG_samp["Δφr"]          # Radial cross function (N)
    Δφθ = KG_samp["Δφθ"]          # Polar cross function (K)

    # Orbital constants
    E = KG_samp["E"]              # Energy (scalar)
    Lz = KG_samp["Lz"]            # z-angular momentum (scalar)

    # Spherical harmonics (length K, polar-dependent)
    S0 = SH_samp["S0"]            # Base value (K)
    S1 = SH_samp["S1"]            # 1st θ-derivative (K)
    S2 = SH_samp["S2"]            # 2nd θ-derivative (K)

    # Radial functions (length N, radial-dependent)
    Y = Y_samp["Y"]               # Radial function (N)
    Yp = Y_samp["Yp"]             # 1st r-derivative (N)
    X = Y_samp["X"]               # Auxiliary radial function (N)
    Binc = Y_samp["Binc"]

    # Constant parameters
    params = Y_samp["params"]
    s = params.s                  # Spin parameter (scalar)
    l = params.l                  # Angular quantum number (scalar)
    m = params.m                  # Magnetic quantum number (scalar)
    a = params.a                  # Black hole spin (scalar)
    ω = params.omega              # Frequency (scalar)
    λ = params.lambda             # Eigenvalue (scalar)

    # Black hole horizons (scalar)
    rp = 1.0 + sqrt(1.0 - a^2)    # Outer horizon
    rm = 1.0 - sqrt(1.0 - a^2)    # Inner horizon

    # --------------------------
    # 2. Generate 2D Grid Matrices (Vectorization Core)
    # --------------------------
    # Create N×K grids for r, rs, θ (broadcast 1D vectors to 2D)
    r_grid = reshape(r_vec, :, 1)  # N×1: broadcast across columns
    θ_grid = reshape(θ_vec, 1, :)  # 1×K: broadcast across rows

    # Precompute θ-dependent grid matrices (K→N×K broadcast)
    sinθ_grid = sin.(θ_grid)
    cosθ_grid = cos.(θ_grid)
    cotθ_grid = cot.(θ_grid)
    inv_sinθ_grid = 1.0 ./ sinθ_grid
    inv_sin2θ_grid = 1.0 ./ (sinθ_grid .^ 2)  # Precompute 1/sin²θ

    # Precompute r-dependent grid matrices (N→N×K broadcast)
    r2_grid = r_grid .^ 2
    Δ_grid = r2_grid .- 2.0 .* r_grid .+ a^2  # Kerr Δ(r) grid

    # --------------------------
    # 3. Compute Complex Parameters (Vectorized)
    # --------------------------
    # ρ and ρbar (N×K, complex-valued)
    ρ_grid = -1.0 ./ (r_grid .- im .* a .* cosθ_grid)
    ρbar_grid = -1.0 ./ (r_grid .+ im .* a .* cosθ_grid)
    inv_ρ_grid = 1.0 ./ ρ_grid  # Precompute 1/ρ

    # --------------------------
    # 4. L-series Operators (Vectorized)
    # --------------------------
    L1_grid = (-m) .* inv_sinθ_grid .+ a .* ω .* sinθ_grid .+ cotθ_grid
    L2_grid = (-m) .* inv_sinθ_grid .+ a .* ω .* sinθ_grid .+ 2.0 .* cotθ_grid
    L2p_grid = m .* cosθ_grid .* inv_sin2θ_grid .+ a .* ω .* cosθ_grid .- 2.0 .* inv_sin2θ_grid

    # --------------------------
    # 5. Spherical Harmonic Combinations (Vectorized)
    # --------------------------
    S0_grid = reshape(S0, 1, :)   # 1×K: broadcast across rows
    S1_grid = reshape(S1, 1, :)   # 1×K: broadcast across rows
    S2_grid = reshape(S2, 1, :)   # 1×K: broadcast across rows

    L2S_grid = S1_grid .+ L2_grid .* S0_grid
    L1Sp_grid = S2_grid .+ L1_grid .* S1_grid
    L1L2S_grid = L1Sp_grid .+ L2p_grid .* S0_grid .+ L2_grid .* S1_grid .+ L1_grid .* L2_grid .* S0_grid

    # L1pL2pS (vectorized complex computation)
    term1 = inv_ρ_grid .* L1L2S_grid
    term2 = 3im .* a .* sinθ_grid .* L1_grid .* S0_grid
    term3 = 3im .* a .* cosθ_grid .* S0_grid
    term4 = 2im .* a .* sinθ_grid .* S1_grid
    term5 = -1im .* a .* sinθ_grid .* L2_grid .* S0_grid
    L1pL2pS_grid = term1 .+ term2 .+ term3 .+ term4 .+ term5

    # --------------------------
    # 6. Phase Terms (Vectorized, using precomputed rs_vec)
    # --------------------------
    # Radial phase (phase[i,j]): use precomputed rs_vec instead of rstar_from_r
    log_term_vec = log.((r_vec .- rp) ./ (r_vec .- rm))  # 1D (N)
    phase_vec = exp.(im .* ω .* rs_vec .- im .* a .* m .* log_term_vec ./ (2.0 .* sqrt(1.0 .- a^2)))  # 1D
    phase_grid = reshape(phase_vec, :, 1)  # N×1: broadcast across columns

    # rphase/θphase (1D, then broadcast to N×K for J terms)
    qr = range(0.0, π, length=N)
    qθ = range(0.0, π, length=K)
    drphase = ω .* (dtr) .- m .* (dφr) .+ n 
    rphaseL = ω * (Δtr[1]) - m * (Δφr[1])
    rphaseR = ω * (Δtr[end]) - m * (Δφr[end]) + n * 1.0pi
    dθphase = ω .* dtθ .- m .* dφθ .+ k
    θphaseL = ω * Δtθ[1] - m * Δφθ[1]
    θphaseR = ω * Δtθ[end] - m * Δφθ[end] + k * 1.0pi

    # --------------------------
    # 7. W-series Quantities (Vectorized)
    # --------------------------
    Y_grid = reshape(Y, :, 1)
    Yp_grid = reshape(Yp, :, 1)
    X_grid = reshape(X, :, 1)

    # Wnn: use precomputed phase_grid
    Wnn_grid = ρbar_grid .* L1pL2pS_grid .* r2_grid .* Y_grid .* phase_grid ./ 2.0

    # Wnmbar
    term_A = L2S_grid .+ im .* a .* sinθ_grid .* (ρ_grid .- ρbar_grid) .* S0_grid
    term_B = 2.0 .* Y_grid .+ r_grid .* Yp_grid
    term_C = L2S_grid .* (ρ_grid .+ ρbar_grid) .* r_grid .* Y_grid
    Wnmbar_grid = (-r_grid .* ρbar_grid .* (term_A .* term_B .+ term_C) ./ ρ_grid) .* phase_grid

    # Wmbarmbar
    term_D = X_grid ./ (2.0 .* sqrt.(r2_grid .+ a^2))
    term_E = (Y_grid .+ 2.0 .* r_grid .* Yp_grid) .* phase_grid
    term_F = ρ_grid .* r_grid .* (2.0 .* Y_grid .+ r_grid .* Yp_grid) .* phase_grid
    Wmbarmbar_grid = S0_grid .* ρbar_grid .* (term_D .+ term_E .+ term_F) ./ ρ_grid

    # --------------------------
    # 8. N/Mbar Terms (Vectorized)
    # --------------------------
    urp_grid = reshape(urp, :, 1)
    urm_grid = reshape(urm, :, 1)
    uθp_grid = reshape(uθp, 1, :)
    uθm_grid = reshape(uθm, 1, :)

    # Np/Nm (using precomputed Δ_grid)
    numerator_N = E .* (r2_grid .+ a^2) .- a .* Lz
    Np_grid = (numerator_N .+ urp_grid) ./ Δ_grid
    Nm_grid = (numerator_N .+ urm_grid) ./ Δ_grid

    # Mbarp/Mbarm
    term_M = im .* sinθ_grid .* (a .* E .- Lz .* inv_sin2θ_grid)
    Mbarp_grid = term_M .+ uθp_grid
    Mbarm_grid = term_M .+ uθm_grid

    # --------------------------
    # 9. J-series Terms (Vectorized)
    # --------------------------
    Np_sq = Np_grid .^ 2
    Nm_sq = Nm_grid .^ 2
    Mbarp_sq = Mbarp_grid .^ 2
    Mbarm_sq = Mbarm_grid .^ 2

    Jpp_grid = (Wnn_grid .* Np_sq .+ Wnmbar_grid .* Np_grid .* Mbarp_grid .+ Wmbarmbar_grid .* Mbarp_sq) ./ (4 .* π .^2)
    Jpm_grid = (Wnn_grid .* Np_sq .+ Wnmbar_grid .* Np_grid .* Mbarm_grid .+ Wmbarmbar_grid .* Mbarm_sq) ./ (4 .* π .^2)
    Jmp_grid = (Wnn_grid .* Nm_sq .+ Wnmbar_grid .* Nm_grid .* Mbarp_grid .+ Wmbarmbar_grid .* Mbarp_sq) ./ (4 .* π .^2)
    Jmm_grid = (Wnn_grid .* Nm_sq .+ Wnmbar_grid .* Nm_grid .* Mbarm_grid .+ Wmbarmbar_grid .* Mbarm_sq) ./ (4 .* π .^2)

    Xi = m * (Δφr_func(qr0) + Δφθ_func(qθ0) - qφ0) - ω * (Δtr_func(qr0) + Δtθ_func(qθ0) - qt0) - k * qθ0 - n * qr0
    prefactor = 4im .* π .* ω .* exp(im * Xi) ./ (Γ .* Binc)

    return Jpp_grid, Jpm_grid, Jmp_grid, Jmm_grid, drphase, dθphase, rphaseL, rphaseR, θphaseL, θphaseR, prefactor
end

function integrand_generic_sample_cheby_p2(KG_samp, Y_samp, SH_samp, n::Int64, k::Int64)
    # --------------------------
    # 1. Extract Parameters (1D Vectors)
    # --------------------------
    # Radial (N) and Polar (K) sampling vectors
    r_vec = KG_samp["r"]          # Length N (radial samples)
    rs_vec = KG_samp["rs"]        # Precomputed tortoise coordinates (length N)
    θ_vec = KG_samp["θ"]          # Length K (polar samples)
    N = KG_samp["N_sample"]       # Grid size N (radial)
    K = KG_samp["K_sample"]       # Grid size K (polar)
    Γ = KG_samp["Γ"]
    qt0, qr0, qθ0, qφ0 = KG_samp["initialPhases"]
    Δtr_func, Δtθ_func, Δφr_func, Δφθ_func = KG_samp["CrossFunction"]
    dtr = KG_samp["dtr"]
    dφr = KG_samp["dφr"]
    dtθ = KG_samp["dtθ"]
    dφθ = KG_samp["dφθ"]

    # Radial velocity (length N)
    urp = KG_samp["ur_fwd"]       # Forward radial velocity
    urm = KG_samp["ur_rev"]       # Reverse radial velocity

    # Polar velocity (length K)
    uθp = KG_samp["uθ_fwd"]       # Forward polar velocity
    uθm = KG_samp["uθ_rev"]       # Reverse polar velocity

    # Cross functions (length N for Δtr/Δφr, length K for Δtθ/Δφθ)
    Δtr = KG_samp["Δtr"]          # Radial cross function (N)
    Δtθ = KG_samp["Δtθ"]          # Polar cross function (K)
    Δφr = KG_samp["Δφr"]          # Radial cross function (N)
    Δφθ = KG_samp["Δφθ"]          # Polar cross function (K)

    # Orbital constants
    E = KG_samp["E"]              # Energy (scalar)
    Lz = KG_samp["Lz"]            # z-angular momentum (scalar)

    # Spherical harmonics (length K, polar-dependent)
    S0 = SH_samp["S0"]            # Base value (K)
    S1 = SH_samp["S1"]            # 1st θ-derivative (K)
    S2 = SH_samp["S2"]            # 2nd θ-derivative (K)

    # Radial functions (length N, radial-dependent)
    Y = Y_samp["Y"]               # Radial function (N)
    Yp = Y_samp["Yp"]             # 1st r-derivative (N)
    X = Y_samp["X"]               # Auxiliary radial function (N)
    Cinc = Y_samp["Cinc"]

    # Constant parameters
    params = Y_samp["params"]
    s = params.s                  # Spin parameter (scalar)
    l = params.l                  # Angular quantum number (scalar)
    m = params.m                  # Magnetic quantum number (scalar)
    a = params.a                  # Black hole spin (scalar)
    ω = params.omega              # Frequency (scalar)
    λ = params.lambda             # Eigenvalue (scalar)

    # Black hole horizons (scalar)
    rp = 1.0 + sqrt(1.0 - a^2)    # Outer horizon
    rm = 1.0 - sqrt(1.0 - a^2)    # Inner horizon

    # --------------------------
    # 2. Generate 2D Grid Matrices (Vectorization Core)
    # --------------------------
    # Create N×K grids for r, rs, θ (broadcast 1D vectors to 2D)
    r_grid = reshape(r_vec, :, 1)  # N×1: broadcast across columns
    θ_grid = reshape(θ_vec, 1, :)  # 1×K: broadcast across rows

    # Precompute θ-dependent grid matrices (K→N×K broadcast)
    sinθ_grid = sin.(θ_grid)
    cosθ_grid = cos.(θ_grid)
    cotθ_grid = cot.(θ_grid)
    inv_sinθ_grid = 1.0 ./ sinθ_grid
    inv_sin2θ_grid = 1.0 ./ (sinθ_grid .^ 2)  # Precompute 1/sin²θ

    # Precompute r-dependent grid matrices (N→N×K broadcast)
    r2_grid = r_grid .^ 2
    Δ_grid = r2_grid .- 2.0 .* r_grid .+ a^2  # Kerr Δ(r) grid

    # --------------------------
    # 3. Compute Complex Parameters (Vectorized)
    # --------------------------
    # ρ and ρbar (N×K, complex-valued)
    ρ_grid = -1.0 ./ (r_grid .- im .* a .* cosθ_grid)
    ρbar_grid = -1.0 ./ (r_grid .+ im .* a .* cosθ_grid)
    inv_ρ_grid = 1.0 ./ ρ_grid  # Precompute 1/ρ

    # --------------------------
    # 4. L-series Operators (Vectorized)
    # --------------------------
    L1_grid = m .* inv_sinθ_grid .- a .* ω .* sinθ_grid .+ cotθ_grid
    L2_grid = m .* inv_sinθ_grid .- a .* ω .* sinθ_grid .+ 2.0 .* cotθ_grid
    L2p_grid = - m .* cosθ_grid .* inv_sin2θ_grid .- a .* ω .* cosθ_grid .- 2.0 .* inv_sin2θ_grid

    # --------------------------
    # 5. Spherical Harmonic Combinations (Vectorized)
    # --------------------------
    S0_grid = reshape(S0, 1, :)   # 1×K: broadcast across rows
    S1_grid = reshape(S1, 1, :)   # 1×K: broadcast across rows
    S2_grid = reshape(S2, 1, :)   # 1×K: broadcast across rows

    L2S_grid = S1_grid .+ L2_grid .* S0_grid
    L1Sp_grid = S2_grid .+ L1_grid .* S1_grid
    L1L2S_grid = L1Sp_grid .+ L2p_grid .* S0_grid .+ L2_grid .* S1_grid .+ L1_grid .* L2_grid .* S0_grid

    # L1pL2pS (vectorized complex computation)
    term1 = inv_ρ_grid .* L1L2S_grid
    term2 = 3im .* a .* sinθ_grid .* L1_grid .* S0_grid
    term3 = 3im .* a .* cosθ_grid .* S0_grid
    term4 = 2im .* a .* sinθ_grid .* S1_grid
    term5 = -1im .* a .* sinθ_grid .* L2_grid .* S0_grid
    L1pL2pS_grid = term1 .+ term2 .+ term3 .+ term4 .+ term5

    # --------------------------
    # 6. Phase Terms (Vectorized, using precomputed rs_vec)
    # --------------------------
    # Radial phase (phase[i,j]): use precomputed rs_vec instead of rstar_from_r
    log_term_vec = log.((r_vec .- rp) ./ (r_vec .- rm))  # 1D (N)
    phase_vec = exp.(.- im .* ω .* rs_vec .+ im .* a .* m .* log_term_vec ./ (2.0 .* sqrt(1.0 .- a^2)))  # 1D
    phase_grid = reshape(phase_vec, :, 1)  # N×1: broadcast across columns

    # rphase/θphase (1D, then broadcast to N×K for J terms)
    qr = range(0.0, π, length=N)
    qθ = range(0.0, π, length=K)
    drphase = ω .* (dtr) .- m .* (dφr) .+ n 
    rphaseL = ω * (Δtr[1]) - m * (Δφr[1])
    rphaseR = ω * (Δtr[end]) - m * (Δφr[end]) + n * 1.0pi
    dθphase = ω .* dtθ .- m .* dφθ .+ k
    θphaseL = ω * Δtθ[1] - m * Δφθ[1]
    θphaseR = ω * Δtθ[end] - m * Δφθ[end] + k * 1.0pi

    # --------------------------
    # 7. W-series Quantities (Vectorized)
    # --------------------------
    Y_grid = reshape(Y, :, 1)
    Yp_grid = reshape(Yp, :, 1)
    X_grid = reshape(X, :, 1)

    # Wll: use precomputed phase_grid
    Wll_grid = ρbar_grid .* L1pL2pS_grid .* r2_grid .* Y_grid .* phase_grid .* 2.0

    # Wlm
    term_A = L2S_grid .+ im .* a .* sinθ_grid .* (ρ_grid .- ρbar_grid) .* S0_grid
    term_B = 2.0 .* Y_grid .+ r_grid .* Yp_grid
    term_C = L2S_grid .* (ρ_grid .+ ρbar_grid) .* r_grid .* Y_grid
    Wlm_grid = 4.0 .* (r_grid .* ρbar_grid .* (term_A .* term_B .+ term_C) ./ ρ_grid) .* phase_grid

    # Wmm
    term_D = X_grid ./ (2.0 .* sqrt.(r2_grid .+ a^2))
    term_E = (Y_grid .+ 2.0 .* r_grid .* Yp_grid) .* phase_grid
    term_F = ρ_grid .* r_grid .* (2.0 .* Y_grid .+ r_grid .* Yp_grid) .* phase_grid
    Wmm_grid = 4.0 .* S0_grid .* ρbar_grid .* (term_D .+ term_E .+ term_F) ./ ρ_grid

    # --------------------------
    # 8. L/M Terms (Vectorized)
    # --------------------------
    urp_grid = reshape(urp, :, 1)
    urm_grid = reshape(urm, :, 1)
    uθp_grid = reshape(uθp, 1, :)
    uθm_grid = reshape(uθm, 1, :)

    # Lp/Lm (using precomputed Δ_grid)
    numerator_L = E .* (r2_grid .+ a^2) .- a .* Lz
    Lp_grid = (numerator_L .- urp_grid) ./ Δ_grid
    Lm_grid = (numerator_L .- urm_grid) ./ Δ_grid

    # Mp/Mm
    term_M = .- im .* sinθ_grid .* (a .* E .- Lz .* inv_sin2θ_grid)
    Mp_grid = term_M .+ uθp_grid
    Mm_grid = term_M .+ uθm_grid

    # --------------------------
    # 9. J-series Terms (Vectorized)
    # --------------------------
    Lp_sq = Lp_grid .^ 2
    Lm_sq = Lm_grid .^ 2
    Mp_sq = Mp_grid .^ 2
    Mm_sq = Mm_grid .^ 2

    Jpp_grid = (Wll_grid .* Lp_sq .+ Wlm_grid .* Lp_grid .* Mp_grid .+ Wmm_grid .* Mp_sq) ./ (4 .* π .^2)
    Jpm_grid = (Wll_grid .* Lp_sq .+ Wlm_grid .* Lp_grid .* Mm_grid .+ Wmm_grid .* Mm_sq) ./ (4 .* π .^2)
    Jmp_grid = (Wll_grid .* Lm_sq .+ Wlm_grid .* Lm_grid .* Mp_grid .+ Wmm_grid .* Mp_sq) ./ (4 .* π .^2)
    Jmm_grid = (Wll_grid .* Lm_sq .+ Wlm_grid .* Lm_grid .* Mm_grid .+ Wmm_grid .* Mm_sq) ./ (4 .* π .^2)

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
    Xi = m * (Δφr_func(qr0) + Δφθ_func(qθ0) - qφ0) - ω * (Δtr_func(qr0) + Δtθ_func(qθ0) - qt0) - k * qθ0 - n * qr0
    prefactor = .- im .* π .* factor .* η .* exp(im * Xi) ./ (Γ .* κ .* Cinc)

    return Jpp_grid, Jpm_grid, Jmp_grid, Jmm_grid, drphase, dθphase, rphaseL, rphaseR, θphaseL, θphaseR, prefactor
end

function integrand_eccentric_sample_m2(KG_samp, Y_samp, SH, n::Int64)
    # --------------------------
    # 1. Extract Parameters (1D Radial, θ as constant)
    # --------------------------
    r_vec = KG_samp["r"]          # Radial sampling points (length N)
    rs_vec = KG_samp["rs"]        # Precomputed tortoise coordinates (length N)
    θ = KG_samp["θ"]              # Fixed at π/2 (scalar)
    N = KG_samp["N_sample"]       # Number of radial grid points
    Γ = KG_samp["Γ"]
    qt0, qr0, _, qφ0 = KG_samp["InitialPhases"]
    Δtr_func, Δφr_func = KG_samp["CrossFunction"]  # Only radial cross functions

    # Radial velocities
    urp = KG_samp["ur_fwd"]       # Forward radial velocity
    urm = KG_samp["ur_rev"]       # Reverse radial velocity

    # Radial cross functions
    Δtr = KG_samp["Δtr"]          # Radial time difference
    Δφr = KG_samp["Δφr"]          # Radial azimuthal difference

    # Orbital constants
    E = KG_samp["E"]              # Energy (scalar)
    Lz = KG_samp["Lz"]            # z-angular momentum (scalar)

    # Spherical harmonics (θ-fixed, retained as scalars)
    S0, S1, S2 = SH

    # Constant parameters
    params = Y_samp["params"]
    s, l, m = params.s, params.l, params.m
    a, ω, λ = params.a, params.omega, params.lambda

    # Black hole horizons
    rp = 1.0 + sqrt(1.0 - a^2)    # Outer horizon
    rm = 1.0 - sqrt(1.0 - a^2)    # Inner horizon

    # --------------------------
    # 2. θ-fixed Terms (retained as non-zero constants)
    # --------------------------

    # --------------------------
    # 3. Complex Parameters (1D radial vectors)
    # --------------------------
    ρ = -1.0 ./ r_vec             # Simplified to -1/r since cosθ=0
    ρbar = -1.0 ./ r_vec          # Equal to ρ due to cosθ=0
    inv_ρ = -r_vec                # 1/ρ = -r

    # --------------------------
    # 4. L-series Operators (θ-fixed, broadcast to 1D)
    # --------------------------
    L1 = (-m) + a*ω               # Simplified from (-m)*1 + aω*1 + 0 (cotθ=0)
    L2 = (-m) + a*ω               # Simplified from (-m)*1 + aω*1 + 0 (2*cotθ=0)
    L2p = -2.0                    # Simplified from 0 + 0 - 2*1 (cosθ=0)

    # --------------------------
    # 5. Spherical Harmonic Combinations (θ-fixed, 1D vectors)
    # --------------------------
    L2S = S1 + L2 * S0
    L1Sp = S2 + L1 * S1
    L1L2S = L1Sp + L2p * S0 + L2 * S1 + L1 * L2 * S0

    # L1pL2pS (retain θ-related constant terms)
    term1 = inv_ρ .* L1L2S
    term2 = 3im * a * L1 * S0     # Simplified since sinθ=1
    term3 = 0.0                   # Result is 0 due to cosθ=0
    term4 = 2im * a * S1          # Simplified since sinθ=1
    term5 = -1im * a * L2 * S0    # Simplified since sinθ=1
    L1pL2pS = term1 .+ term2 .+ term3 .+ term4 .+ term5

    # --------------------------
    # 6. Phase Terms (radial-only)
    # --------------------------
    log_term = log.((r_vec .- rp) ./ (r_vec .- rm))
    phase = exp.(im .* ω .* rs_vec .- im .* a .* m .* log_term ./ (2 .* sqrt.(1 .- a .^2)))

    qr = range(0, π; length=N)
    rphase = ω .* Δtr .- m .* Δφr .+ n .* qr  # Radial-only phase

    # --------------------------
    # 7. W-series Quantities (1D radial)
    # --------------------------
    Y = Y_samp["Y"]
    Yp = Y_samp["Yp"]
    X = Y_samp["X"]
    Binc = Y_samp["Binc"]

    Wnn = ρbar .* L1pL2pS .* (r_vec .^ 2) .* Y .* phase ./ 2.0

    term_A = L2S .+ im * a * (ρ .- ρbar) .* S0  # Simplified since sinθ=1
    term_B = 2Y .+ r_vec .* Yp
    term_C = L2S .* (ρ .+ ρbar) .* r_vec .* Y
    Wnmbar = (-r_vec .* ρbar .* (term_A .* term_B .+ term_C) ./ ρ) .* phase

    term_D = X ./ (2 .* sqrt.(r_vec .^ 2 .+ a^2))
    term_E = (Y .+ 2r_vec .* Yp) .* phase
    term_F = ρ .* r_vec .* (2Y .+ r_vec .* Yp) .* phase
    Wmbarmbar = S0 .* ρbar .* (term_D .+ term_E .+ term_F) ./ ρ

    # --------------------------
    # 8. N and Mbar Terms (retain non-zero Mbar constant)
    # --------------------------
    Δ = r_vec .^ 2 .- 2 .* r_vec .+ a.^2  # Kerr metric Δ(r)
    numerator_N = E .* (r_vec .^ 2 .+ a.^2) .- a .* Lz
    Np = (numerator_N .+ urp) ./ Δ
    Nm = (numerator_N .+ urm) ./ Δ

    # Core Mbar term (uθ=0 for fixed θ, but base term is non-zero)
    term_M = im * (a*E - Lz)       # Simplified since sinθ=1, inv_sin2θ=1
    Mbarp = term_M
    Mbarm = term_M

    # --------------------------
    # 9. J-series Terms (1D radial, corrected normalization)
    # --------------------------
    Np_sq = Np .^ 2
    Nm_sq = Nm .^ 2
    Mbarp_sq = Mbarp .^ 2
    Mbarm_sq = Mbarm .^ 2

    exp_p = exp.(im .* rphase)
    exp_m = exp.(-im .* rphase)

    # 1D Fourier transform: normalization factor 1/(2π) (radial-only)
    Jp = (Wnn .* Np_sq .+ Wnmbar .* Np .* Mbarp .+ Wmbarmbar .* Mbarp_sq) .* exp_p
    Jm = (Wnn .* Nm_sq .+ Wnmbar .* Nm .* Mbarm .+ Wmbarmbar .* Mbarm_sq) .* exp_m
    Jlmn = (Jp .+ Jm) ./ (2π)  # Corrected: 1D Fourier normalization (original 2D: 1/(4π²))

    # --------------------------
    # 10. Final Phase Correction and Result
    # --------------------------
    Xi = m*(Δφr_func(qr0) - qφ0) - ω*(Δtr_func(qr0) - qt0) - n*qr0  # No θ-related terms

    return 4im .* π .* ω .* Jlmn .* exp.(im .* Xi) ./ (Γ .* Binc)
end

function integrand_eccentric_sample_p2(KG_samp, Y_samp, SH, n::Int64)
    # --------------------------
    # 1. Extract Parameters (1D Radial, θ as constant)
    # --------------------------
    r_vec = KG_samp["r"]          # Radial sampling points (length N)
    rs_vec = KG_samp["rs"]        # Precomputed tortoise coordinates (length N)
    θ = KG_samp["θ"]              # Fixed at π/2 (scalar)
    N = KG_samp["N_sample"]       # Number of radial grid points
    Γ = KG_samp["Γ"]
    qt0, qr0, _, qφ0 = KG_samp["InitialPhases"]
    Δtr_func, Δφr_func = KG_samp["CrossFunction"]  # Only radial cross functions

    # Radial velocities
    urp = KG_samp["ur_fwd"]       # Forward radial velocity
    urm = KG_samp["ur_rev"]       # Reverse radial velocity

    # Radial cross functions
    Δtr = KG_samp["Δtr"]          # Radial time difference
    Δφr = KG_samp["Δφr"]          # Radial azimuthal difference

    # Orbital constants
    E = KG_samp["E"]              # Energy (scalar)
    Lz = KG_samp["Lz"]            # z-angular momentum (scalar)

    # Spherical harmonics (θ-fixed, retained as scalars)
    S0, S1, S2 = SH

    # Constant parameters
    params = Y_samp["params"]
    s, l, m = params.s, params.l, params.m
    a, ω, λ = params.a, params.omega, params.lambda

    # Black hole horizons
    rp = 1.0 + sqrt(1.0 - a^2)    # Outer horizon
    rm = 1.0 - sqrt(1.0 - a^2)    # Inner horizon

    # --------------------------
    # 2. θ-fixed Terms (retained as non-zero constants)
    # --------------------------

    # --------------------------
    # 3. Complex Parameters (1D radial vectors)
    # --------------------------
    ρ = -1.0 ./ r_vec             # Simplified to -1/r since cosθ=0
    ρbar = -1.0 ./ r_vec          # Equal to ρ due to cosθ=0
    inv_ρ = -r_vec                # 1/ρ = -r

    # --------------------------
    # 4. L-series Operators (θ-fixed, broadcast to 1D)
    # --------------------------
    L1 = m - a*ω              
    L2 = m - a*ω              
    L2p = -2.0                    # Simplified from 0 + 0 - 2*1 (cosθ=0)

    # --------------------------
    # 5. Spherical Harmonic Combinations (θ-fixed, 1D vectors)
    # --------------------------
    L2S = S1 + L2 * S0
    L1Sp = S2 + L1 * S1
    L1L2S = L1Sp + L2p * S0 + L2 * S1 + L1 * L2 * S0

    # L1pL2pS (retain θ-related constant terms)
    term1 = inv_ρ .* L1L2S
    term2 = 3im * a * L1 * S0     # Simplified since sinθ=1
    term3 = 0.0                   # Result is 0 due to cosθ=0
    term4 = 2im * a * S1          # Simplified since sinθ=1
    term5 = -1im * a * L2 * S0    # Simplified since sinθ=1
    L1pL2pS = term1 .+ term2 .+ term3 .+ term4 .+ term5

    # --------------------------
    # 6. Phase Terms (radial-only)
    # --------------------------
    log_term = log.((r_vec .- rp) ./ (r_vec .- rm))
    phase = exp.(.- im .* ω .* rs_vec .+ im .* a .* m .* log_term ./ (2 .* sqrt.(1 .- a .^2)))

    qr = range(0, π; length=N)
    rphase = ω .* Δtr .- m .* Δφr .+ n .* qr  # Radial-only phase

    # --------------------------
    # 7. W-series Quantities (1D radial)
    # --------------------------
    Y = Y_samp["Y"]
    Yp = Y_samp["Yp"]
    X = Y_samp["X"]
    Cinc = Y_samp["Cinc"]

    Wll = 2.0 .* ρbar .* L1pL2pS .* (r_vec .^ 2) .* Y .* phase

    term_A = L2S .+ im * a * (ρ .- ρbar) .* S0
    term_B = 2Y .+ r_vec .* Yp
    term_C = L2S .* (ρ .+ ρbar) .* r_vec .* Y
    Wlm = 4.0 .* (r_vec .* ρbar .* (term_A .* term_B .+ term_C) ./ ρ) .* phase

    term_D = X ./ (2 .* sqrt.(r_vec .^ 2 .+ a^2))
    term_E = (Y .+ 2r_vec .* Yp) .* phase
    term_F = ρ .* r_vec .* (2Y .+ r_vec .* Yp) .* phase
    Wmm = 4.0 * S0 .* ρbar .* (term_D .+ term_E .+ term_F) ./ ρ

    # --------------------------
    # 8. L and M Terms (retain non-zero M constant)
    # --------------------------
    Δ = r_vec .^ 2 .- 2 .* r_vec .+ a.^2  # Kerr metric Δ(r)
    numerator_L = E .* (r_vec .^ 2 .+ a.^2) .- a .* Lz
    Lp = (numerator_L .- urp) ./ Δ
    Lm = (numerator_L .- urm) ./ Δ

    # Core M term (uθ=0 for fixed θ, but base term is non-zero)
    term_M = - im * (a*E - Lz)  
    Mp = term_M
    Mm = term_M

    # --------------------------
    # 9. J-series Terms (1D radial, corrected normalization)
    # --------------------------
    Lp_sq = Lp .^ 2
    Lm_sq = Lm .^ 2
    Mp_sq = Mp .^ 2
    Mm_sq = Mm .^ 2

    exp_p = exp.(im .* rphase)
    exp_m = exp.(-im .* rphase)

    # 1D Fourier transform: normalization factor 1/(2π) (radial-only)
    Jp = (Wll .* Lp_sq .+ Wlm .* Lp .* Mp .+ Wmm .* Mp_sq) .* exp_p
    Jm = (Wll .* Lm_sq .+ Wlm .* Lm .* Mm .+ Wmm .* Mm_sq) .* exp_m
    Jlmn = (Jp .+ Jm) ./ (2π)  # Corrected: 1D Fourier normalization (original 2D: 1/(4π²))

    # --------------------------
    # 10. Final Phase Correction and Result
    # --------------------------
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
    Xi = m*(Δφr_func(qr0) - qφ0) - ω*(Δtr_func(qr0) - qt0) - n*qr0  # No θ-related terms

    return .- im .* π .* factor .* η .* Jlmn .* exp.(im .* Xi) ./ (Γ .* Cinc .* κ)
end

function integrand_eccentric_sample_cheby_m2(KG_samp, Y_samp, SH, n::Int64)
    # --------------------------
    # 1. Extract Parameters (1D Radial, θ as constant)
    # --------------------------
    r_vec = KG_samp["r"]          # Radial sampling points (length N)
    rs_vec = KG_samp["rs"]        # Precomputed tortoise coordinates (length N)
    θ = KG_samp["θ"]              # Fixed at π/2 (scalar)
    N = KG_samp["N_sample"]       # Number of radial grid points
    Γ = KG_samp["Γ"]
    qt0, qr0, _, qφ0 = KG_samp["InitialPhases"]
    Δtr_func, Δφr_func = KG_samp["CrossFunction"]  # Only radial cross functions
    dtr = KG_samp["dtr"]
    dφr = KG_samp["dφr"]

    # Radial velocities
    urp = KG_samp["ur_fwd"]       # Forward radial velocity
    urm = KG_samp["ur_rev"]       # Reverse radial velocity

    # Radial cross functions
    Δtr = KG_samp["Δtr"]          # Radial time difference
    Δφr = KG_samp["Δφr"]          # Radial azimuthal difference

    # Orbital constants
    E = KG_samp["E"]              # Energy (scalar)
    Lz = KG_samp["Lz"]            # z-angular momentum (scalar)

    # Spherical harmonics (θ-fixed, retained as scalars)
    S0, S1, S2 = SH

    # Constant parameters
    params = Y_samp["params"]
    s, l, m = params.s, params.l, params.m
    a, ω, λ = params.a, params.omega, params.lambda

    # Black hole horizons
    rp = 1.0 + sqrt(1.0 - a^2)    # Outer horizon
    rm = 1.0 - sqrt(1.0 - a^2)    # Inner horizon

    # --------------------------
    # 2. θ-fixed Terms (retained as non-zero constants)
    # --------------------------

    # --------------------------
    # 3. Complex Parameters (1D radial vectors)
    # --------------------------
    ρ = -1.0 ./ r_vec             # Simplified to -1/r since cosθ=0
    ρbar = -1.0 ./ r_vec          # Equal to ρ due to cosθ=0
    inv_ρ = -r_vec                # 1/ρ = -r

    # --------------------------
    # 4. L-series Operators (θ-fixed, broadcast to 1D)
    # --------------------------
    L1 = (-m) + a*ω               # Simplified from (-m)*1 + aω*1 + 0 (cotθ=0)
    L2 = (-m) + a*ω               # Simplified from (-m)*1 + aω*1 + 0 (2*cotθ=0)
    L2p = -2.0                    # Simplified from 0 + 0 - 2*1 (cosθ=0)

    # --------------------------
    # 5. Spherical Harmonic Combinations (θ-fixed, 1D vectors)
    # --------------------------
    L2S = S1 + L2 * S0
    L1Sp = S2 + L1 * S1
    L1L2S = L1Sp + L2p * S0 + L2 * S1 + L1 * L2 * S0

    # L1pL2pS (retain θ-related constant terms)
    term1 = inv_ρ .* L1L2S
    term2 = 3im * a * L1 * S0     # Simplified since sinθ=1
    term3 = 0.0                   # Result is 0 due to cosθ=0
    term4 = 2im * a * S1          # Simplified since sinθ=1
    term5 = -1im * a * L2 * S0    # Simplified since sinθ=1
    L1pL2pS = term1 .+ term2 .+ term3 .+ term4 .+ term5

    # --------------------------
    # 6. Phase Terms (radial-only)
    # --------------------------
    log_term = log.((r_vec .- rp) ./ (r_vec .- rm))
    phase = exp.(im .* ω .* rs_vec .- im .* a .* m .* log_term ./ (2 .* sqrt.(1 .- a .^2)))

    drphase = ω .* (dtr) .- m .* (dφr) .+ n 
    rphaseL = ω * (Δtr[1]) - m * (Δφr[1])
    rphaseR = ω * (Δtr[end]) - m * (Δφr[end]) + n * 1.0pi

    # --------------------------
    # 7. W-series Quantities (1D radial)
    # --------------------------
    Y = Y_samp["Y"]
    Yp = Y_samp["Yp"]
    X = Y_samp["X"]
    Binc = Y_samp["Binc"]

    Wnn = ρbar .* L1pL2pS .* (r_vec .^ 2) .* Y .* phase ./ 2.0

    term_A = L2S .+ im * a * (ρ .- ρbar) .* S0  # Simplified since sinθ=1
    term_B = 2Y .+ r_vec .* Yp
    term_C = L2S .* (ρ .+ ρbar) .* r_vec .* Y
    Wnmbar = (-r_vec .* ρbar .* (term_A .* term_B .+ term_C) ./ ρ) .* phase

    term_D = X ./ (2 .* sqrt.(r_vec .^ 2 .+ a^2))
    term_E = (Y .+ 2r_vec .* Yp) .* phase
    term_F = ρ .* r_vec .* (2Y .+ r_vec .* Yp) .* phase
    Wmbarmbar = S0 .* ρbar .* (term_D .+ term_E .+ term_F) ./ ρ

    # --------------------------
    # 8. N and Mbar Terms (retain non-zero Mbar constant)
    # --------------------------
    Δ = r_vec .^ 2 .- 2 .* r_vec .+ a.^2  # Kerr metric Δ(r)
    numerator_N = E .* (r_vec .^ 2 .+ a.^2) .- a .* Lz
    Np = (numerator_N .+ urp) ./ Δ
    Nm = (numerator_N .+ urm) ./ Δ

    # Core Mbar term (uθ=0 for fixed θ, but base term is non-zero)
    term_M = im * (a*E - Lz)       # Simplified since sinθ=1, inv_sin2θ=1
    Mbarp = term_M
    Mbarm = term_M

    # --------------------------
    # 9. J-series Terms (1D radial, corrected normalization)
    # --------------------------
    Np_sq = Np .^ 2
    Nm_sq = Nm .^ 2
    Mbarp_sq = Mbarp .^ 2
    Mbarm_sq = Mbarm .^ 2

    # 1D Fourier transform: normalization factor 1/(2π) (radial-only)
    Jp = (Wnn .* Np_sq .+ Wnmbar .* Np .* Mbarp .+ Wmbarmbar .* Mbarp_sq) ./ (2.0pi)
    Jm = (Wnn .* Nm_sq .+ Wnmbar .* Nm .* Mbarm .+ Wmbarmbar .* Mbarm_sq) ./ (2.0pi)

    Xi = m*(Δφr_func(qr0) - qφ0) - ω*(Δtr_func(qr0) - qt0) - n*qr0 

    prefactor = 4im * π * ω * exp(im * Xi) / (Γ * Binc)

    return Jp, Jm, drphase, rphaseL, rphaseR, prefactor
end

function integrand_eccentric_sample_cheby_p2(KG_samp, Y_samp, SH, n::Int64)
    # --------------------------
    # 1. Extract Parameters (1D Radial, θ as constant)
    # --------------------------
    r_vec = KG_samp["r"]          # Radial sampling points (length N)
    rs_vec = KG_samp["rs"]        # Precomputed tortoise coordinates (length N)
    θ = KG_samp["θ"]              # Fixed at π/2 (scalar)
    N = KG_samp["N_sample"]       # Number of radial grid points
    Γ = KG_samp["Γ"]
    qt0, qr0, _, qφ0 = KG_samp["InitialPhases"]
    Δtr_func, Δφr_func = KG_samp["CrossFunction"]  # Only radial cross functions
    dtr = KG_samp["dtr"]
    dφr = KG_samp["dφr"]

    # Radial velocities
    urp = KG_samp["ur_fwd"]       # Forward radial velocity
    urm = KG_samp["ur_rev"]       # Reverse radial velocity

    # Radial cross functions
    Δtr = KG_samp["Δtr"]          # Radial time difference
    Δφr = KG_samp["Δφr"]          # Radial azimuthal difference

    # Orbital constants
    E = KG_samp["E"]              # Energy (scalar)
    Lz = KG_samp["Lz"]            # z-angular momentum (scalar)

    # Spherical harmonics (θ-fixed, retained as scalars)
    S0, S1, S2 = SH

    # Constant parameters
    params = Y_samp["params"]
    s, l, m = params.s, params.l, params.m
    a, ω, λ = params.a, params.omega, params.lambda

    # Black hole horizons
    rp = 1.0 + sqrt(1.0 - a^2)    # Outer horizon
    rm = 1.0 - sqrt(1.0 - a^2)    # Inner horizon

    # --------------------------
    # 2. θ-fixed Terms (retained as non-zero constants)
    # --------------------------

    # --------------------------
    # 3. Complex Parameters (1D radial vectors)
    # --------------------------
    ρ = -1.0 ./ r_vec             # Simplified to -1/r since cosθ=0
    ρbar = -1.0 ./ r_vec          # Equal to ρ due to cosθ=0
    inv_ρ = -r_vec                # 1/ρ = -r

    # --------------------------
    # 4. L-series Operators (θ-fixed, broadcast to 1D)
    # --------------------------
    L1 = m - a*ω             
    L2 = m - a*ω            
    L2p = -2.0                    # Simplified from 0 + 0 - 2*1 (cosθ=0)

    # --------------------------
    # 5. Spherical Harmonic Combinations (θ-fixed, 1D vectors)
    # --------------------------
    L2S = S1 + L2 * S0
    L1Sp = S2 + L1 * S1
    L1L2S = L1Sp + L2p * S0 + L2 * S1 + L1 * L2 * S0

    # L1pL2pS (retain θ-related constant terms)
    term1 = inv_ρ .* L1L2S
    term2 = 3im * a * L1 * S0     # Simplified since sinθ=1
    term3 = 0.0                   # Result is 0 due to cosθ=0
    term4 = 2im * a * S1          # Simplified since sinθ=1
    term5 = -1im * a * L2 * S0    # Simplified since sinθ=1
    L1pL2pS = term1 .+ term2 .+ term3 .+ term4 .+ term5

    # --------------------------
    # 6. Phase Terms (radial-only)
    # --------------------------
    log_term = log.((r_vec .- rp) ./ (r_vec .- rm))
    phase = exp.(.- im .* ω .* rs_vec .+ im .* a .* m .* log_term ./ (2 .* sqrt.(1 .- a .^2)))

    drphase = ω .* (dtr) .- m .* (dφr) .+ n 
    rphaseL = ω * (Δtr[1]) - m * (Δφr[1])
    rphaseR = ω * (Δtr[end]) - m * (Δφr[end]) + n * 1.0pi

    # --------------------------
    # 7. W-series Quantities (1D radial)
    # --------------------------
    Y = Y_samp["Y"]
    Yp = Y_samp["Yp"]
    X = Y_samp["X"]
    Cinc = Y_samp["Cinc"]

    Wll = 2.0 .* ρbar .* L1pL2pS .* (r_vec .^ 2) .* Y .* phase

    term_A = L2S .+ im * a * (ρ .- ρbar) .* S0
    term_B = 2Y .+ r_vec .* Yp
    term_C = L2S .* (ρ .+ ρbar) .* r_vec .* Y
    Wlm = 4.0 .* (r_vec .* ρbar .* (term_A .* term_B .+ term_C) ./ ρ) .* phase

    term_D = X ./ (2 .* sqrt.(r_vec .^ 2 .+ a^2))
    term_E = (Y .+ 2r_vec .* Yp) .* phase
    term_F = ρ .* r_vec .* (2Y .+ r_vec .* Yp) .* phase
    Wmm = 4.0 * S0 .* ρbar .* (term_D .+ term_E .+ term_F) ./ ρ

    # --------------------------
    # 8. L and M Terms (retain non-zero M constant)
    # --------------------------
    Δ = r_vec .^ 2 .- 2 .* r_vec .+ a.^2  # Kerr metric Δ(r)
    numerator_L = E .* (r_vec .^ 2 .+ a.^2) .- a .* Lz
    Lp = (numerator_L .- urp) ./ Δ
    Lm = (numerator_L .- urm) ./ Δ

    # Core M term (uθ=0 for fixed θ, but base term is non-zero)
    term_M = - im * (a*E - Lz)  
    Mp = term_M
    Mm = term_M

    # --------------------------
    # 9. J-series Terms (1D radial, corrected normalization)
    # --------------------------
    Lp_sq = Lp .^ 2
    Lm_sq = Lm .^ 2
    Mp_sq = Mp .^ 2
    Mm_sq = Mm .^ 2

    # 1D Fourier transform: normalization factor 1/(2π) (radial-only)
    Jp = (Wll .* Lp_sq .+ Wlm .* Lp .* Mp .+ Wmm .* Mp_sq) ./ (2.0pi)
    Jm = (Wll .* Lm_sq .+ Wlm .* Lm .* Mm .+ Wmm .* Mm_sq) ./ (2.0pi)

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
    Xi = m*(Δφr_func(qr0) - qφ0) - ω*(Δtr_func(qr0) - qt0) - n*qr0 

    prefactor = - im * π * factor * η * exp(im * Xi) / (Γ * Cinc * κ)

    return Jp, Jm, drphase, rphaseL, rphaseR, prefactor
end

function integrand_inclined_sample_m2(KG_samp, Y_dic, SH_samp, k::Int64)
    # --------------------------
    # 1. Extract Parameters (1D Polar, r fixed as constant)
    # --------------------------
    r = KG_samp["r"]              # Fixed radial coordinate (scalar)
    rs = KG_samp["rs"]            # Precomputed tortoise coordinate (scalar)
    θ_vec = KG_samp["θ"]          # Polar sampling points (length K)
    K = KG_samp["K_sample"]       # Number of polar grid points
    Γ = KG_samp["Γ"]
    qt0, _, qθ0, qφ0 = KG_samp["InitialPhases"]  # Omit radial initial phase
    Δtθ_func, Δφθ_func = KG_samp["CrossFunction"]  # Only polar cross functions

    # Polar velocities
    uθp = KG_samp["uθ_fwd"]       # Forward polar velocity
    uθm = KG_samp["uθ_rev"]       # Reverse polar velocity

    # Polar cross functions
    Δtθ = KG_samp["Δtθ"]          # Polar time difference
    Δφθ = KG_samp["Δφθ"]          # Polar azimuthal difference

    # Orbital constants
    E = KG_samp["E"]              # Energy (scalar)
    Lz = KG_samp["Lz"]            # z-angular momentum (scalar)

    # Spherical harmonics (polar-dependent, 1D vectors)
    S0 = SH_samp["S0"]            # Base value (length K)
    S1 = SH_samp["S1"]            # 1st θ-derivative (length K)
    S2 = SH_samp["S2"]            # 2nd θ-derivative (length K)

    # Constant parameters
    params = Y_dic["params"]
    s, l, m = params.s, params.l, params.m
    a, ω, λ = params.a, params.omega, params.lambda

    # Black hole horizons (scalar)
    rp = 1.0 + sqrt(1.0 - a^2)    # Outer horizon
    rm = 1.0 - sqrt(1.0 - a^2)    # Inner horizon

    # --------------------------
    # 2. r-fixed Terms (retained as constants)
    # --------------------------
    r2 = r^2                      # r² (fixed)
    Δ = r2 - 2r + a^2             # Kerr metric Δ(r) (fixed)
    sqrt_r2_a2 = sqrt(r2 + a^2)   # √(r² + a²) (fixed)

    # --------------------------
    # 3. Complex Parameters (polar 1D vectors)
    # --------------------------
    cosθ = cos.(θ_vec)
    ρ = -1.0 ./ (r .- im .* a .* cosθ)  # ρ(θ)
    ρbar = -1.0 ./ (r .+ im .* a .* cosθ)  # ρ̄(θ)
    inv_ρ = 1.0 ./ ρ              # 1/ρ(θ)

    # --------------------------
    # 4. L-series Operators (polar 1D vectors)
    # --------------------------
    sinθ = sin.(θ_vec)
    inv_sinθ = 1.0 ./ sinθ
    inv_sin2θ = 1.0 ./ (sinθ .^ 2)
    cotθ = cosθ ./ sinθ

    L1 = (-m) .* inv_sinθ .+ a .* ω .* sinθ .+ cotθ
    L2 = (-m) .* inv_sinθ .+ a .* ω .* sinθ .+ 2.0 .* cotθ
    L2p = m .* cosθ .* inv_sin2θ .+ a .* ω .* cosθ .- 2.0 .* inv_sin2θ

    # --------------------------
    # 5. Spherical Harmonic Combinations (polar 1D vectors)
    # --------------------------
    L2S = S1 .+ L2 .* S0
    L1Sp = S2 .+ L1 .* S1
    L1L2S = L1Sp .+ L2p .* S0 .+ L2 .* S1 .+ L1 .* L2 .* S0

    # L1pL2pS (polar-dependent complex combinations)
    term1 = inv_ρ .* L1L2S
    term2 = 3im .* a .* sinθ .* L1 .* S0
    term3 = 3im .* a .* cosθ .* S0
    term4 = 2im .* a .* sinθ .* S1
    term5 = -1im .* a .* sinθ .* L2 .* S0
    L1pL2pS = term1 .+ term2 .+ term3 .+ term4 .+ term5

    # --------------------------
    # 6. Phase Terms (polar-only)
    # --------------------------
    log_term = log((r - rp) / (r - rm))  # Scalar (r-fixed)
    phase = exp(im * ω * rs - im * a * m * log_term / (2.0 * sqrt(1.0 - a^2)))  # 1D (K)

    qθ = range(0.0, π; length=K)
    θphase = ω .* Δtθ .- m .* Δφθ .+ k .* qθ  # Polar-only phase (1D K)

    # --------------------------
    # 7. W-series Quantities (polar 1D vectors)
    # --------------------------
    Y = Y_dic["Y"]               # Radial function (scalar, r-fixed)
    Yp = Y_dic["Yp"]             # Radial derivative (scalar)
    X = Y_dic["X"]               # Auxiliary radial function (scalar)
    Binc = Y_dic["Binc"]

    # Wnn: radial terms are scalar (r-fixed)
    Wnn = ρbar .* L1pL2pS .* r2 .* Y .* phase ./ 2.0

    term_A = L2S .+ im .* a .* sinθ .* (ρ .- ρbar) .* S0
    term_B = 2.0 .* Y .+ r .* Yp    # Scalar broadcast to 1D
    term_C = L2S .* (ρ .+ ρbar) .* r .* Y
    Wnmbar = (-r .* ρbar .* (term_A .* term_B .+ term_C) ./ ρ) .* phase

    term_D = X ./ (2.0 .* sqrt_r2_a2)  # Scalar broadcast to 1D
    term_E = (Y .+ 2.0 .* r .* Yp) .* phase
    term_F = ρ .* r .* (2.0 .* Y .+ r .* Yp) .* phase
    Wmbarmbar = S0 .* ρbar .* (term_D .+ term_E .+ term_F) ./ ρ

    # --------------------------
    # 8. N and Mbar Terms (retain N as constant, Mbar polar-dependent)
    # --------------------------
    N = (E * (r2 + a^2) - a * Lz) / Δ  # ur_fwd=0 for fixed r → Np = numerator_N / Δ

    # Polar velocity terms (Mbarp/Mbarm are 1D)
    term_M = im .* sinθ .* (a .* E .- Lz .* inv_sin2θ)
    Mbarp = term_M .+ uθp
    Mbarm = term_M .+ uθm

    # --------------------------
    # 9. J-series Terms (1D polar, 1D Fourier normalization)
    # --------------------------             # Same as Np_sq for fixed r
    Mbarp_sq = Mbarp .^ 2
    Mbarm_sq = Mbarm .^ 2

    exp_p = exp.(im .* θphase)
    exp_m = exp.(-im .* θphase)

    # 1D Fourier transform: normalization factor 1/(2π) (polar-only)
    Jp = (Wnn .* N^2 .+ Wnmbar .* N .* Mbarp .+ Wmbarmbar .* Mbarp_sq) .* exp_p
    Jm = (Wnn .* N^2 .+ Wnmbar .* N .* Mbarm .+ Wmbarmbar .* Mbarm_sq) .* exp_m
    Jlmk = (Jp .+ Jm) ./ (2π)  # Corrected: 1D Fourier normalization

    # --------------------------
    # 10. Final Phase Correction and Result
    # --------------------------
    Xi = m * (Δφθ_func(qθ0) - qφ0) - ω * (Δtθ_func(qθ0) - qt0) - k * qθ0  # No radial terms

    return 4im .* π .* ω .* Jlmk .* exp.(im .* Xi) ./ (Γ .* Binc)
end

function integrand_inclined_sample_p2(KG_samp, Y_dic, SH_samp, k::Int64)
    # --------------------------
    # 1. Extract Parameters (1D Polar, r fixed as constant)
    # --------------------------
    r = KG_samp["r"]              # Fixed radial coordinate (scalar)
    rs = KG_samp["rs"]            # Precomputed tortoise coordinate (scalar)
    θ_vec = KG_samp["θ"]          # Polar sampling points (length K)
    K = KG_samp["K_sample"]       # Number of polar grid points
    Γ = KG_samp["Γ"]
    qt0, _, qθ0, qφ0 = KG_samp["InitialPhases"]  # Omit radial initial phase
    Δtθ_func, Δφθ_func = KG_samp["CrossFunction"]  # Only polar cross functions

    # Polar velocities
    uθp = KG_samp["uθ_fwd"]       # Forward polar velocity
    uθm = KG_samp["uθ_rev"]       # Reverse polar velocity

    # Polar cross functions
    Δtθ = KG_samp["Δtθ"]          # Polar time difference
    Δφθ = KG_samp["Δφθ"]          # Polar azimuthal difference

    # Orbital constants
    E = KG_samp["E"]              # Energy (scalar)
    Lz = KG_samp["Lz"]            # z-angular momentum (scalar)

    # Spherical harmonics (polar-dependent, 1D vectors)
    S0 = SH_samp["S0"]            # Base value (length K)
    S1 = SH_samp["S1"]            # 1st θ-derivative (length K)
    S2 = SH_samp["S2"]            # 2nd θ-derivative (length K)

    # Constant parameters
    params = Y_dic["params"]
    s, l, m = params.s, params.l, params.m
    a, ω, λ = params.a, params.omega, params.lambda

    # Black hole horizons (scalar)
    rp = 1.0 + sqrt(1.0 - a^2)    # Outer horizon
    rm = 1.0 - sqrt(1.0 - a^2)    # Inner horizon

    # --------------------------
    # 2. r-fixed Terms (retained as constants)
    # --------------------------
    r2 = r^2                      # r² (fixed)
    Δ = r2 - 2r + a^2             # Kerr metric Δ(r) (fixed)
    sqrt_r2_a2 = sqrt(r2 + a^2)   # √(r² + a²) (fixed)

    # --------------------------
    # 3. Complex Parameters (polar 1D vectors)
    # --------------------------
    cosθ = cos.(θ_vec)
    ρ = -1.0 ./ (r .- im .* a .* cosθ)  # ρ(θ)
    ρbar = -1.0 ./ (r .+ im .* a .* cosθ)  # ρ̄(θ)
    inv_ρ = 1.0 ./ ρ              # 1/ρ(θ)

    # --------------------------
    # 4. L-series Operators (polar 1D vectors)
    # --------------------------
    sinθ = sin.(θ_vec)
    inv_sinθ = 1.0 ./ sinθ
    inv_sin2θ = 1.0 ./ (sinθ .^ 2)
    cotθ = cosθ ./ sinθ

    L1 = m .* inv_sinθ .- a .* ω .* sinθ .+ cotθ
    L2 = m .* inv_sinθ .- a .* ω .* sinθ .+ 2.0 .* cotθ
    L2p = .- m .* cosθ .* inv_sin2θ .- a .* ω .* cosθ .- 2.0 .* inv_sin2θ

    # --------------------------
    # 5. Spherical Harmonic Combinations (polar 1D vectors)
    # --------------------------
    L2S = S1 .+ L2 .* S0
    L1Sp = S2 .+ L1 .* S1
    L1L2S = L1Sp .+ L2p .* S0 .+ L2 .* S1 .+ L1 .* L2 .* S0

    # L1pL2pS (polar-dependent complex combinations)
    term1 = inv_ρ .* L1L2S
    term2 = 3im .* a .* sinθ .* L1 .* S0
    term3 = 3im .* a .* cosθ .* S0
    term4 = 2im .* a .* sinθ .* S1
    term5 = -1im .* a .* sinθ .* L2 .* S0
    L1pL2pS = term1 .+ term2 .+ term3 .+ term4 .+ term5

    # --------------------------
    # 6. Phase Terms (polar-only)
    # --------------------------
    log_term = log((r - rp) / (r - rm))  # Scalar (r-fixed)
    phase = exp(- im * ω * rs + im * a * m * log_term / (2.0 * sqrt(1.0 - a^2)))  # 1D (K)

    qθ = range(0.0, π; length=K)
    θphase = ω .* Δtθ .- m .* Δφθ .+ k .* qθ  # Polar-only phase (1D K)

    # --------------------------
    # 7. W-series Quantities (polar 1D vectors)
    # --------------------------
    Y = Y_dic["Y"]               # Radial function (scalar, r-fixed)
    Yp = Y_dic["Yp"]             # Radial derivative (scalar)
    X = Y_dic["X"]               # Auxiliary radial function (scalar)
    Cinc = Y_dic["Cinc"]

    # Wll: radial terms are scalar (r-fixed)
    Wll = 2.0 .* ρbar .* L1pL2pS .* r2 .* Y .* phase

    term_A = L2S .+ im .* a .* sinθ .* (ρ .- ρbar) .* S0
    term_B = 2.0 .* Y .+ r .* Yp    # Scalar broadcast to 1D
    term_C = L2S .* (ρ .+ ρbar) .* r .* Y
    Wlm = 4.0 .* (r .* ρbar .* (term_A .* term_B .+ term_C) ./ ρ) .* phase

    term_D = X ./ (2.0 .* sqrt_r2_a2)  # Scalar broadcast to 1D
    term_E = (Y .+ 2.0 .* r .* Yp) .* phase
    term_F = ρ .* r .* (2.0 .* Y .+ r .* Yp) .* phase
    Wmm = 4.0 .* S0 .* ρbar .* (term_D .+ term_E .+ term_F) ./ ρ

    # --------------------------
    # 8. L and M Terms
    # --------------------------
    L = (E * (r2 + a^2) - a * Lz) / Δ

    # Polar velocity terms (Mp/Mm are 1D)
    term_M = - im .* sinθ .* (a .* E .- Lz .* inv_sin2θ)
    Mp = term_M .+ uθp
    Mm = term_M .+ uθm

    # --------------------------
    # 9. J-series Terms (1D polar, 1D Fourier normalization)
    # --------------------------           
    Mp_sq = Mp .^ 2
    Mm_sq = Mm .^ 2

    exp_p = exp.(im .* θphase)
    exp_m = exp.(-im .* θphase)

    # 1D Fourier transform: normalization factor 1/(2π) (polar-only)
    Jp = (Wll .* L^2 .+ Wlm .* L .* Mp .+ Wmm .* Mp_sq) .* exp_p
    Jm = (Wll .* L^2 .+ Wlm .* L .* Mm .+ Wmm .* Mm_sq) .* exp_m
    Jlmk = (Jp .+ Jm) ./ (2π)  # Corrected: 1D Fourier normalization

    # --------------------------
    # 10. Final Phase Correction and Result
    # --------------------------
    Xi = m * (Δφθ_func(qθ0) - qφ0) - ω * (Δtθ_func(qθ0) - qt0) - k * qθ0  # No radial terms
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

    return .- im .* π .* η .* Jlmk .* exp.(im .* Xi) .* factor ./ (Γ .* Cinc .* κ)
end

function integrand_inclined_sample_cheby_m2(KG_samp, Y_dic, SH_samp, k::Int64)
    # --------------------------
    # 1. Extract Parameters (1D Polar, r fixed as constant)
    # --------------------------
    r = KG_samp["r"]              # Fixed radial coordinate (scalar)
    rs_vals = KG_samp["rs"]       # Precomputed tortoise coordinates (length K)
    θ_vec = KG_samp["θ"]          # Polar sampling points (length K)
    K = KG_samp["K_sample"]       # Number of polar grid points
    Γ = KG_samp["Γ"]
    qt0, _, qθ0, qφ0 = KG_samp["InitialPhases"]  # Omit radial initial phase
    Δtθ_func, Δφθ_func = KG_samp["CrossFunction"]  # Only polar cross functions
    dtθ = KG_samp["dtθ"]
    dφθ = KG_samp["dφθ"]

    # Polar velocities
    uθp = KG_samp["uθ_fwd"]       # Forward polar velocity
    uθm = KG_samp["uθ_rev"]       # Reverse polar velocity

    # Polar cross functions
    Δtθ = KG_samp["Δtθ"]          # Polar time difference
    Δφθ = KG_samp["Δφθ"]          # Polar azimuthal difference

    # Orbital constants
    E = KG_samp["E"]              # Energy (scalar)
    Lz = KG_samp["Lz"]            # z-angular momentum (scalar)

    # Spherical harmonics (polar-dependent, 1D vectors)
    S0 = SH_samp["S0"]            # Base value (length K)
    S1 = SH_samp["S1"]            # 1st θ-derivative (length K)
    S2 = SH_samp["S2"]            # 2nd θ-derivative (length K)

    # Constant parameters
    params = Y_dic["params"]
    s, l, m = params.s, params.l, params.m
    a, ω, λ = params.a, params.omega, params.lambda

    # Black hole horizons (scalar)
    rp = 1.0 + sqrt(1.0 - a^2)    # Outer horizon
    rm = 1.0 - sqrt(1.0 - a^2)    # Inner horizon

    # --------------------------
    # 2. r-fixed Terms (retained as constants)
    # --------------------------
    r2 = r^2                      # r² (fixed)
    Δ = r2 - 2r + a^2             # Kerr metric Δ(r) (fixed)
    sqrt_r2_a2 = sqrt(r2 + a^2)   # √(r² + a²) (fixed)

    # --------------------------
    # 3. Complex Parameters (polar 1D vectors)
    # --------------------------
    cosθ = cos.(θ_vec)
    ρ = -1.0 ./ (r .- im .* a .* cosθ)  # ρ(θ)
    ρbar = -1.0 ./ (r .+ im .* a .* cosθ)  # ρ̄(θ)
    inv_ρ = 1.0 ./ ρ              # 1/ρ(θ)

    # --------------------------
    # 4. L-series Operators (polar 1D vectors)
    # --------------------------
    sinθ = sin.(θ_vec)
    inv_sinθ = 1.0 ./ sinθ
    inv_sin2θ = 1.0 ./ (sinθ .^ 2)
    cotθ = cosθ ./ sinθ

    L1 = (-m) .* inv_sinθ .+ a .* ω .* sinθ .+ cotθ
    L2 = (-m) .* inv_sinθ .+ a .* ω .* sinθ .+ 2.0 .* cotθ
    L2p = m .* cosθ .* inv_sin2θ .+ a .* ω .* cosθ .- 2.0 .* inv_sin2θ

    # --------------------------
    # 5. Spherical Harmonic Combinations (polar 1D vectors)
    # --------------------------
    L2S = S1 .+ L2 .* S0
    L1Sp = S2 .+ L1 .* S1
    L1L2S = L1Sp .+ L2p .* S0 .+ L2 .* S1 .+ L1 .* L2 .* S0

    # L1pL2pS (polar-dependent complex combinations)
    term1 = inv_ρ .* L1L2S
    term2 = 3im .* a .* sinθ .* L1 .* S0
    term3 = 3im .* a .* cosθ .* S0
    term4 = 2im .* a .* sinθ .* S1
    term5 = -1im .* a .* sinθ .* L2 .* S0
    L1pL2pS = term1 .+ term2 .+ term3 .+ term4 .+ term5

    # --------------------------
    # 6. Phase Terms (polar-only)
    # --------------------------
    log_term = log((r - rp) / (r - rm))  # Scalar (r-fixed)
    phase = exp.(im .* ω .* rs_vals .- im .* a .* m .* log_term ./ (2.0 .* sqrt(1.0 .- a^2)))  # 1D (K)

    qθ = range(0.0, π; length=K)
    dθphase = ω .* dtθ .- m .* dφθ .+ k
    θphaseL = ω * Δtθ[1] - m * Δφθ[1]
    θphaseR = ω * Δtθ[end] - m * Δφθ[end] + k * 1.0pi

    # --------------------------
    # 7. W-series Quantities (polar 1D vectors)
    # --------------------------
    Y = Y_dic["Y"]               # Radial function (scalar, r-fixed)
    Yp = Y_dic["Yp"]             # Radial derivative (scalar)
    X = Y_dic["X"]               # Auxiliary radial function (scalar)
    Binc = Y_dic["Binc"]

    # Wnn: radial terms are scalar (r-fixed)
    Wnn = ρbar .* L1pL2pS .* r2 .* Y .* phase ./ 2.0

    term_A = L2S .+ im .* a .* sinθ .* (ρ .- ρbar) .* S0
    term_B = 2.0 .* Y .+ r .* Yp    # Scalar broadcast to 1D
    term_C = L2S .* (ρ .+ ρbar) .* r .* Y
    Wnmbar = (-r .* ρbar .* (term_A .* term_B .+ term_C) ./ ρ) .* phase

    term_D = X ./ (2.0 .* sqrt_r2_a2)  # Scalar broadcast to 1D
    term_E = (Y .+ 2.0 .* r .* Yp) .* phase
    term_F = ρ .* r .* (2.0 .* Y .+ r .* Yp) .* phase
    Wmbarmbar = S0 .* ρbar .* (term_D .+ term_E .+ term_F) ./ ρ

    # --------------------------
    # 8. N and Mbar Terms (retain N as constant, Mbar polar-dependent)
    # --------------------------
    # Radial velocity terms (Np/Nm are scalar for fixed r)
    numerator_N = E .* (r2 + a^2) .- a .* Lz  # Scalar
    Np = (numerator_N) ./ Δ  # ur_fwd=0 for fixed r → Np = numerator_N / Δ
    Nm = (numerator_N) ./ Δ  # ur_rev=0 for fixed r → Nm = Np

    # Polar velocity terms (Mbarp/Mbarm are 1D)
    term_M = im .* sinθ .* (a .* E .- Lz .* inv_sin2θ)
    Mbarp = term_M .+ uθp
    Mbarm = term_M .+ uθm

    # --------------------------
    # 9. J-series Terms (1D polar, 1D Fourier normalization)
    # --------------------------
    Np_sq = Np^2                  # Scalar squared
    Nm_sq = Nm^2                  # Same as Np_sq for fixed r
    Mbarp_sq = Mbarp .^ 2
    Mbarm_sq = Mbarm .^ 2

    # 1D Fourier transform: normalization factor 1/(2π) (polar-only)
    Jp = (Wnn .* Np_sq .+ Wnmbar .* Np .* Mbarp .+ Wmbarmbar .* Mbarp_sq) ./ (2.0pi)
    Jm = (Wnn .* Nm_sq .+ Wnmbar .* Nm .* Mbarm .+ Wmbarmbar .* Mbarm_sq) ./ (2.0pi)

    # --------------------------
    # 10. Final Phase Correction and Result
    # --------------------------
    Xi = m * (Δφθ_func(qθ0) - qφ0) - ω * (Δtθ_func(qθ0) - qt0) - k * qθ0  # No radial terms
    prefactor = 4im * π * ω * exp(im * Xi) / (Γ * Binc)

    return Jp, Jm, dθphase, θphaseL, θphaseR, prefactor
end

function integrand_inclined_sample_cheby_p2(KG_samp, Y_dic, SH_samp, k::Int64)
    # --------------------------
    # 1. Extract Parameters (1D Polar, r fixed as constant)
    # --------------------------
    r = KG_samp["r"]              # Fixed radial coordinate (scalar)
    rs_vals = KG_samp["rs"]       # Precomputed tortoise coordinates (length K)
    θ_vec = KG_samp["θ"]          # Polar sampling points (length K)
    K = KG_samp["K_sample"]       # Number of polar grid points
    Γ = KG_samp["Γ"]
    qt0, _, qθ0, qφ0 = KG_samp["InitialPhases"]  # Omit radial initial phase
    Δtθ_func, Δφθ_func = KG_samp["CrossFunction"]  # Only polar cross functions
    dtθ = KG_samp["dtθ"]
    dφθ = KG_samp["dφθ"]

    # Polar velocities
    uθp = KG_samp["uθ_fwd"]       # Forward polar velocity
    uθm = KG_samp["uθ_rev"]       # Reverse polar velocity

    # Polar cross functions
    Δtθ = KG_samp["Δtθ"]          # Polar time difference
    Δφθ = KG_samp["Δφθ"]          # Polar azimuthal difference

    # Orbital constants
    E = KG_samp["E"]              # Energy (scalar)
    Lz = KG_samp["Lz"]            # z-angular momentum (scalar)

    # Spherical harmonics (polar-dependent, 1D vectors)
    S0 = SH_samp["S0"]            # Base value (length K)
    S1 = SH_samp["S1"]            # 1st θ-derivative (length K)
    S2 = SH_samp["S2"]            # 2nd θ-derivative (length K)

    # Constant parameters
    params = Y_dic["params"]
    s, l, m = params.s, params.l, params.m
    a, ω, λ = params.a, params.omega, params.lambda

    # Black hole horizons (scalar)
    rp = 1.0 + sqrt(1.0 - a^2)    # Outer horizon
    rm = 1.0 - sqrt(1.0 - a^2)    # Inner horizon

    # --------------------------
    # 2. r-fixed Terms (retained as constants)
    # --------------------------
    r2 = r^2                      # r² (fixed)
    Δ = r2 - 2*r + a^2             # Kerr metric Δ(r) (fixed)
    sqrt_r2_a2 = sqrt(r2 + a^2)   # √(r² + a²) (fixed)

    # --------------------------
    # 3. Complex Parameters (polar 1D vectors)
    # --------------------------
    cosθ = cos.(θ_vec)
    ρ = -1.0 ./ (r .- im .* a .* cosθ)  # ρ(θ)
    ρbar = -1.0 ./ (r .+ im .* a .* cosθ)  # ρ̄(θ)
    inv_ρ = 1.0 ./ ρ              # 1/ρ(θ)

    # --------------------------
    # 4. L-series Operators (polar 1D vectors)
    # --------------------------
    sinθ = sin.(θ_vec)
    inv_sinθ = 1.0 ./ sinθ
    inv_sin2θ = 1.0 ./ (sinθ .^ 2)
    cotθ = cosθ ./ sinθ

    L1 = m .* inv_sinθ .- a .* ω .* sinθ .+ cotθ
    L2 = m .* inv_sinθ .- a .* ω .* sinθ .+ 2.0 .* cotθ
    L2p = - m .* cosθ .* inv_sin2θ .- a .* ω .* cosθ .- 2.0 .* inv_sin2θ

    # --------------------------
    # 5. Spherical Harmonic Combinations (polar 1D vectors)
    # --------------------------
    L2S = S1 .+ L2 .* S0
    L1Sp = S2 .+ L1 .* S1
    L1L2S = L1Sp .+ L2p .* S0 .+ L2 .* S1 .+ L1 .* L2 .* S0

    # L1pL2pS (polar-dependent complex combinations)
    term1 = inv_ρ .* L1L2S
    term2 = 3im .* a .* sinθ .* L1 .* S0
    term3 = 3im .* a .* cosθ .* S0
    term4 = 2im .* a .* sinθ .* S1
    term5 = -1im .* a .* sinθ .* L2 .* S0
    L1pL2pS = term1 .+ term2 .+ term3 .+ term4 .+ term5

    # --------------------------
    # 6. Phase Terms (polar-only)
    # --------------------------
    log_term = log((r - rp) / (r - rm))  # Scalar (r-fixed)
    phase = exp.(.- im .* ω .* rs_vals .+ im .* a .* m .* log_term ./ (2.0 .* sqrt(1.0 .- a^2)))  # 1D (K)

    qθ = range(0.0, π; length=K)
    dθphase = ω .* dtθ .- m .* dφθ .+ k
    θphaseL = ω * Δtθ[1] - m * Δφθ[1]
    θphaseR = ω * Δtθ[end] - m * Δφθ[end] + k * 1.0pi

    # --------------------------
    # 7. W-series Quantities (polar 1D vectors)
    # --------------------------
    Y = Y_dic["Y"]               # Radial function (scalar, r-fixed)
    Yp = Y_dic["Yp"]             # Radial derivative (scalar)
    X = Y_dic["X"]               # Auxiliary radial function (scalar)
    Cinc = Y_dic["Cinc"]

    # Wll: radial terms are scalar (r-fixed)
    Wll = 2.0 .* ρbar .* L1pL2pS .* r2 .* Y .* phase

    term_A = L2S .+ im .* a .* sinθ .* (ρ .- ρbar) .* S0
    term_B = 2.0 .* Y .+ r .* Yp    # Scalar broadcast to 1D
    term_C = L2S .* (ρ .+ ρbar) .* r .* Y
    Wlm = 4.0 .* (r .* ρbar .* (term_A .* term_B .+ term_C) ./ ρ) .* phase

    term_D = X ./ (2.0 .* sqrt_r2_a2)  # Scalar broadcast to 1D
    term_E = (Y .+ 2.0 .* r .* Yp) .* phase
    term_F = ρ .* r .* (2.0 .* Y .+ r .* Yp) .* phase
    Wmm = 4.0 .* S0 .* ρbar .* (term_D .+ term_E .+ term_F) ./ ρ

    # --------------------------
    # 8. L and M Terms
    # --------------------------
    # Radial velocity terms
    L = (E .* (r2 .+ a^2) .- a .* Lz) ./ Δ 

    # Polar velocity terms (Mp/Mm are 1D)
    term_M = .- im .* sinθ .* (a .* E .- Lz .* inv_sin2θ)
    Mp = term_M .+ uθp
    Mm = term_M .+ uθm

    # --------------------------
    # 9. J-series Terms (1D polar, 1D Fourier normalization)
    # --------------------------             
    Mp_sq = Mp .^ 2
    Mm_sq = Mm .^ 2

    # 1D Fourier transform: normalization factor 1/(2π) (polar-only)
    Jp = (Wll .* L .^ 2 .+ Wlm .* L .* Mp .+ Wmm .* Mp_sq) ./ (2.0pi)
    Jm = (Wll .* L .^ 2 .+ Wlm .* L .* Mm .+ Wmm .* Mm_sq) ./ (2.0pi)

    # --------------------------
    # 10. Final Phase Correction and Result
    # --------------------------
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
    Xi = m * (Δφθ_func(qθ0) - qφ0) - ω * (Δtθ_func(qθ0) - qt0) - k * qθ0 
    prefactor = .- im .* π .* η .* exp(im * Xi) .* factor ./ (Γ .* Cinc .* κ)

    return Jp, Jm, dθphase, θphaseL, θphaseR, prefactor
end

function carter_ingredients_sample(KG_samp, a, m, omega)
    θ_vec = KG_samp["θ"]
    E = KG_samp["E"]
    Lz = KG_samp["Lz"]

    Lmnk = (m .* (cot.(θ_vec)) .^ 2 .* Lz .- a .^ 2 .* omega .* (cos.(θ_vec)) .^ 2 .* E) ./ 1.0pi

    return Lmnk
end

function kerr_geo_generic_sample_dense(KG::Dict, N_interval::Int64, K_interval::Int64)
    ispow2(N_interval) || throw(ArgumentError("N_interval must be a power of 2"))
    ispow2(K_interval) || throw(ArgumentError("K_interval must be a power of 2"))
    N_points = N_interval + 1
    K_points = K_interval + 1

    a = KG["a"]
    p = KG["p"]
    _, r, θ, _ = KG["Trajectory"]
    _, ur, uθ, _ = KG["FourVelocity"]
    Δtr, Δtθ, Δφr, Δφθ = KG["CrossFunction"]
    qt0, qr0, qθ0, qφ0 = KG["InitialPhases"]
    Frequencies = KG["Frequencies"]
    ϒr = Frequencies["ϒr"]
    ϒθ = Frequencies["ϒθ"]
    Γ = Frequencies["ϒt"]

    rq(qr) = r((qr - qr0) / ϒr)
    θq(qθ) = θ((qθ - qθ0) / ϒθ)
    urq(qr) = (r((qr - qr0) / ϒr)^2 + a^2 * cos(θ((qr - qr0) / ϒr))^2) * ur((qr - qr0) / ϒr)
    uθq(qθ) = (r((qθ - qθ0) / ϒθ)^2 + a^2 * cos(θ((qθ - qθ0) / ϒθ))^2) * uθ((qθ - qθ0) / ϒθ)

    qr_full = collect(range(0, π; length = N_points))
    qθ_full = collect(range(0, π; length = K_points))

    rq_vals = [rq(q) for q in qr_full]
    rs_vals = rstar_from_r.(a, rq_vals)
    θ_vals = [θq(q) for q in qθ_full]
    Δtr_vals = [Δtr(q) for q in qr_full]
    Δφr_vals = [Δφr(q) for q in qr_full]
    Δtθ_vals = [Δtθ(q) for q in qθ_full]
    Δφθ_vals = [Δφθ(q) for q in qθ_full]
    ur_fwd = [urq(q) for q in qr_full]
    uθ_fwd = [uθq(q) for q in qθ_full]

    return Dict(
        "a" => a,
        "p" => p,
        "r" => rq_vals,
        "rs" => rs_vals,
        "θ" => θ_vals,
        "Δtr" => Δtr_vals,
        "Δtθ" => Δtθ_vals,
        "Δφr" => Δφr_vals,
        "Δφθ" => Δφθ_vals,
        "ur_fwd" => ur_fwd,
        "ur_rev" => -ur_fwd,
        "uθ_fwd" => uθ_fwd,
        "uθ_rev" => -uθ_fwd,
        "N_interval" => N_interval,
        "K_interval" => K_interval,
        "N_sample" => N_points,
        "K_sample" => K_points,
        "E" => KG["Energy"],
        "Lz" => KG["AngularMomentum"],
        "Γ" => Γ,
        "Frequencies" => Frequencies,
        "Trajectory" => KG["Trajectory"],
        "FourVelocity" => KG["FourVelocity"],
        "CrossFunction" => KG["CrossFunction"],
        "InitialPhases" => KG["InitialPhases"],
        "initialPhases" => KG["InitialPhases"],
    )
end

function subsample_generic_sample(KG_sample::Dict, N_interval::Int64, K_interval::Int64)
    ispow2(N_interval) || throw(ArgumentError("N_interval must be a power of 2"))
    ispow2(K_interval) || throw(ArgumentError("K_interval must be a power of 2"))
    full_N = length(KG_sample["r"])
    full_K = length(KG_sample["θ"])
    target_N = N_interval + 1
    target_K = K_interval + 1
    full_N_minus = full_N - 1
    full_K_minus = full_K - 1
    target_N_minus = target_N - 1
    target_K_minus = target_K - 1
    full_N_minus % target_N_minus == 0 || throw(ArgumentError("N_interval does not divide the stored generic grid"))
    full_K_minus % target_K_minus == 0 || throw(ArgumentError("K_interval does not divide the stored generic grid"))
    stepN = div(full_N_minus, target_N_minus)
    stepK = div(full_K_minus, target_K_minus)
    idxN = 1:stepN:full_N
    idxK = 1:stepK:full_K
    sample = Dict{Any, Any}()
    for key in ("r", "rs", "Δtr", "Δφr", "ur_fwd", "ur_rev")
        sample[key] = KG_sample[key][idxN]
    end
    for key in ("θ", "Δtθ", "Δφθ", "uθ_fwd", "uθ_rev")
        sample[key] = KG_sample[key][idxK]
    end
    for key in ("a", "p", "E", "Lz", "Γ", "Frequencies", "Trajectory", "FourVelocity", "CrossFunction", "InitialPhases", "initialPhases")
        sample[key] = KG_sample[key]
    end
    sample["N_interval"] = N_interval
    sample["K_interval"] = K_interval
    sample["N_sample"] = target_N
    sample["K_sample"] = target_K
    return sample
end

function refine_generic_y_sample(Y_soln, Y_sample::Dict, next_sample::Dict)
    return _refine_y_sample(Y_soln, Y_sample, next_sample)
end

function refine_generic_y_sample_threaded(Y_soln, Y_sample::Dict, next_sample::Dict)
    return _refine_y_sample(Y_soln, Y_sample, next_sample; threaded = true)
end

refine_generic_swsh_sample(SH, prev_SHsamp::Dict, next_sample::Dict) = refine_inclined_swsh_sample(SH, prev_SHsamp, next_sample)

function refine_generic_swsh_sample_threaded(SH, prev_SHsamp::Dict, next_sample::Dict)
    prev_points = length(prev_SHsamp["S0"])
    next_points = length(next_sample["θ"])
    next_points >= prev_points || throw(ArgumentError("next_sample must not be coarser than prev_SHsamp"))
    next_points == prev_points && return prev_SHsamp
    refined = Dict{String, Vector{ComplexF64}}(
        "S0" => Vector{ComplexF64}(undef, next_points),
        "S1" => Vector{ComplexF64}(undef, next_points),
        "S2" => Vector{ComplexF64}(undef, next_points),
    )
    s = SH.params.s
    m = SH.params.m
    c = SH.params.c
    λ = SH.lambda
    Threads.@threads :static for i in 1:next_points
        if isodd(i)
            j = (i + 1) >>> 1
            refined["S0"][i] = prev_SHsamp["S0"][j]
            refined["S1"][i] = prev_SHsamp["S1"][j]
            refined["S2"][i] = prev_SHsamp["S2"][j]
        else
            s0, s1, s2 = _swsh_values(SH, next_sample["θ"][i], s, m, c, λ)
            refined["S0"][i] = s0
            refined["S1"][i] = s1
            refined["S2"][i] = s2
        end
    end
    return refined
end

end
