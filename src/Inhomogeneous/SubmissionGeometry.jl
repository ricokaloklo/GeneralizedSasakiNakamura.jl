# Included inside ConvolutionIntegrals. These owners contain one submission's
# orbit/grid data only. Every mode cache retains its own phase factors and scratch.
# Arrays below are read-only by contract while the submission is active. Returning
# them through the scientific API does not transfer permission to deep-clear them.

abstract type AbstractSubmissionGeometryOwner end

struct GenericM2Geometry
    r::Vector{Float64}
    rs::Vector{Float64}
    theta::Vector{Float64}
    dtr::Vector{Float64}
    dttheta::Vector{Float64}
    dphir::Vector{Float64}
    dphitheta::Vector{Float64}
    cross::Vector{Function}
    initial_phases::NTuple{4, Float64}
    wr::Vector{Float64}
    wtheta::Vector{Float64}
    qr::Vector{Float64}
    qtheta::Vector{Float64}
    Np::Vector{Float64}
    Nm::Vector{Float64}
    Np2::Vector{Float64}
    Nm2::Vector{Float64}
    Lp::Vector{Float64}
    Lm::Vector{Float64}
    Lp2::Vector{Float64}
    Lm2::Vector{Float64}
    st::Vector{Float64}
    ct::Vector{Float64}
    invst::Vector{Float64}
    invst2::Vector{Float64}
    termM::Vector{ComplexF64}
    uthetap::Vector{Float64}
    uthetam::Vector{Float64}
    rho::Matrix{ComplexF64}
    rhobar::Matrix{ComplexF64}
    invrho::Matrix{ComplexF64}
    rho_minus::Matrix{ComplexF64}
    rho_plus::Matrix{ComplexF64}
    N::Int
    K::Int
    Gamma::Float64
    E::Float64
    Lz::Float64
    a::Float64
end

mutable struct GenericGeometryOwner <: AbstractSubmissionGeometryOwner
    master::Union{Nothing, Dict}
    trap_samples::Dict{Tuple{Int, Int}, Dict}
    cheby_samples::Dict{Tuple{Int, Int}, Dict}
    adaptive_segments::Dict{Tuple, Dict}
    adaptive_weights::Dict{Tuple, Vector{Float64}}
    geometries::IdDict{Dict, GenericM2Geometry}
    lock::ReentrantLock
    context::NamedTuple
    stats::Dict{String, Int}
    observer::Union{Nothing, Function}
    closed::Bool
end

mutable struct InclinedGeometryOwner <: AbstractSubmissionGeometryOwner
    master::Union{Nothing, Dict}
    samples::Dict{Int, Dict}
    cheby_samples::Dict{Int, Dict}
    levin_trap_samples::Dict{Int, Dict}
    lock::ReentrantLock
    context::NamedTuple
    stats::Dict{String, Int}
    observer::Union{Nothing, Function}
    closed::Bool
end

mutable struct EccentricGeometryOwner <: AbstractSubmissionGeometryOwner
    master::Union{Nothing, Dict}
    samples::Dict{Int, Dict}
    cheby_samples::Dict{Int, Dict}
    adaptive_segments::Dict{Tuple{Int, Int, Int}, Any}
    lock::ReentrantLock
    context::NamedTuple
    stats::Dict{String, Int}
    observer::Union{Nothing, Function}
    closed::Bool
end

# Numerical context captured when an owner is built; the owner refuses use once it changes.
function _submission_geometry_context()
    return (precision_bits = precision(BigFloat), rounding_mode = rounding(BigFloat),
        kg_module = KerrGeodesics, kg_version = Base.pkgversion(KerrGeodesics),
        kg_path = pathof(KerrGeodesics))
end

GenericGeometryOwner(master::Dict; observer = nothing) = GenericGeometryOwner(master,
    Dict{Tuple{Int, Int}, Dict}(), Dict{Tuple{Int, Int}, Dict}(), Dict{Tuple, Dict}(),
    Dict{Tuple, Vector{Float64}}(), IdDict{Dict, GenericM2Geometry}(), ReentrantLock(),
    _submission_geometry_context(), Dict{String, Int}(), observer, false)
InclinedGeometryOwner(master::Dict; observer = nothing) = InclinedGeometryOwner(master, Dict{Int, Dict}(),
    Dict{Int, Dict}(), Dict{Int, Dict}(), ReentrantLock(),
    _submission_geometry_context(), Dict{String, Int}(), observer, false)
EccentricGeometryOwner(master::Dict; observer = nothing) = EccentricGeometryOwner(master, Dict{Int, Dict}(),
    Dict{Int, Dict}(), Dict{Tuple{Int, Int, Int}, Any}(), ReentrantLock(),
    _submission_geometry_context(), Dict{String, Int}(), observer, false)

function _check_submission_geometry_unlocked(owner::AbstractSubmissionGeometryOwner, master)
    owner.closed && throw(ArgumentError("Submission geometry owner is closed"))
    owner.master === master || throw(ArgumentError("Mode cache belongs to another orbit master"))
    owner.context == _submission_geometry_context() ||
        throw(ArgumentError("Submission precision/rounding/dependency context changed"))
    return nothing
end

function check_submission_geometry(owner::AbstractSubmissionGeometryOwner, master)
    return lock(owner.lock) do
        _check_submission_geometry_unlocked(owner, master)
    end
end

function _submission_geometry_increment!(owner, key::String)
    owner.stats[key] = get(owner.stats, key, 0) + 1
    return nothing
end

function _submission_geometry_lookup!(build::Function, owner::AbstractSubmissionGeometryOwner,
        field::Symbol, key, master)
    # The lookup, miss, construction and publication share one lock. Failed builds
    # publish no value. No unlocked get!/haskey/read races with another builder.
    # Optional diagnostic observer is fixed at construction, before tasks start.
    # It may notify/gate, but must not retain science keys or reenter owner APIs.
    observer = owner.observer
    observer === nothing || Base.invokelatest(observer, :lookup_requested, field, key)
    return lock(owner.lock) do
        _check_submission_geometry_unlocked(owner, master)
        memo = getfield(owner, field)
        prefix = string(field)
        _submission_geometry_increment!(owner, prefix * "_lookups")
        if haskey(memo, key)
            _submission_geometry_increment!(owner, prefix * "_hits")
            return memo[key]
        end
        _submission_geometry_increment!(owner, prefix * "_misses")
        observer === nothing || Base.invokelatest(observer, :build_started, field, key)
        value = build()
        memo[key] = value
        _submission_geometry_increment!(owner, prefix * "_builds")
        observer === nothing || Base.invokelatest(observer, :published, field, key)
        return value
    end
end

function submission_geometry_stats(owner::AbstractSubmissionGeometryOwner)
    return lock(owner.lock) do
        copy(owner.stats)
    end
end
generic_geometry_stats(owner::GenericGeometryOwner) = submission_geometry_stats(owner)
inclined_geometry_stats(owner::InclinedGeometryOwner) = submission_geometry_stats(owner)
eccentric_geometry_stats(owner::EccentricGeometryOwner) = submission_geometry_stats(owner)

function _release_submission_geometry!(owner::AbstractSubmissionGeometryOwner)
    # Admission must already be stopped and submitted tasks joined by the caller.
    # Only unlink this owner's containers; scientific API payloads remain valid.
    lock(owner.lock) do
        owner.closed = true
        for field in fieldnames(typeof(owner))
            field in (:stats, :master) && continue
            value = getfield(owner, field)
            value isa AbstractDict && empty!(value)
        end
        owner.master = nothing
        owner.observer = nothing
    end
    return nothing
end
release_generic_geometry!(owner::GenericGeometryOwner) = _release_submission_geometry!(owner)
release_inclined_geometry!(owner::InclinedGeometryOwner) = _release_submission_geometry!(owner)
release_eccentric_geometry!(owner::EccentricGeometryOwner) = _release_submission_geometry!(owner)

function _submission_public_trajectory(a, p, e, x, trajectory)
    trajectory === nothing && return kerr_geo_orbit(a, p, e, x)
    trajectory isa Dict || throw(ArgumentError("trajectory must be a KG orbit Dict"))
    for (key, expected) in (("a", a), ("p", p), ("e", e), ("Cosθ_inc", x))
        haskey(trajectory, key) && isequal(trajectory[key], expected) ||
            throw(ArgumentError("trajectory orbit parameter mismatch: " * key))
    end
    return trajectory
end

function _submission_public_cache(cache, geometry_owner, family::Symbol)
    expected = family === :generic ? GenericM2FluxCache :
        family === :inclined ? InclinedFluxCache : family === :eccentric ? EccentricFluxCache : Nothing
    cache === nothing || cache isa expected || throw(ArgumentError("Mode cache orbit-family mismatch"))
    if geometry_owner !== nothing
        cache !== nothing && cache.geometry_owner === geometry_owner ||
            throw(ArgumentError("geometry_owner requires its own task-private mode cache"))
    end
    if cache !== nothing && cache.geometry_owner !== nothing
        check_submission_geometry(cache.geometry_owner, cache.geometry_owner.master)
    end
    return cache
end

function _submission_public_master(cache, trajectory)
    cache === nothing && return nothing
    owner = cache.geometry_owner
    owner === nothing && return nothing
    check_submission_geometry(owner, owner.master)
    if trajectory !== nothing
        owner.master["Trajectory"] === trajectory["Trajectory"] &&
            owner.master["Frequencies"] === trajectory["Frequencies"] ||
            throw(ArgumentError("Mode cache and trajectory have different orbit owners"))
    end
    return owner.master
end

function _generic_submission_trap_sample(owner, master::Dict, N::Int, K::Int)
    owner === nothing && return GridSampling.subsample_generic_sample(master, N, K)
    owner isa GenericGeometryOwner || throw(ArgumentError("Expected generic geometry owner"))
    return _submission_geometry_lookup!(owner, :trap_samples, (N, K), master) do
        GridSampling.subsample_generic_sample(master, N, K)
    end
end

function _inclined_levin_owner_sample(cache, KG::Dict, K::Int, cheby::Bool)
    owner = cache === nothing ? nothing : cache.geometry_owner
    if owner === nothing
        return cheby ? kerr_geo_inclined_sample_cheby(KG, K) : kerr_geo_inclined_sample(KG, K)
    end
    field = cheby ? :cheby_samples : :levin_trap_samples
    return _submission_geometry_lookup!(owner, field, K, owner.master) do
        KG["Trajectory"] === owner.master["Trajectory"] &&
            KG["Frequencies"] === owner.master["Frequencies"] ||
            throw(ArgumentError("Inclined sample belongs to another orbit"))
        cheby ? kerr_geo_inclined_sample_cheby(KG, K) : kerr_geo_inclined_sample(KG, K)
    end
end

# Classify (e, x) into an orbit family (1e-12 tolerance), check that the mode cache, trajectory and geometry
# owner belong together, and return (KG, master, cache, family); master is nothing without a geometry owner.
function _submission_public_route(a, p, e, x, cache, trajectory, geometry_owner)
    circular = isapprox(e, 0.0; atol=1e-12) && isapprox(abs(x), 1.0; atol=1e-12)
    family = circular ? :circular : isapprox(e, 0.0; atol=1e-12) ? :inclined :
        isapprox(abs(x), 1.0; atol=1e-12) ? :eccentric : :generic
    cache = _submission_public_cache(cache, geometry_owner, family)
    if trajectory === nothing && cache !== nothing && cache.geometry_owner !== nothing
        master = cache.geometry_owner.master
        trajectory = haskey(master, "Energy") ? master : _kg_from_presampled_master(master)
    end
    KG = _submission_public_trajectory(a, p, e, x, trajectory)
    master = _submission_public_master(cache, KG)
    return (KG=KG, master=master, cache=cache, family=family)
end

# m-independent geometry of one trapezoidal (N, K) generic sample: trapezoidal weights, the tetrad
# factors N±, L± along the radial grid, polar trigonometric factors and the rho = -1/(r - i a cos θ) tables.
# The m-dependent phases and the per-mode scratch are added in _generic_m2_m_cache.
function _generic_m2_geometry(KG_samp::Dict)
    r = KG_samp["r"]::Vector{Float64}
    rs = KG_samp["rs"]::Vector{Float64}
    theta = KG_samp["θ"]::Vector{Float64}
    N = KG_samp["N_sample"]::Int
    K = KG_samp["K_sample"]::Int
    a = KG_samp["a"]::Float64
    E = KG_samp["E"]::Float64
    Lz = KG_samp["Lz"]::Float64
    dqr = pi / (N - 1)
    dqtheta = pi / (K - 1)
    wr = _trap_weights_1d(N, dqr)
    wtheta = _trap_weights_1d(K, dqtheta)
    qr = [(i - 1) * dqr for i in 1:N]
    qtheta = [(j - 1) * dqtheta for j in 1:K]
    dphir = KG_samp["Δφr"]::Vector{Float64}
    dphitheta = KG_samp["Δφθ"]::Vector{Float64}
    Np = Vector{Float64}(undef, N)
    Nm = Vector{Float64}(undef, N)
    Np2 = Vector{Float64}(undef, N)
    Nm2 = Vector{Float64}(undef, N)
    Lp = Vector{Float64}(undef, N)
    Lm = Vector{Float64}(undef, N)
    Lp2 = Vector{Float64}(undef, N)
    Lm2 = Vector{Float64}(undef, N)
    urp = KG_samp["ur_fwd"]::Vector{Float64}
    urm = KG_samp["ur_rev"]::Vector{Float64}
    @inbounds for i in 1:N
        ri = r[i]
        r2 = ri * ri
        Delta = r2 - 2.0 * ri + a^2
        numer = E * (r2 + a^2) - a * Lz
        Np[i] = (numer + urp[i]) / Delta
        Nm[i] = (numer + urm[i]) / Delta
        Np2[i] = Np[i] * Np[i]
        Nm2[i] = Nm[i] * Nm[i]
        Lp[i] = (numer - urp[i]) / Delta
        Lm[i] = (numer - urm[i]) / Delta
        Lp2[i] = Lp[i] * Lp[i]
        Lm2[i] = Lm[i] * Lm[i]
    end
    st = Vector{Float64}(undef, K)
    ct = Vector{Float64}(undef, K)
    invst = Vector{Float64}(undef, K)
    invst2 = Vector{Float64}(undef, K)
    termM = Vector{ComplexF64}(undef, K)
    uthetap = KG_samp["uθ_fwd"]::Vector{Float64}
    uthetam = KG_samp["uθ_rev"]::Vector{Float64}
    @inbounds for j in 1:K
        stj = sin(theta[j])
        ctj = cos(theta[j])
        invstj = 1.0 / stj
        invst2j = invstj * invstj
        st[j] = stj
        ct[j] = ctj
        invst[j] = invstj
        invst2[j] = invst2j
        termM[j] = im * stj * (a * E - Lz * invst2j)
    end
    rho = Matrix{ComplexF64}(undef, N, K)
    rhobar = Matrix{ComplexF64}(undef, N, K)
    invrho = Matrix{ComplexF64}(undef, N, K)
    rho_minus = Matrix{ComplexF64}(undef, N, K)
    rho_plus = Matrix{ComplexF64}(undef, N, K)
    @inbounds for j in 1:K
        ctj = ct[j]
        for i in 1:N
            rhoij = -1.0 / (r[i] - im * a * ctj)
            rhobarij = -1.0 / (r[i] + im * a * ctj)
            rho[i, j] = rhoij
            rhobar[i, j] = rhobarij
            invrho[i, j] = 1.0 / rhoij
            rho_minus[i, j] = rhoij - rhobarij
            rho_plus[i, j] = rhoij + rhobarij
        end
    end
    return GenericM2Geometry(r, rs, theta, KG_samp["Δtr"], KG_samp["Δtθ"], dphir, dphitheta,
        KG_samp["CrossFunction"], KG_samp["initialPhases"], wr, wtheta, qr, qtheta,
        Np, Nm, Np2, Nm2, Lp, Lm, Lp2, Lm2, st, ct, invst, invst2, termM,
        uthetap, uthetam, rho, rhobar, invrho, rho_minus, rho_plus, N, K, KG_samp["Γ"], E, Lz, a)
end
