# Included inside ModeSummation. One owner is restricted to one orbit and grid.
# Clearing owners removes references; it never empties scientific return objects.

export PublicSubmission, PublicModeOwner, PublicSubmissionCancelled
export create_public_submission, ensure_public_submission, with_public_submission
export submission_orbit!, submission_master!, submission_geometry!
export with_submission_mode, with_submission_operation, request_submission_cancel!, release_public_submission!
export spawn_submission_task!

# Thrown when a mode, operation or task asks for admission after the submission closed admission; `reason` is the
# cancel reason.
struct PublicSubmissionCancelled <: Exception
    reason::String
end
Base.showerror(io::IO, e::PublicSubmissionCancelled) = print(io, "Submission admission closed: ", e.reason)

# Per-mode handle passed to lifecycle observers: the mode's private cache, its (s, l, m, n, k) tag and whether the
# cache has been cleared.
mutable struct PublicModeOwner
    cache::Any
    mode::Any
    closed::Bool
end

# One orbit's shared read-only trajectory, master sample and geometry, plus the admission, cancellation, resource and
# task bookkeeping of a point-particle submission. Create with create_public_submission, run work through
# with_submission_operation / with_submission_mode, and release with release_public_submission!.
mutable struct PublicSubmission
    parameters::Tuple
    family::Symbol
    Nmax::Int
    Kmax::Int
    context::NamedTuple
    orbit::Any
    master::Any
    geometry::Any
    alternate_masters::Dict{Symbol, Any}
    alternate_geometries::Dict{Symbol, Any}
    lock::ReentrantLock
    condition::Threads.Condition
    admission_closed::Bool
    cancel_requested::Bool
    cancel_reason::String
    active_modes::Int
    active_operations::Int
    active_by_task::IdDict{Task, Int}
    observer_depth::IdDict{Task, Int}
    accepted_tasks::IdDict{Task, Bool}
    registered_tasks::Vector{Task}
    resources::Vector{Any}
    stop_requested::Any
    lifecycle_observer::Any
    geometry_observer::Any
    releasing::Bool
    releasing_task::Union{Nothing, Task}
    released::Bool
    modes_admitted::Int
    modes_cleared::Int
    tasks_admitted::Int
    tasks_joined::Int
end

_public_family(e, x) = (x == 1 || x == -1) ? (e == 0 ? :circular : :eccentric) : (e == 0 ? :inclined : :generic)
_public_context() = (precision_bits = precision(BigFloat), rounding_mode = rounding(BigFloat),
    kg_module = KerrGeodesics, kg_version = Base.pkgversion(KerrGeodesics), kg_path = pathof(KerrGeodesics))

# New submission for one orbit (a, p, e, x) and sampling-grid caps Nmax, Kmax. Orbit, dense sample and
# geometry owner are built lazily on first use.
function create_public_submission(a, p, e, x; Nmax::Int = 2^14, Kmax::Int = 2^12,
        stop_requested = nothing, lifecycle_observer = nothing, geometry_observer = nothing)
    Nmax > 0 && Kmax > 0 || throw(ArgumentError("Submission grid caps must be positive"))
    guard = ReentrantLock()
    return PublicSubmission((a, p, e, x), _public_family(e, x), Nmax, Kmax,
        _public_context(), nothing, nothing, nothing, Dict{Symbol, Any}(), Dict{Symbol, Any}(), guard, Threads.Condition(guard),
        false, false, "closed", 0, 0, IdDict{Task, Int}(), IdDict{Task, Int}(), IdDict{Task, Bool}(), Task[], Any[],
        stop_requested, lifecycle_observer, geometry_observer, false, nothing, false, 0, 0, 0, 0)
end

function _public_check_context(scope::PublicSubmission)
    scope.context == _public_context() || throw(ArgumentError("Submission precision/rounding/dependency context changed"))
    scope.released && throw(ArgumentError("Submission is released"))
    return nothing
end

# Return (scope, owned): a new submission (owned = true) when `submission` is nothing, otherwise the given
# one after checking that it belongs to the same orbit and grid caps (owned = false).
function ensure_public_submission(submission, a, p, e, x; Nmax::Int = submission === nothing ? 2^14 : submission.Nmax, Kmax::Int = submission === nothing ? 2^12 : submission.Kmax,
        stop_requested = nothing, lifecycle_observer = nothing, geometry_observer = nothing)
    submission === nothing && return (create_public_submission(a, p, e, x;
        Nmax, Kmax, stop_requested, lifecycle_observer, geometry_observer), true)
    submission isa PublicSubmission || throw(ArgumentError("submission must be a PublicSubmission"))
    lock(submission.lock) do
        _public_check_context(submission)
        isequal(submission.parameters, (a, p, e, x)) || throw(ArgumentError("Submission belongs to another orbit"))
        (submission.Nmax == Nmax && submission.Kmax == Kmax) ||
            throw(ArgumentError("Submission belongs to another sampling grid"))
        stop_requested !== nothing && (submission.stop_requested = stop_requested)
        lifecycle_observer !== nothing && (submission.lifecycle_observer = lifecycle_observer)
        geometry_observer !== nothing && (submission.geometry_observer = geometry_observer)
    end
    return (submission, false)
end

# Called under scope.lock. Normal release drains accepted work, including its
# nested modes. Explicit cancellation rejects its next mode admission instead.
_public_admission_allowed(scope, task) = !scope.admission_closed ||
    (!scope.cancel_requested && (haskey(scope.accepted_tasks, task) || get(scope.active_by_task, task, 0) > 0))

# The lifecycle event sites below invoke callbacks after leaving scope.lock.
# This helper only holds a short lock for callback/depth bookkeeping. Geometry
# callbacks may still run under their geometry owner lock. Task-keyed depth is
# removed on normal return and error; it cannot become a permanent Task root.
function _public_observer_call(scope::PublicSubmission, select_observer, args...)
    task = current_task()
    observer = lock(scope.lock) do
        callback = select_observer()
        callback === nothing || (scope.observer_depth[task] = get(scope.observer_depth, task, 0) + 1)
        callback
    end
    observer === nothing && return nothing
    try
        Base.invokelatest(observer, args...)
    finally
        lock(scope.lock) do
            depth = scope.observer_depth[task] - 1
            depth == 0 ? delete!(scope.observer_depth, task) : (scope.observer_depth[task] = depth)
            notify(scope.condition; all = true)
        end
    end
    return nothing
end

function _public_observe(scope, event::Symbol, owner = nothing)
    return _public_observer_call(scope, () -> scope.lifecycle_observer, event, scope, owner)
end

function _public_geometry_observer(scope::PublicSubmission)
    # The geometry owner fixes its observer at construction; wrap it so its callbacks are
    # tracked in observer_depth like lifecycle callbacks.
    observer = lock(scope.lock) do
        scope.geometry_observer
    end
    observer === nothing && return nothing
    return (args...) -> _public_observer_call(scope, () -> observer, args...)
end

# Called under scope.lock, before admission or a wait can occur. Cancellation
# remains allowed from an observer; release/spawn would create self-wait cycles.
function _public_check_observer_reentry(scope::PublicSubmission, action::String)
    get(scope.observer_depth, current_task(), 0) == 0 ||
        throw(ArgumentError("Submission " * action * " cannot run from an observer callback"))
    return nothing
end

# Close admission and mark the submission cancelled; later mode admissions throw PublicSubmissionCancelled.
function request_submission_cancel!(scope::PublicSubmission; reason = "cooperative cancellation")
    changed = lock(scope.lock) do
        was_open = !scope.admission_closed
        scope.admission_closed = true
        scope.cancel_requested = true
        scope.cancel_reason = string(reason)
        notify(scope.condition; all = true)
        was_open
    end
    changed && _public_observe(scope, :admission_closed)
    return nothing
end

function _public_check_stop(scope)
    callback = scope.stop_requested
    if callback !== nothing
        requested = try
            callback()
        catch
            request_submission_cancel!(scope; reason = "stop callback raised an exception")
            rethrow()
        end
        requested === true && request_submission_cancel!(scope)
    end
    return nothing
end

# Lazily built KerrGeodesics orbit shared by all modes of the submission.
function submission_orbit!(scope::PublicSubmission)
    return lock(scope.lock) do
        _public_check_context(scope)
        !_public_admission_allowed(scope, current_task()) && get(scope.active_by_task, current_task(), 0) == 0 &&
            throw(PublicSubmissionCancelled(scope.cancel_reason))
        scope.orbit === nothing && (scope.orbit = kerr_geo_orbit(scope.parameters...))
        scope.orbit
    end
end

function _public_check_family(family::Symbol)
    family in (:circular, :eccentric, :inclined, :generic) ||
        throw(ArgumentError("Unknown submission orbit family"))
    return family
end

function _public_build_master(scope::PublicSubmission, family::Symbol)
    KG = submission_orbit!(scope)
    x = scope.parameters[4]
    return family === :generic ? GridSampling.kerr_geo_generic_sample_dense(KG, scope.Nmax, scope.Kmax) :
        family === :eccentric ? GridSampling.kerr_geo_eccentric_sample_dense(KG, scope.Nmax) :
        family === :inclined ? GridSampling.kerr_geo_inclined_sample_dense(KG, x, scope.Kmax) : KG
end

# Lazily built dense sample (master grid) of the orbit for `family`; non-default families are kept separately.
function submission_master!(scope::PublicSubmission; family::Symbol = scope.family)
    _public_check_family(family)
    return lock(scope.lock) do
        _public_check_context(scope)
        !_public_admission_allowed(scope, current_task()) && get(scope.active_by_task, current_task(), 0) == 0 &&
            throw(PublicSubmissionCancelled(scope.cancel_reason))
        if family === scope.family
            scope.master === nothing && (scope.master = _public_build_master(scope, family))
            return scope.master
        end
        return get!(scope.alternate_masters, family) do
            _public_build_master(scope, family)
        end
    end
end

function _public_build_geometry(scope::PublicSubmission, family::Symbol)
    master = submission_master!(scope; family)
    observer = _public_geometry_observer(scope)
    return family === :generic ? ConvolutionIntegrals.GenericGeometryOwner(master; observer) :
        family === :eccentric ? ConvolutionIntegrals.EccentricGeometryOwner(master; observer) :
        ConvolutionIntegrals.InclinedGeometryOwner(master; observer)
end

# Lazily built geometry owner for `family` (nothing for circular orbits), shared read-only by the mode caches.
function submission_geometry!(scope::PublicSubmission; family::Symbol = scope.family)
    _public_check_family(family)
    return lock(scope.lock) do
        _public_check_context(scope)
        !_public_admission_allowed(scope, current_task()) && get(scope.active_by_task, current_task(), 0) == 0 &&
            throw(PublicSubmissionCancelled(scope.cancel_reason))
        family === :circular && return nothing
        if family === scope.family
            scope.geometry === nothing && (scope.geometry = _public_build_geometry(scope, family))
            return scope.geometry
        end
        return get!(scope.alternate_geometries, family) do
            _public_build_geometry(scope, family)
        end
    end
end

function _public_private_cache(scope; family::Symbol = scope.family, share_geometry::Bool = true)
    _public_check_family(family)
    geometry = share_geometry ? submission_geometry!(scope; family) : nothing
    return family === :generic ? (geometry === nothing ? ConvolutionIntegrals.GenericFluxCache() : ConvolutionIntegrals.GenericFluxCache(geometry)) :
        family === :eccentric ? (geometry === nothing ? ConvolutionIntegrals.EccentricFluxCache() : ConvolutionIntegrals.EccentricFluxCache(geometry)) :
        family === :inclined ? (geometry === nothing ? ConvolutionIntegrals.InclinedFluxCache() : ConvolutionIntegrals.InclinedFluxCache(geometry)) : nothing
end

function _public_clear_cache!(scope, cache; family::Symbol = scope.family)
    cache === nothing && return nothing
    family === :generic ? ConvolutionIntegrals.clear_generic_mode_cache!(cache) :
        family === :eccentric ? ConvolutionIntegrals.clear_eccentric_mode_cache!(cache) :
        ConvolutionIntegrals.clear_inclined_mode_cache!(cache)
    return nothing
end

# Run f() as one admitted operation of the submission; resources it registered are closed on exit.
function with_submission_operation(f, scope::PublicSubmission)
    _public_check_stop(scope)
    task = current_task()
    depth = lock(scope.lock) do
        _public_check_context(scope)
        _public_admission_allowed(scope, task) || throw(PublicSubmissionCancelled(scope.cancel_reason))
        scope.active_operations += 1
        scope.active_by_task[task] = get(scope.active_by_task, task, 0) + 1
        scope.active_by_task[task]
    end
    primary_error = nothing
    try
        return f()
    catch error
        primary_error = error
        rethrow()
    finally
        cleanup_error = nothing
        try
            # File handles belong to this call, including a setup failure before
            # its inner try. Do not close another concurrent caller's handles.
            _public_close_operation_resources!(scope, task, depth)
        catch error
            cleanup_error = error
        finally
            lock(scope.lock) do
                scope.active_operations -= 1
                count = scope.active_by_task[task] - 1
                count == 0 ? delete!(scope.active_by_task, task) : (scope.active_by_task[task] = count)
                notify(scope.condition; all = true)
            end
        end
        if cleanup_error !== nothing
            primary_error === nothing ? throw(cleanup_error) : throw(CompositeException(Any[primary_error, cleanup_error]))
        end
    end
end

# Admit one mode, run f(cache) with a task-private mode cache, then clear that cache. Returns f's result.
function with_submission_mode(f, scope::PublicSubmission; mode = nothing, family::Symbol = scope.family, share_geometry::Bool = true)
    _public_check_family(family)
    _public_check_stop(scope)
    task = current_task()
    lock(scope.lock) do
        _public_check_context(scope)
        _public_admission_allowed(scope, task) || throw(PublicSubmissionCancelled(scope.cancel_reason))
        scope.active_modes += 1
        scope.active_by_task[task] = get(scope.active_by_task, task, 0) + 1
        scope.modes_admitted += 1
    end
    owner = PublicModeOwner(nothing, mode, false)
    succeeded = false
    primary_error = nothing
    result = nothing
    try
        owner.cache = _public_private_cache(scope; family, share_geometry)
        _public_observe(scope, :mode_admitted, owner)
        result = f(owner.cache)
        succeeded = true
        _public_observe(scope, :mode_complete, owner)
        return result
    catch error
        primary_error = error
        rethrow()
    finally
        # A completed kernel result can still be unreturned if its observer
        # throws. Detach this local reference; preserve the evaluated return
        # value on normal return and never modify Task.result or scientific data.
        result = nothing
        cleanup_errors = Any[]
        try
            succeeded || _public_observe(scope, :mode_failed, owner)
        catch error
            push!(cleanup_errors, error)
        end
        try
            _public_clear_cache!(scope, owner.cache; family)
            owner.closed = true
            _public_observe(scope, :mode_cache_cleared, owner)
        catch error
            push!(cleanup_errors, error)
        finally
            owner.cache = nothing
            lock(scope.lock) do
                scope.active_modes -= 1
                count = scope.active_by_task[task] - 1
                count == 0 ? delete!(scope.active_by_task, task) : (scope.active_by_task[task] = count)
                owner.closed && (scope.modes_cleared += 1)
                notify(scope.condition; all = true)
            end
        end
        if !isempty(cleanup_errors)
            primary_error === nothing || pushfirst!(cleanup_errors, primary_error)
            throw(CompositeException(cleanup_errors))
        end
    end
end

Base.@noinline function _public_task_worker(payload::Base.RefValue{Any})
    bundle = payload[]
    try
        return bundle.f()
    finally
        scope = bundle.scope
        lock(scope.lock) do
            delete!(scope.accepted_tasks, current_task())
            scope.active_operations -= 1
            notify(scope.condition; all = true)
        end
        scope = nothing
        bundle = nothing
        payload[] = nothing
    end
end

# Spawn f() as a task registered with the submission; release joins it before clearing caches.
function spawn_submission_task!(f, scope::PublicSubmission)
    lock(scope.lock) do
        _public_check_observer_reentry(scope, "spawn")
    end
    _public_check_stop(scope)
    task = lock(scope.lock) do
        _public_check_observer_reentry(scope, "spawn")
        _public_check_context(scope)
        scope.admission_closed && throw(PublicSubmissionCancelled(scope.cancel_reason))
        payload = Ref{Any}((f = f, scope = scope))
        task = Threads.@spawn _public_task_worker($payload)
        push!(scope.registered_tasks, task)
        scope.accepted_tasks[task] = true
        scope.active_operations += 1
        scope.tasks_admitted += 1
        task
    end
    # Registration is committed before notification. If the observer throws,
    # release still owns/drains this task even though spawn did not return it.
    _public_observe(scope, :task_admitted)
    return task
end

function _public_track_resource!(scope, resource)
    resource === nothing && return nothing
    lock(scope.lock) do
        task = current_task()
        push!(scope.resources, (task = task, depth = get(scope.active_by_task, task, 0), resource = resource))
    end
    return resource
end

function _public_close_resource!(scope, resource)
    resource === nothing && return nothing
    close(resource)
    lock(scope.lock) do
        filter!(item -> item.resource !== resource, scope.resources)
    end
    return nothing
end

function _public_close_operation_resources!(scope, task, depth)
    resources = lock(scope.lock) do
        [item.resource for item in scope.resources if item.task === task && item.depth == depth]
    end
    errors = Any[]
    for resource in resources
        try
            _public_close_resource!(scope, resource)
        catch error
            push!(errors, error)
        end
    end
    isempty(errors) || throw(CompositeException(errors))
    return nothing
end

function _public_drain_task!(task::Task, errors)
    # Even an interruption of the parent's wait cannot authorize early release.
    # A failed scientific Task is done; preserve its error after partner drain.
    while true
        try
            fetch(task)
            return nothing
        catch error
            push!(errors, error)
            istaskdone(task) && return nothing
        end
    end
end

function _public_wait_operations!(scope, errors)
    lock(scope.lock) do
        while scope.active_modes != 0 || scope.active_operations != 0 || !isempty(scope.observer_depth)
            try
                wait(scope.condition)
            catch error
                push!(errors, error)
            end
        end
    end
    return nothing
end

# Close admission, join registered tasks, wait for active operations, then release geometry owners and
# resources. Idempotent; errors collected on the way are rethrown together as a CompositeException.
function release_public_submission!(scope::PublicSubmission)
    tasks = lock(scope.lock) do
        _public_check_observer_reentry(scope, "release")
        scope.released && return nothing
        get(scope.active_by_task, current_task(), 0) == 0 ||
            throw(ArgumentError("Release must run outside an active mode; its parent must join workers"))
        any(task -> task === current_task(), scope.registered_tasks) &&
            throw(ArgumentError("A registered worker cannot release its own submission"))
        if scope.releasing
            scope.releasing_task === current_task() && throw(ArgumentError("Reentrant submission release cannot wait on itself"))
            while !scope.released
                wait(scope.condition)
            end
            return nothing
        end
        scope.releasing = true
        scope.releasing_task = current_task()
        scope.admission_closed = true
        copy(scope.registered_tasks)
    end
    tasks === nothing && return nothing
    errors = Any[]
    geometries = nothing
    geometry = nothing
    task = nothing
    try
        try
            _public_observe(scope, :admission_closed)
        catch error
            push!(errors, error)
        end
        # Fetch every admitted task even if the first one failed. Scientific
        # failures remain visible only after all partners are finished.
        for task in tasks
            try
                _public_drain_task!(task, errors)
            finally
                lock(scope.lock) do
                    scope.tasks_joined += 1
                end
                try
                    _public_observe(scope, :task_joined)
                catch error
                    push!(errors, error)
                end
            end
        end
        _public_wait_operations!(scope, errors)
        geometries = lock(scope.lock) do
            owners = Pair{Symbol, Any}[]
            scope.geometry === nothing || push!(owners, scope.family => scope.geometry)
            append!(owners, collect(scope.alternate_geometries))
            owners
        end
        for (family, geometry) in geometries
            try
                family === :generic ? ConvolutionIntegrals.release_generic_geometry!(geometry) :
                    family === :eccentric ? ConvolutionIntegrals.release_eccentric_geometry!(geometry) :
                    ConvolutionIntegrals.release_inclined_geometry!(geometry)
            catch error
                push!(errors, error)
            end
        end
        for item in copy(scope.resources)
            try
                _public_close_resource!(scope, item.resource)
            catch error
                push!(errors, error)
            end
        end
    finally
        # Retry drain here if an unforeseen wrapper failure skipped its normal
        # path. There is no kernel interruption or early shared-owner clearing.
        for task in tasks
            istaskdone(task) || _public_drain_task!(task, errors)
        end
        _public_wait_operations!(scope, errors)
        # Parents/tasks may legitimately retain return data. Remove this owner's
        # references only, after drain; never deep-empty KG, YSolution or SWSH.
        try
            lock(scope.lock) do
                scope.geometry = nothing
                scope.master = nothing
                scope.orbit = nothing
                empty!(scope.alternate_masters)
                empty!(scope.alternate_geometries)
                empty!(scope.resources)
                empty!(scope.registered_tasks)
                empty!(scope.accepted_tasks)
                scope.released = true
                scope.releasing = false
                scope.releasing_task = nothing
                notify(scope.condition; all = true)
            end
        finally
            # Drain is finished. Drop this wrapper's local Task/owner lists even
            # if cleanup or the eventual CompositeException does not return.
            # User-owned Task.result and scientific payloads remain untouched.
            task = nothing
            tasks = nothing
            geometry = nothing
            geometries = nothing
        end
    end
    try
        _public_observe(scope, :submission_released)
    catch error
        push!(errors, error)
    end
    isempty(errors) || throw(CompositeException(errors))
    return nothing
end

# Create a submission, call f(scope) and release the submission on exit.
function with_public_submission(f, a, p, e, x; kwargs...)
    scope = create_public_submission(a, p, e, x; kwargs...)
    return _with_public_scope(scope, true) do
        f(scope)
    end
end

function _with_public_scope(f, scope, owned::Bool)
    primary_error = nothing
    try
        return f()
    catch error
        primary_error = error
        rethrow()
    finally
        if owned
            try
                release_public_submission!(scope)
            catch cleanup_error
                primary_error === nothing ? rethrow() : throw(CompositeException(Any[primary_error, cleanup_error]))
            end
        end
    end
end
