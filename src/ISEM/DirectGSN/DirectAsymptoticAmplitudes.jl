module DirectAsymptoticAmplitudes

using ..DirectMatching: DirectRoute, DirectConjugatedRoute

export direct_asymptotic_amplitudes, direct_amplitude_pair
export branch_pair_relative_error

@inline _finite_complex(value) =
    isfinite(real(value)) && isfinite(imag(value))

# Incidence and reflection amplitudes divided by the route's internal solution
# scale (unit GSN transmission); `nothing` when either is not finite in Float64.
function _try_unit_pair(route::DirectRoute)
    scale = route.solution_scale
    isfinite(scale) && scale > 0 ||
        error("invalid direct GSN internal solution scale")
    incidence = route.incidence / scale
    reflection = route.reflection / scale
    _finite_complex(incidence) && _finite_complex(reflection) || return nothing
    return ComplexF64(incidence), ComplexF64(reflection)
end

function _try_unit_pair(route::DirectConjugatedRoute)
    pair = _try_unit_pair(route.route)
    pair === nothing && return nothing
    return conj(pair[1]), conj(pair[2])
end

function _unit_pair(route::Union{DirectRoute,DirectConjugatedRoute})
    pair = _try_unit_pair(route)
    pair === nothing &&
        error("unit-transmission direct GSN amplitudes exceed Float64 range")
    return pair
end

# Unit-transmission amplitudes as a named tuple, labelled Binc/Bref for the IN
# branch and Cinc/Cref for the UP branch.
function direct_asymptotic_amplitudes(
    route::Union{DirectRoute,DirectConjugatedRoute},
)
    incidence, reflection = _unit_pair(route)
    if route.branch == :IN
        return (
            branch=:IN,
            convention=:unit_gsn_transmission,
            incidence_name=:Binc,
            reflection_name=:Bref,
            incidence=incidence,
            reflection=reflection,
            Binc=incidence,
            Bref=reflection,
        )
    elseif route.branch == :UP
        return (
            branch=:UP,
            convention=:unit_gsn_transmission,
            incidence_name=:Cinc,
            reflection_name=:Cref,
            incidence=incidence,
            reflection=reflection,
            Cinc=incidence,
            Cref=reflection,
        )
    end
    throw(ArgumentError("unsupported direct GSN route branch $(route.branch)."))
end

direct_amplitude_pair(route::Union{DirectRoute,DirectConjugatedRoute}) =
    _unit_pair(route)

# (incidence, reflection) from a route, an amplitude named tuple, or any 2-element collection.
function direct_amplitude_pair(value)
    if value isa NamedTuple && haskey(value, :incidence) && haskey(value, :reflection)
        return (value.incidence, value.reflection)
    end
    return Tuple(value)
end

# Largest componentwise difference of two amplitude pairs, relative to the
# largest amplitude magnitude in either pair.
function branch_pair_relative_error(candidate, reference)
    cand = direct_amplitude_pair(candidate)
    ref = direct_amplitude_pair(reference)
    scale = maximum(abs, (cand..., ref...))
    scale == 0 && return 0.0
    return maximum(abs.(cand .- ref)) / scale
end

end
