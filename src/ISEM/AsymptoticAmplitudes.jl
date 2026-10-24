module AsymptoticAmplitudes

using ..Matching: Route, ConjugatedRoute

@inline _finite_complex(value) =
    isfinite(real(value)) && isfinite(imag(value))

# Incidence and reflection amplitudes divided by the route's internal solution
# scale (unit GSN transmission); `nothing` when either is not finite in Float64.
function _try_unit_pair(route::Route)
    scale = route.solution_scale
    isfinite(scale) && scale > 0 ||
        error("invalid ISEM GSN internal solution scale")
    incidence = route.incidence / scale
    reflection = route.reflection / scale
    _finite_complex(incidence) && _finite_complex(reflection) || return nothing
    return ComplexF64(incidence), ComplexF64(reflection)
end

function _try_unit_pair(route::ConjugatedRoute)
    pair = _try_unit_pair(route.route)
    pair === nothing && return nothing
    return conj(pair[1]), conj(pair[2])
end

end
