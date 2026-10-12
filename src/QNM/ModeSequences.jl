"""
    qnm_sequence(spins, s, l, m, n, branch=ordinary; guess=nothing, kwargs...)
    qnm_sequence(mode::QNMMode, spins; guess=nothing, kwargs...)

Solve a labeled QNM along the supplied spin sequence. Each accepted result is
the explicit seed for the next spin; the sequence stops at the first failure.
Keywords are those of [`qnm_frequency`](@ref) (`convention`, `segment`, ...).
"""
function qnm_sequence(spins, s::Integer, l::Integer, m::Integer, n::Integer,
        branch=ordinary; guess=nothing, kwargs...)
    mode = _qnm_flat_mode(s, l, m, n, branch)
    results = LeaverResult[]
    current_guess = guess
    for spin in collect(spins)
        result = current_guess === nothing ?
            qnm_frequency(spin, s, l, m, n, branch; kwargs...) :
            qnm_frequency(mode, spin; guess=current_guess, kwargs...)
        push!(results, result)
        result.status == :accepted || break
        current_guess = result.omega
    end
    return results
end

function qnm_sequence(mode::QNMMode, spins; guess=nothing, kwargs...)
    values = collect(spins)
    isempty(values) && return LeaverResult[]
    results = LeaverResult[]
    current_guess = guess
    for spin in values
        result = qnm_frequency(mode, spin; guess=current_guess, kwargs...)
        push!(results, result)
        result.status == :accepted || break
        current_guess = result.omega
    end
    return results
end
