# Reader for the BHPT reference files written by bhpt_acceptance_reference.wls ("case|key|value").
# `load_bhpt_refs(path)` returns Dict{String,Dict{String,Complex{BigFloat}}} (case tag -> key -> value);
# values are parsed at the current BigFloat precision from the printed digits (precision marks removed).
function _parse_mathematica_number(text::AbstractString)
    # "0``acc" is a number known only to be below 10^-acc in absolute value (the paclet's zero with an
    # accuracy mark), not an exact zero: no usable reference.
    occursin(r"(^|[^0-9.])0``", strip(text)) && return nothing
    t = replace(strip(text), r"`+-?[0-9.]*" => "")     # drop precision marks  `41.7, ``55.0, ``-6.8
    t = replace(t, "*^" => "e")
    t = replace(t, " " => "")
    occursin("Indeterminate", t) && return nothing
    occursin("ComplexInfinity", t) && return nothing
    occursin("ABORTED", t) && return nothing
    real_part = "0"; imag_part = "0"
    if endswith(t, "*I")
        body = t[1:end-2]
        # split at the last sign that does not belong to an exponent
        split_at = 0
        for i in length(body):-1:2
            if (body[i] == '+' || body[i] == '-') && body[i-1] != 'e'
                split_at = i
                break
            end
        end
        if split_at == 0
            imag_part = body
        else
            real_part = body[1:split_at-1]
            imag_part = body[split_at:end]
        end
    else
        real_part = t
    end
    tobig(s) = begin
        s = s == "+" ? "1" : s == "-" ? "-1" : s
        if occursin("/", s)
            num, den = split(s, "/")
            BigFloat(parse(BigInt, num)) / BigFloat(parse(BigInt, den))
        else
            BigFloat(s)
        end
    end
    return Complex{BigFloat}(tobig(real_part), tobig(imag_part))
end

function load_bhpt_refs(path::AbstractString)
    refs = Dict{String,Dict{String,Union{Nothing,Complex{BigFloat}}}}()
    for line in eachline(path)
        parts = split(line, "|")
        length(parts) == 3 || continue
        tag, key, value = parts
        key == "seconds" && continue
        d = get!(refs, String(tag), Dict{String,Union{Nothing,Complex{BigFloat}}}())
        d[String(key)] = _parse_mathematica_number(value)
    end
    return refs
end

# tag of a case, as the .wls writes it (a, omega as Mathematica InputForm of exact rationals)
function bhpt_tag(a_str::AbstractString, l::Int, m::Int, omega_str::AbstractString)
    return "a=$a_str;l=$l;m=$m;omega=$omega_str"
end

digits_agreed(x, ref) = ref === nothing ? NaN : (x == ref ? Inf : -log10(Float64(abs(x - ref) / abs(ref))))
