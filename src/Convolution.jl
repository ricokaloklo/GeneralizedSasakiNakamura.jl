module Convolution

using QuadGK

function convolution_integral(s::Int, a, R_homo, sourceterm_without_phasing, sourceterm_phasing, lower_limit, upper_limit)
    f(r) = R_homo(r) * Delta(a, r)^s * sourceterm_without_phasing(r)
    q(r) = sourceterm_phasing(r)

    integral, _ = quadgk(r -> f(r)*exp(1im*q(r)), lower_limit, upper_limit)
    return integral
end

end