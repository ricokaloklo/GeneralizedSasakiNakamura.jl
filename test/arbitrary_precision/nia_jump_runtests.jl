# The jump of the incidence amplitude across the negative imaginary axis, from the two one-sided limits
# of the arbitrary-precision route (omega = +-1e-30 - i sigma; the -delta side through the mirror mode
# (m, omega) -> (-m, -conj omega)), against the analytic MST branch strength q of the Float64 module
# (`mano_suzuki_takasugi_negative_imaginary_axis_jump`), in the Teukolsky unit-transmission normalization:
#     B_inc(-delta) = B_inc(+delta) + i q B_ref(+delta),      B_ref continuous.
# The strength is a binary64 quantity certified to about 1e-13, so the jump relation is required to that
# accuracy (100 times the module's own error estimate, at least 1e-13), the continuity of the reflection to
# 1e-28 at 113 bits. Cases with an uncertified strength (sigma = 3) are reported, not tested. Opt-in.
#
#     GSN_SOURCE=<checkout> julia --project=<env> test/arbitrary_precision/nia_jump_runtests.jl

using Test
using Printf

const _AP_SOURCE = get(ENV, "GSN_SOURCE", "")
if isempty(_AP_SOURCE)
    using GeneralizedSasakiNakamura
else
    include(joinpath(_AP_SOURCE, "src", "GeneralizedSasakiNakamura.jl"))
    using .GeneralizedSasakiNakamura
end
const GSN = GeneralizedSasakiNakamura
const MST = GSN.ISEM.ManoSuzukiTakasugiInfinity

relative(x, y) = Float64(abs(x - y) / abs(y))

@testset "negative-imaginary-axis jump vs MST branch strength" begin
    setprecision(BigFloat, 113) do
        delta = BigFloat("1e-30")
        for (a, l, m) in ((2 // 5, 2, 2), (2 // 5, 3, 1), (9 // 10, 2, 2)), sigma in (1 // 10, 1 // 2, 1 // 1, 3 // 1)
            ab = BigFloat(a)
            sb = BigFloat(sigma)
            t0 = time()
            plus = GSN.Teukolsky_radial(-2, l, m, ab, Complex{BigFloat}(delta, -sb), GSN.IN)
            minus = GSN.Teukolsky_radial(-2, l, m, ab, Complex{BigFloat}(-delta, -sb), GSN.IN)
            lambda = GSN.ISEM.Parameters.gsn_parameters(-2, l, m, ab, Complex{BigFloat}(0, -sb)).lambda
            jump = MST.mano_suzuki_takasugi_negative_imaginary_axis_jump(-2, l, m, Float64(ab), ComplexF64(0, -Float64(sb)),
                ComplexF64(lambda), ComplexF64(plus.incidence_amplitude), ComplexF64(plus.reflection_amplitude);
                require_certificate=false)
            q = jump.up_jump_coefficient / im
            predicted = plus.incidence_amplitude + im * q * plus.reflection_amplitude
            jump_error = relative(minus.incidence_amplitude, predicted)
            continuity = relative(minus.reflection_amplitude, plus.reflection_amplitude)
            allowed = 100 * max(jump.estimated_relerr, 1e-13)
            @printf("a=%-4s l=%d m=%d sigma=%-4s  B_inc(-) vs B_inc(+) + i q B_ref: %.1e (allowed %.1e, q certified %s)  B_ref continuity %.1e  %5.1f s\n",
                string(a), l, m, string(sigma), jump_error, allowed, jump.certificate_accepted, continuity, time() - t0)
            @test continuity <= 1e-28
            if jump.certificate_accepted
                @test jump_error <= allowed
            else
                @test_skip jump_error <= allowed
            end
        end
    end
end
