# Arbitrary precision

The ISEM routes solve the homogeneous GSN/Teukolsky radial problem in the floating type of their
inputs. With `BigFloat` spin and frequency the IN and UP solutions, their asymptotic amplitudes and
the radial functions are computed to the working precision of `BigFloat`, for real and complex
frequencies; the Float64 inputs keep the Float64 routes unchanged.

```julia
using GeneralizedSasakiNakamura
setprecision(BigFloat, 256)
a = BigFloat("0.4"); omega = BigFloat("0.5")               # 4/10 and 1/2 to 256 bits
Rin, Rup = Teukolsky_radial(-2, 2, 2, a, omega)             # unit Teukolsky transmission
Rin.incidence_amplitude                                       # Complex{BigFloat}, ~74 digits
R, dRdr, d2Rdr2, err = Rin.Teukolsky_solution(BigFloat("3.5"))
Xin = GSN_radial(-2, 2, 2, a, omega, IN)                     # unit GSN transmission
Xin(BigFloat(10))                                            # (X, dX/dr*, 0) at r* = 10
Rc = Teukolsky_radial(-2, 2, 2, a, Complex{BigFloat}(BigFloat("0.5"), BigFloat("-0.3")), IN)
```

A binary64 value among the inputs is used exactly: `0.4` is the binary double closest to 4/10,
`BigFloat("0.4")` is 4/10 to the current precision. The separation constant is taken from
`SpinWeightedSpheroidalHarmonics` at the same precision, which requires a version of that package
that accepts `BigFloat` at complex `c` (v1.4.1 or later). `data_type=Complex{BigFloat}` with
`method="auto"` or `"ISEM"` promotes Float64 inputs to `BigFloat`. Quasinormal frequencies and the
frequency derivative of the incidence amplitude at arbitrary precision are provided by
`QNM.qnm_arbitrary_precision` (see [Quasinormal modes](QNM.md)).

## Fixed-precision wide types (Double64)

The routes also run in a fixed-precision type wider than binary64, such as `Double64` of
DoubleFloats.jl (106 bits, 31.9 digits), which is not a dependency of the package: `Teukolsky_radial(-2, 2, 2,
Double64(2) / 5, Double64(1) / 2)` returns `Complex{Double64}` amplitudes and radial values. The separation
constant comes from SpinWeightedSpheroidalHarmonics in that type, the series and the Taylor patches are
summed in it, the MST data of a complex frequency are built in BigFloat at the type's precision plus the
MST guard (the recurrences of the MST sums need more bits than the type has) and converted on evaluation,
and the asymptotic series is the infinity side at real frequencies (20 times faster than MST here). The
working precision cannot be raised in such a type, so the route is built once with a tolerance of 64 ulps
(``3 \times 10^{-30}``) and no guard; when the measured cancellation of the matching and of the series sums
exceeds 12 bits (an exponentially small reflection, a strongly damped frequency, the horizon chain at
``a = 0.9``, ``\omega = 3i``), the route is rebuilt in BigFloat at the type's precision with the full rebuild
loop, and that route is returned, with BigFloat values. Over the 178 reference cases (`GSN_AP_BITS=double64`)
the amplitudes agree with the paclet to 26.4-32.9 digits and the radial values to 27.2-32.7, all passing
the ``N\log_{10}2 - 6 = 25.9`` requirement; a real-frequency route costs 10-15 ms, a complex one 0.2-1.2 s
(the MST data in BigFloat), against 0.3-0.7 s for BigFloat at 113 bits and 0.4-0.8 ms for Float64.

## What is computed and how

The route is the ISEM construction in the compact coordinate ``x = (r - r_+)/(r - r_-)``: the
horizon Frobenius series, Taylor patches along the real axis, and a representation of the two
solutions at infinity, matched at an interior point. All series are the exact generated coefficient
tables of the Float64 routes, evaluated in the working type.

- **Tolerance and orders.** The relative truncation tolerance of every series is 64 ulps of the
  inputs' precision unless `tol` is given; the Taylor order of the patches is 2.2 per requested
  decimal digit (`N` overrides it), the horizon series uses twice and the asymptotic series three
  times that order. The patch step is at most half the empirical radius of the solution series
  (ratio of its last coefficients) and is shrunk until the last four terms are below the tolerance.
- **Working precision.** The BigFloat precision is raised internally above the bits needed for the
  tolerance by 24 guard bits, by the cancellation of a Taylor step (the tail criterion lets a step
  reach ``k h \approx N\,\mathrm{tol}^{1/N}/e`` for a solution of local wavenumber ``k``, whose terms
  grow to about ``e^{kh}`` before they cancel to the value: 13 bits at 113 bits, 31 bits at 256 bits),
  and by the continuation loss below when the asymptotic seed is used at a complex frequency. The
  endpoint series are seeded where their own cancellation stays below the same bound. The derived
  parameters ``\kappa = \sqrt{1 - a^2}`` and ``\nu = \arctan(a/\kappa)`` are recomputed at the working
  precision: the coefficient tables take ``\nu``, the local solutions ``\kappa`` and ``a``, and rounded
  separately at the input precision they describe two slightly different equations, so that every
  solution carries a generic error of that size, which the matching below amplifies. Evaluations of
  the stored series and the GSN-to-Teukolsky conversion run at the working precision as well (their
  sums cancel by the same factors), and the returned amplitudes and radial values keep it.
- **Real frequencies.** Two representations of the infinity side exist: the asymptotic series
  ``e^{\sigma/y} y^{\rho} \sum d_n y^n`` (``y = 1 - x``), optimally truncated at the largest ``y`` where a
  term falls below the tolerance and continued inward by Taylor patches (`infinity=:series`), and the
  MST series below (`infinity=:mst`). Both give the same digits (at 256 bits the IN and UP amplitude
  ratios of ``a = 0.4, s = -2, l = m = 2, \omega = 1/2`` agree with the Black Hole Perturbation Toolkit
  Teukolsky paclet to 74 digits either way); the MST side is 1.5-4 times faster at 256 bits once its
  series factors are stored, so `:auto` takes it at every frequency and keeps the asymptotic series as
  the fallback at a real frequency.
- **Complex frequencies.** The asymptotic series cannot define a single branch at infinity when
  ``\operatorname{Im}\omega \ne 0``: the truncated series leaves an admixture of the recessive
  solution of relative size ``e^{-2|\omega| r_0}`` at the seed radius ``r_0``, while the recessive
  solution is ``e^{-2|\operatorname{Im}\omega| r_0}`` smaller than the dominant one, so the branch
  coefficient is wrong by ``e^{-2(|\omega| - |\operatorname{Im}\omega|) r_0}``; continuing the
  inward-decaying branch from ``r_0`` to the matching radius ``r_m`` then loses a factor
  ``e^{2|\operatorname{Im}\omega|(r_0 - r_m)}``. With ``r_0 \sim \ln(1/\mathrm{tol})/(2|\omega|)``
  the result is only ``\mathrm{tol}^{1 - |\operatorname{Im}\omega|/|\omega|}`` accurate, and nothing
  helps on the negative imaginary axis. The route therefore uses the Mano-Suzuki-Takasugi (MST)
  series for the UP and DOWN solutions at complex frequency: the renormalized angular momentum is
  solved in the working type from the Float64 value, the series coefficients come from the
  continued fractions, the Tricomi functions ``c^n U(a + n, b + 2n, c)`` are generated by their
  three-term recurrence in ``n`` at a precision raised by 96 bits, re-anchored by direct Kummer-series
  evaluations (themselves at a precision raised by ``|c|/\ln 2`` bits for the cancellation of the
  series) whenever the accumulated condition of the recurrence exceeds 56 bits, and the two Teukolsky
  solutions are converted to GSN states at the matching point ``x = 1/2``. At 256 bits the IN and UP
  amplitude ratios at ``\omega = 1/2 - 0.3i`` agree with the paclet to the certified digits of the
  reference; at ``2 - 2i`` to 46-72 digits, where the paclet's own 120-digit run is the limit.
- **Matching conditioning.** The connection coefficients are ratios of ``2 \times 2`` determinants
  whose terms cancel when the two basis solutions differ in size at the matching point: at complex
  frequency the infinity branches differ by ``e^{2|\operatorname{Im}\omega| r^*_m}``, which is why the
  matching stays at ``x = 1/2`` rather than far out, and an exponentially small reflection costs
  ``|B_{\rm inc}/B_{\rm ref}|`` (``10^{54}`` for the IN solution at ``\omega = 10``, ``a = 0.4``,
  ``l = m = 2``, where ``|B_{\rm ref}/B_{\rm trans}| = 3 \times 10^{-53}``). The cancellation actually
  met, ``(|u| + |v|)/|u - v|``, is recorded as `matching_condition` in the solution's metadata, the
  largest cancellation of any series sum as `series_condition`, and the route is rebuilt with the
  series tolerance tightened by the measured bits (plus 8) and the working precision raised to cover
  both. The measurement is repeated after every rebuild (`rebuilds` counts them), because a build
  cannot see a cancellation beyond its own noise, about ``1/\mathrm{tol}``: at ``\omega = 10`` the
  first build at 113 bits measures ``10^{34}``, the second ``4 \times 10^{55}``, and the third, at a
  tolerance of ``10^{-90}`` and 361 working bits, returns the reflection to 32 digits in 5 s. Without
  the loop the reflection came out as ``10^{-31}``, pure noise, and its measured sensitivity to the
  inputs was the noise as well.
- **Negative real part and the negative imaginary axis.** The principal-branch MST formulas are
  not on the physical sheet for ``\operatorname{Re}\omega < 0``; those frequencies are solved for the
  mirror mode ``(m, \omega) \to (-m, -\bar\omega)`` and complex-conjugated, as the Float64 routes do
  for real negative frequencies. The negative imaginary axis is a branch cut of the physical
  amplitudes: the two-sided limits ``\omega = \pm\delta - i\sigma`` are obtained directly (the
  ``+\delta`` side and the axis point itself on the principal branch, the ``-\delta`` side through the
  mirror mode), and they differ by O(1) already at ``\sigma = 0.1``. The Float64 rotated-contour
  route converges to both limits as ``\delta \to 0``; the paclet's principal branch, which is
  continuous across the axis, gives the ``+\delta`` value on both sides, so its values at
  ``\operatorname{Re}\omega < 0`` are not references. The jump itself is the analytic MST one: with the
  branch strength ``q`` of `mano_suzuki_takasugi_negative_imaginary_axis_jump` (a binary64 quantity
  certified to about ``10^{-13}``), ``B_{\rm inc}(-\delta) = B_{\rm inc}(+\delta) + i q B_{\rm ref}(+\delta)``
  in Teukolsky unit-transmission normalization holds to ``10^{-13}``-``10^{-16}`` and the reflection is
  continuous to ``10^{-29}`` at 113 bits, for ``a = 0.4, 0.9``, ``(l, m) = (2, 2), (3, 1)``,
  ``\sigma = 0.1, 0.5, 1`` (test/arbitrary_precision/nia_jump_runtests.jl); at ``\sigma = 3`` the Float64
  strength is uncertified (``10^{-5}``-``10^{-7}``) and the relation holds to that accuracy.
- **Horizon side.** The IN solution is always the horizon Frobenius series continued by Taylor
  patches (the "hybrid" of the design brief): at ``x = 1/2`` that side costs a few patches and a
  small fraction of the build, so a pure-MST IN solution (the hypergeometric series, whose convergence
  radius ends exactly at ``x = 1/2``, i.e. ``r = 1 + 3\kappa``) was not implemented.
- **Not covered.** The large-frequency and near-extremal regimes of the Float64 routes rest on
  binary64 tables and are not used; the Taylor patches cover ``|\omega| \le 10`` at the cost of more
  patches. The superradiant threshold (logarithmic horizon companion) is not built for the UP route.

- **Evaluation outside the matching point.** With the MST side the two branches are continued
  outward by Taylor patches to ``r_m + 25/|\operatorname{Im}\omega|`` (at most 120), built at a
  precision raised by ``2|\operatorname{Im}\omega|(r^*_{\rm end} - r^*_m)/\ln 2`` bits so that the
  branch that decays outward keeps the tolerance despite the admixture of the growing one that
  rounding leaves; a request beyond the chain extends it on demand, to twice the distance from the
  matching radius, at a precision raised by the growth over the new stretch, as long as the point lies
  within 60 radians of phase of the matching radius and the extension costs at most 64 extra bits;
  farther out the MST series are summed directly at the point (about 0.1 s per point at 256 bits).
  Inside the chain a radial function is evaluated in a fraction of a millisecond per point.
- **Conditioning with respect to the inputs.** `BigFloat("0.4")` at ``N`` bits is not 4/10; the
  route solves the problem defined by the rounded inputs exactly to its tolerance, and a quantity
  that is ill-conditioned with respect to ``a`` and ``\omega`` (an incidence amplitude in the strongly
  damped region near the quasinormal frequencies, a small reflection amplitude) differs from the
  exact-rational answer by its condition number times ``2^{-N}``. The acceptance tests measure this
  loss by one-ulp perturbations when a quantity misses its plain requirement, and report it. The
  measurement is meaningful only once the route's own cancellations are covered: a value that is
  noise changes completely under any perturbation.

## Measured

Against the paclet's 120- and 200-digit runs (200 and 400 digits for ``\omega = 10``), which agree
with each other to 85-99 digits on the amplitudes and 54-60 digits on the radial values, for
``a = 0.4`` and ``0.9``, ``(l, m) = (2, 2)`` and ``(3, 1)``, ``s = -2``, 800 checks
(test/arbitrary_precision, `GSN_AP_CASES=all`; digits agreed, range over the cases):

| quantity | 113 bits (34.0 digits) | 256 bits (77.1 digits) |
|---|---|---|
| Inc/Trans, Ref/Trans of IN and UP, ``\omega \in \{1/20, 1/2, 2\}`` | 31.2-33.4 | 73.5-75.5 |
| same, ``\omega \in \{1/2 \pm 0.3i, 0.3 - 1.5i, 2 - 2i\}`` | 31.1-34.2 | 74.4-77.3 |
| same, ``\omega = \delta \pm i\sigma`` (``\delta`` from ``10^{-6}`` to ``10^{-30}``, ``\sigma = 1/10, 1/2, 1, 3``) and ``i/10`` | 30.4-35.3 | 73.6-78.1 |
| same, ``\omega = 10`` (reflections ``10^{-53}`` to ``10^{-62}``) | 29.4-34.3 | 73.7-76.8 |
| ``R(r')``, ``dR(r')``, all of the above | 29.6-35.8 | 72.4-77.9 |

The whole suite was also run on an independent machine (Windows 11, Julia 1.12.6, this commit, the
`ArbitraryPrecision` SWSH at 7cbc9b1) at 64, 113, 192 and 256 bits: 1424 checks per precision (2848
for the 113/256 run), all passing, with digits agreed of 14.1-19.4 at 64 bits (19.3 digits),
29.4-35.2 at 113, 53.3-59.2 at 192 (57.8 digits) and 72.4-78.4 at 256, the lowest values being the
``\omega = 10`` reflections and the ``\sigma = 3`` axis cases; the 22 quasinormal roots passed at 113 and
256 bits (88 checks). The ``-\delta`` side of the negative axis is the mirror of the ``+\delta`` side
of the mode ``-m`` and is not in the paclet references (see above); the exact point ``-i/10`` is a
branch-cut point whose paclet value is a principal-value diagnostic, not a reference.

The default tolerances are 64 ulps, i.e. 32.2 and 75.1 digits. With the MST side a route costs
0.3-1.2 s at 113 bits (4-6 s at ``\omega = 10`` with its three builds) and 0.6-2.4 s at 256 bits, real or
complex frequency alike (the asymptotic-series side costs 1.5-4 times more at 256 bits: 2.7-5.7 s); a
Float64 route costs 0.4-0.8 ms, so the factor is about 500 at 113 bits, the price of BigFloat arithmetic
rather than of the algorithm. Evaluating a radial function costs 0.2-0.5 ms per point inside the
matching radius; the first evaluation beyond it builds the outward chain (1-5 s). A quasinormal-mode
polish takes 15-40 s at 113 and 256 bits. The quasinormal frequencies of the ``(2, 2, 0)`` modes at
``a = 0.4`` and ``0.9`` polished at 128 bits agree with the 192-bit ones to 34 digits.

## Precision bookkeeping

The acceptance tests in `test/arbitrary_precision` compare every stored paclet case at the requested
working bits and print, per quantity, the digits agreed, the digits expected (``N\log_{10}2 - 6``),
the input conditioning where it was measured, and the certified accuracy of the reference (the
agreement between two paclet working precisions). Reported precision is always the agreement with
that independent reference, never a self-consistency figure between tolerances.
