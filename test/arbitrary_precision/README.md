# Arbitrary-precision acceptance tests

Opt-in tests of the arbitrary-precision ISEM route (BigFloat spins and frequencies) against the Black
Hole Perturbation Toolkit Teukolsky paclet (MST method). They are not part of `test/runtests.jl`.

## References

`references/bhpt_<prec>_<group>.txt` hold the paclet output of `bhpt_acceptance_reference.wls`
(`wolframscript -file bhpt_acceptance_reference.wls <workingprecision> <group>`, Teukolsky paclet
1.1.1), one record per line, `tag|key|value`, with Mathematica's `InputForm` digits and precision
marks. For every case the normalization-independent quantities are stored: the separation constant,
the renormalized angular momentum, `Incidence/Transmission`, `Reflection/Transmission`,
`R(rprime)/Transmission` and `dR(rprime)/Transmission` of the In and Up solutions at
`rprime = 3.1126842514806725` (the binary double).

Each group exists at several working precisions (60, 120 and 200 digits; 400 for `omega = 10`, where
the paclet needs at least 200). A case is tested against the first pair of `REFERENCE_PAIRS` in
`runtests.jl` that holds it: 200/400 for `omega = 10`, 120/200 for the other real and complex cases,
120/200 of `bhpt_<prec>_axis.txt` (omega = +-delta -+ i sigma, delta from 1e-6 to 1e-30, sigma =
1/10, 1/2, 1 and 3, both axes, and the exact point i/10; the paclet's one-hour runs at the other exact
axis points -i sigma and i sigma timed out, see `bhpt_axis_report.txt`; the records at the exact point
-i/10 are principal-value diagnostics whose Up radial values are on another branch of the normalization
and are not tested) for the axis, and 60/120 for
the few axis cases whose 200-digit run did not finish. The digits on which the two agree are the certified
accuracy of the reference; the tests require agreement with the higher-precision run to
`min(N log10(2) - 6, certified digits - 1)` at `N` working bits. A reference without a counterpart is
reported but not tested. The 120/200 pairs certify the amplitudes to 85-99 digits and the radial values
to 54-60; the files were produced on another machine by the job described in
`references/codex_task_specification.txt` with `bhpt_acceptance_reference_batch.wls`.

`references/qnm_reference.txt` (report in `qnm_reference_report.txt`, generator `qnm_reference_worker.wls`)
holds 22 quasinormal frequencies of the s = -2, (2, 2) modes at a = 2/5 and 9/10 (n = 0..9 and the
second n = 8 root "8_1"), each obtained twice independently of this package, by the MST
incidence-amplitude zero of the QNM paclet at 120 digits and by a Leaver continued fraction at 100
digits (the two agree to 80-89 digits), with the frequency derivative of the Teukolsky incidence
amplitude by central differences at 120 digits (two step sizes agree to 38-45 digits).
`qnm_runtests.jl` polishes the reference's binary64 seed with `QNM.qnm_arbitrary_precision` and
compares the frequency and the derivative with the same pass rule.

`nia_jump_runtests.jl` checks the jump across the negative imaginary axis: the two one-sided limits of
the route at `omega = +-1e-30 - i sigma` against the analytic MST branch strength of the Float64 module,
`B_inc(-delta) = B_inc(+delta) + i q B_ref(+delta)` with the reflection continuous (no paclet data
involved; the strength is certified to about 1e-13, which sets the pass level).

## Running

```sh
GSN_AP_BITS=113,256 julia --project=<env> test/arbitrary_precision/runtests.jl
GSN_AP_BITS=53,113,192,256 GSN_AP_CASES=all julia --project=<env> test/arbitrary_precision/runtests.jl
GSN_AP_BITS=113,256 GSN_AP_QNM_CASES=all julia --project=<env> test/arbitrary_precision/qnm_runtests.jl
```

`GSN_AP_CASES` is `quick` (default), `real`, `complex`, `axis` or `all`; `GSN_SOURCE=<checkout>` loads the
package by `include` as in `test/extended`. The environment must provide a
`SpinWeightedSpheroidalHarmonics` that accepts BigFloat at complex `c` (the `ArbitraryPrecision`
branch, v1.4.1 or later); the registered 1.4.0 rejects it.

Frequencies with negative real part are not in the references on purpose: the principal-branch MST
formulas are not on the physical sheet there, and the package obtains that half plane from the
symmetry `(m, omega) -> (-m, -conj(omega))` with complex conjugation.
