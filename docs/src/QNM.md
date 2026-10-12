# Quasinormal modes

The QNM interface combines a Leaver continued-fraction root with an independent
GSN-ISEM radial evaluation. Units are $M=1$.

## Compact interface

```julia
using GeneralizedSasakiNakamura

ordinary_mode = qnm(0.68, -2, 2, 2, 0)
mirror_mode = qnm(0.68, -2, 2, 2, 0, mirror)
```

The five numerical inputs are `(a, s, l, m, n)`. `ordinary` is the default
branch. `mirror` polishes the negative-real root and evaluates its amplitudes
and frequency derivative. The radial solver may use the exact conjugation
symmetry; the QNM result is not copied from the ordinary branch.

A successful compact result exposes:

```julia
ordinary_mode.omega
ordinary_mode.lambda
ordinary_mode.incidence_amplitude
ordinary_mode.reflection_amplitude
ordinary_mode.incidence_derivative
ordinary_mode.excitation_factor
```

For a QNM frequency $\omega_q$, define

```math
\alpha^{{\rm GSN}}_{\ell m\omega_q}
=\left.\frac{\partial A^{{\rm GSN},{\rm inc}}_{\ell m\omega}}
{\partial\omega}\right|_{\omega=\omega_q}.
```

The GSN excitation factor returned by `qnm` is

```math
B^{{\rm GSN}}_{\ell m\omega_q}
=\frac{A^{{\rm GSN},{\rm ref}}_{\ell m\omega_q}}
{2\omega_q\alpha^{{\rm GSN}}_{\ell m\omega_q}}.
```

The derivative is obtained from Richardson-extrapolated frequency stencils.
The incidence zero and derivative convergence are checked before acceptance.
Wronskian checks apply only to routes that compute a Wronskian.
Richardson convergence does not measure systematic radial-amplitude error;
passing these checks does not guarantee a total `1e-8` error in the derivative
or excitation factor.

## Root-only interface

Use `qnm_frequency` when amplitudes are not needed:

```julia
root = qnm_frequency(0.68, -2, 2, 2, 0)                  # (a, s, l, m, n), as for qnm
mirror_root = qnm_frequency(0.68, -2, 2, 2, 0, mirror)

sequence = qnm_sequence([0.0, 0.3, 0.68], -2, 2, 2, 0)   # (spins, s, l, m, n)
```

Both functions take the same five numbers as [`qnm`](@ref) and the keywords
`convention` (how `n` is counted, below) and `segment` (which part of a
sequence that leaves and re-enters the negative imaginary axis, below). The
older form `qnm_frequency(QNMMode(s, l, m, n), a)` still works.

### Overtone labels: `convention`

Overtones are numbered by damping, but two numberings are in use. They differ
only for `s = -2`, `l = 2`, `|m| = 2` co-rotating modes (`a*m > 0` on the
ordinary branch), where the Schwarzschild algebraically special frequency
`-2i` (Schwarzschild `n = 8`) splits at `a ≠ 0` into a pair of roots, which
Cook calls `8₀` and `8₁`:

| root | `convention=:overtone` (default; Python `qnm`, Cook & Zalutskiy) | `convention=:complete_spectrum` (Forteza & Mourier, 2107.11829) |
| :--- | :--- | :--- |
| `n = 0, …, 7` | `n` | `n` |
| `8₀` (less damped of the pair) | `n = 8` | `n = 8` |
| `8₁` (more damped of the pair) | not numbered | `n = 9` |
| next roots | `n = 9, 10, …` | `n = 10, 11, …` |

`:overtone` follows each Schwarzschild root to finite spin, so `8₁`, which has
no Schwarzschild partner of its own, is skipped and every higher overtone keeps
its Schwarzschild number. `:complete_spectrum` numbers every root in order of
damping, so from `8₁` on its labels are one higher. For all other `(s, l, m)`
the two columns are identical. A third convention, `:leaver`, is only an index
of the continued-fraction inversion: it does not label a mode by itself and
requires an explicit `guess` (see below).

The default `convention=:overtone` treats `n` as a physical overtone label. The
labels follow the Python [`qnm`](https://github.com/duetosymmetry/qnm) package
(Stein 2019) and the Cook & Zalutskiy catalogue: overtone `n` is the Kerr
continuation of the `n`-th Schwarzschild root, the same mode that
[`qnm_sequence`](@ref) follows from `a = 0`. For `s = -2`, `l = 2`, `|m| = 2`
co-rotating modes the Schwarzschild algebraically special frequency `-2i` splits
at `a > 0` into a pair, which Cook labels `8₀` and `8₁`. Overtone `n = 8`
returns `8₀` (the less damped one); `8₁` is available as `n = 9` with
`convention=:complete_spectrum`, which numbers the pair consecutively:

```julia
qnm_frequency(0.68, -2, 2, 2, 8)                                    # 8₀
qnm_frequency(0.68, -2, 2, 2, 9; convention=:complete_spectrum)     # 8₁
qnm_frequency(0.68, -2, 2, 2, 9)                                    # overtone 9 = complete-spectrum 10
```

Overtones `n ≥ 9` therefore carry the same labels as the Python `qnm` package
and the Cook catalogue. For the co-rotating `(2, 2)` modes the correspondence
with tables that number the pair consecutively is

| this package and Python `qnm` / Cook | Forteza & Mourier (2107.11829), `convention=:complete_spectrum` |
| :--- | :--- |
| `n = 8` (`8₀`) | `n = 8` |
| `8₁` | `n = 9` |
| `n ≥ 9` | `n + 1` |

Counter-rotating modes have no such pair, and their labels coincide in these
conventions. The mirror branch is counter-rotating when `a*m > 0`; for signed
inputs with `a*m < 0` it is co-rotating and inherits the split pair and label
conversion from its ordinary `(a, -m)` partner.

### Sequences that end on the negative imaginary axis

Some sequences reach the negative imaginary axis (NIA) at finite spin. On the
axis the Leaver continued fraction has no minimal solution, so a mode there is a
polynomial solution of the confluent Heun equation (Cook & Zalutskiy, PRD 94,
104074). For `s = -2`, `l = 2`, `m = 0` and `9 ≤ n ≤ 26`, the sequence
attached to `a = 0` (Cook's `n₀`) ends on the axis; after a spin interval
without a mode, a second segment `n₁` leaves the axis and continues to larger
spin, touching the axis again at points that are not modes. The two segments
are separate sequences, and values are never interpolated across the gap.

A *segment* is one continuous piece of such a sequence in spin. Segment `0` is
the piece attached to the Schwarzschild root at `a = 0` (Cook's `n₀`); segment
`1` is the piece that leaves the axis again at larger spin (`n₁`). For most
modes only segment `0` exists and nothing needs to be set. `segment` selects the
piece:

| `segment` | returns |
| :--- | :--- |
| `0` (default) | `n₀`; fails with `:negative_imaginary_axis_endpoint` beyond its end |
| `1` | `n₁` (currently `m = 0` only); fails with `:negative_imaginary_axis_gap` before it starts |
| `:auto` | whichever piece exists at `a` |

```julia
qnm_frequency(0.35, -2, 2, 0, 10)                  # 10₀
qnm_frequency(0.5, -2, 2, 0, 10)                   # fails: 10₀ ended at a ≈ 0.391144
qnm_frequency(0.5, -2, 2, 0, 10; segment=1)        # 10₁
qnm_frequency(0.5, -2, 2, 0, 10; segment=:auto)    # first segment that exists at a
qnm(0.5, -2, 2, 0, 10; segment=1)                  # amplitudes of 10₁
```

When segment 0 reaches the axis the result fails with
`stop_reason=:negative_imaginary_axis_endpoint`. Its provenance records
`endpoint_spin` and `endpoint_frequency` (extrapolated from the accepted path to
`Re ω = 0`), the bracket `endpoint_spin_bracket`, the accepted `root_path`, and
the original solver stop reason `continuation_stop_reason`. With `segment=1`
the solver scans for roots leaving the axis at spins above the endpoint,
continues each candidate back down to the axis, and attributes the first event
whose axis frequency lies between the endpoints of overtones `n` and `n+1`.
The successful result records `reemergence` (spin and frequency where `n₁`
leaves the axis), `segment_endpoints`, and `branch_label = :kerr_overtone_n10_1`.

Before declaring a later segment ended on the NIA, continuation attempts a
restart in the same direction at `Re ω = ±1e-3`, using the last accepted
imaginary frequency and the existing root and branch checks. Successful
crossings are counted in `provenance.continuation_axis_touches`;
`nia_restart_rows` retains restart attempts and the original failures.
These records survive automatic precision fallback. Segment 0 keeps its
original continuation path.

A continuation step is not accepted when the interpolant of `Re ω` or of
`Re ω²` through the last accepted points vanishes inside the step: the path
reaches the axis there, and a root found beyond it belongs to another
sequence. Close to `a = 1` the touches of `n₁` accumulate. If a re-emerged
segment stalls on the axis after the restarts, the result fails with
`:negative_imaginary_axis_tangent_accumulation` rather than claiming an end of
the sequence. For `s = -2`, `l = 2`, `m = 0` and `n ≥ 17` this
happens on the high-spin part of segment 1: for `n = 26` from about `a ≈ 0.81`,
and for `n = 21`–`25` from about `a ≈ 0.91`–`0.97`. There the call fails with
`:negative_imaginary_axis_tangent_accumulation` and returns no frequency.
A Float64 root of segment 1 within two derivative steps of the
axis is accepted only when a BigFloat re-polish at the same spin reproduces it
to the stability gate; the check is recorded in
`provenance.near_axis_precision_check`. A BigFloat (at least 128-bit) final root
of segment 1 in that region must also pass a position certificate: the
continued-fraction derivative from three stencils (step `h`, `h/2` and an
imaginary step) must agree, and the Newton correction must not exceed the
stability gate (`provenance.near_axis_position_check`). If the 128-bit
certificate fails on its derivative, the root is solved again at the same spin
in 192-bit arithmetic and must pass the certificate there
(`provenance.near_axis_position_recovery`); the continuation itself is never
repeated in higher precision.

| `stop_reason` | meaning |
| :--- | :--- |
| `:negative_imaginary_axis_endpoint` | the requested segment ended on the NIA before `a` |
| `:negative_imaginary_axis_gap` | `a` lies between the end of `n₀` and the start of `n₁` |
| `:no_further_segment` | `segment=1` was requested where segment 0 exists |
| `:negative_imaginary_axis_tangent_accumulation` | segment 1 stalled on the axis after its birth (touches accumulate toward `a = 1`) |
| `:near_axis_precision_unverified` | a segment-1 Float64 root next to the axis was not reproduced by the BigFloat re-polish |
| `:near_axis_position_derivative_unverified`, `:near_axis_position_cf_unverified`, `:near_axis_position_stencil_unverified` | the BigFloat position certificate could not be established |
| `:near_axis_position_recovery_unverified` | the 192-bit same-spin re-solve did not pass the certificate |
| `:near_axis_position_branch_guard` | the certified root moved off the continuation branch |
| `:segment_attribution_unavailable` | `m ≠ 0`: re-emergence events are listed in `provenance.reemergence_events`, but their labels are not assigned |
| `:segment_unsupported_at_extremal` | `segment ≥ 1` at `abs(a) = 1` |
| `:initial_guess_off_branch` | a Newton start on or beyond the NIA for its branch |

For `m ≠ 0` sequences end near, not on, the axis, and Cook notes that the
re-emerged labels are heuristic, so segment 1 is not assigned. The mirror of
segment `k` is `-conj` of the ordinary segment `k`. `segment` is not available
for the `8₀`/`8₁` pair and is limited to `0`, `1`, and `:auto`.

The continued-fraction inversion convention is available separately and
requires both an inversion index and a root guess:

```julia
leaver_root = qnm_frequency(
    0.68, -2, 2, 2, 0;
    convention=:leaver,
    inversion_index=root.inversion_index,
    guess=root.omega,
)
```

`convention=:complete_spectrum` numbers the complete spectrum consecutively in
the algebraically special neighborhood (the pair as 8 and 9, then `n ≥ 10`),
as in Forteza & Mourier.

## Detailed radial solutions

The compact call avoids constructing callable radial functions. Request the
full bridge only when it is needed:

```julia
detailed = qnm(0.68, -2, 2, 2, 0; detailed=true)

r = 10.0
rstar = rstar_from_r(0.68, r)
detailed.X(rstar)
detailed.Y(r)
detailed.R(r)
```

`X` uses the tortoise coordinate $r_*$; `Y` and `R` use Boyer-Lindquist $r$.
The detailed normalization check compares GSN and Teukolsky conventions;
passing it is not an independent accuracy bound for either solution.

`detailed.teukolsky.excitation_factor` uses the convention

```math
B^{\mathrm{T}}_{\ell m\omega_q}
=\frac{(B^{\mathrm{ref}}/B^{\mathrm{trans}})_{\ell m\omega_q}}
{2i\omega_q\left.\partial_\omega
(B^{\mathrm{inc}}/B^{\mathrm{trans}})_{\ell m\omega}\right|_{\omega_q}}.
```

A reference defined with `2 omega` instead of `2 i omega` must therefore be
divided by `i` before comparison. This phase convention is distinct from a
radial normalization error.

## Result types

- `QNMResult`: all required root, radial, normalization, and derivative gates
  pass.
- `QNMEstimate`: finite diagnostic result whose limiting gate remains recorded.
- `QNMFailure`: the first failed stage and reason are returned without
  promoting nonfinite amplitudes.
- `QNMEndpointResult`: an exact extremal branch endpoint that is not treated as
  an ordinary isolated simple pole.

At the Schwarzschild algebraically-special frequency and at a synchronous
extremal accumulation endpoint, the ordinary simple-zero excitation formula
is not applied.

## Exact extremal spin

For `a = 1` or `a = -1` (exact equality) each mode is first followed to the
near-extremal spins `a = sqrt(1 - κ²)`, `κ = 0.04, …, 0.0025`. A branch whose
distance from `m a/2` halves with `κ` is a zero-damping branch and returns the
synchronous endpoint `ω = m a/2` as a `QNMEndpointResult`. Every other branch is
a damped mode: its frequency is the root of the exact-extremal radial
recurrence of Richartz (Phys. Rev. D 93, 064062 (2016)), seeded by the
near-extremal limit. The root is accepted only when it is stable under doubling
of the recurrence depth and lies on the branch of the near-extremal ladder.
The `P_3` coefficient of that paper is used with `-2im(1 + 6iω)`; the printed
`-2im(1 - 6iω)` does not satisfy the Teukolsky equation for `m ≠ 0`.

The radial amplitudes at `a = ±1` come from integrating out of an irregular
singular horizon. For co-rotating damped modes this is ill-conditioned in
`Float64`; `qnm` then keeps the exact frequency but reports the failed incidence
gate instead of an excitation factor.

## Lower-level functions

```@docs
QNMMode
QNMBranch
QNMResult
QNMEstimate
QNMFailure
QNMEndpointResult
QNMPairResult
LeaverResult
ISEMValidationResult
ExcitationFactorResult
angular_A_to_lambda
lambda_to_angular_A
qnm
qnm_pair
qnm_frequency
qnm_sequence
validate_qnm_with_isem
qnm_excitation_factor
```
