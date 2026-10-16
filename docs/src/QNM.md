# Quasinormal modes

A QNM frequency is computed as a root of Leaver's continued fraction (a
*Leaver root*); its scattering amplitudes and excitation factor are then
evaluated at that frequency by solving the GSN equation (see below). Units are
$M=1$.

## Compact interface

```julia
using GeneralizedSasakiNakamura

ordinary_mode = qnm(0.68, -2, 2, 2, 0)
mirror_mode = qnm(0.68, -2, 2, 2, 0, mirror)
```

The five numerical inputs are `(a, s, l, m, n)`. `ordinary` (the default) is the
mode with positive real frequency and `mirror` the one with negative real
frequency. By default the mirror frequency is the exact reflection `-conj(ω)` of
the ordinary `(s, l, -m, n)` root at the same spin; pass `mirror_guess` to
polish it directly instead. Its amplitudes and frequency derivative are
evaluated at the mirror frequency.

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

For `s = -2, …, 2`, `|s| ≤ l ≤ 15` and `n ≤ 10` the derivative comes from a
single radial solve: the GSN equation is integrated together with its frequency
derivative (forward-mode dual numbers) along a complex contour on which the
ingoing and outgoing solutions at infinity stay balanced. The solve is done at
two resolutions. It is accepted when both agree to `1e-8`, the incidence
amplitude vanishes to `1e-8` relative to the reflection amplitude, and the
Wronskian is constant to `1e-8`. Otherwise, and for all other modes, the
derivative is Richardson-extrapolated from frequency stencils. A result whose
checks do not reach `1e-8` is returned as a `QNMEstimate` that records the
limiting check.

## Root-only interface

Use `qnm_frequency` when amplitudes are not needed:

```julia
root = qnm_frequency(0.68, -2, 2, 2, 0)                  # (a, s, l, m, n), as for qnm
mirror_root = qnm_frequency(0.68, -2, 2, 2, 0, mirror)

sequence = qnm_sequence([0.0, 0.3, 0.68], -2, 2, 2, 0)   # (spins, s, l, m, n)
```

### Packaged roots

For `s = -2, …, 2`, `|s| ≤ l ≤ 15` and `n ≤ nmax(l)` the package ships the
roots on 60 spins `a = 1 - u²` (`u` from `1` to `0.01`) and at `a = 1`:

| `l` | 0–5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 | 13 | 14 | 15 |
| :-- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `nmax(l)` | 10 | 9 | 8 | 8 | 7 | 6 | 5 | 4 | 4 | 3 | 2 |

The grid runs from `a = 0` to `a = 0.9999`. The tabulated roots are starting
points only: the returned root is polished with the Leaver continued fraction
at the requested spin. For `0.999 < a < 1` a damped mode of the table is
instead evaluated from the near-extremal expansion
`ω = c₀ + c₁x + … + c₈x⁸` in `x = 1 - a²` (ordinary branch, `Float64` spin,
default keywords). Its coefficients are solved order by order from the
Teukolsky equation at `a = 1`, starting from the packaged `a = 1` root, on first
use of the mode (0.03–0.3 s) and then cached. The expansion is used where the
term it omits, of order `exp(-2π Re(2ω - m)/√x)`, is below `1e-16`.

Both functions take the same numbers as [`qnm`](@ref) (`qnm_sequence` takes a list of
spins in place of `a`) and the keywords
`convention` (how `n` is counted, below) and `segment` (which part of a
sequence that leaves and re-enters the negative imaginary axis, below). The
form `qnm_frequency(QNMMode(s, l, m, n), a)` is equivalent.

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

At `a = 0` the pair has merged into `-2i`. Both `:complete_spectrum` `n = 8` and
`n = 9` return that frequency, as `:overtone` `n = 8` does, and
`:complete_spectrum` `n ≥ 10` is the overtone `n - 1`. There `qnm` returns a
`QNMFailure`: the incidence amplitude has no simple zero at `-2i`, so the
excitation factor is not defined.

The `:overtone` labels follow the Python [`qnm`](https://github.com/duetosymmetry/qnm)
package (Stein 2019) and the Cook & Zalutskiy catalogue: overtone `n` is the Kerr
continuation of the `n`-th Schwarzschild root, the same mode that
[`qnm_sequence`](@ref) follows from `a = 0`. For example, at `a = 0.68`:

```julia
qnm_frequency(0.68, -2, 2, 2, 8)                                    # 8₀
qnm_frequency(0.68, -2, 2, 2, 9; convention=:complete_spectrum)     # 8₁
qnm_frequency(0.68, -2, 2, 2, 9)                                    # overtone 9 = complete-spectrum 10
```

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
leaves the axis), `segment_endpoints`, and a `branch_label` such as `:kerr_overtone_n10_1`.

Near the axis the continuation restarts a stalled step, and segment-1 roots
next to the axis are confirmed in BigFloat arithmetic; the corresponding checks
are recorded in the result's `provenance`. Close to `a = 1` the touches of `n₁`
accumulate, and a re-emerged segment that stalls on the axis fails with
`:negative_imaginary_axis_tangent_accumulation` instead of returning a
frequency. For `s = -2`, `l = 2`, `m = 0` this happens on the high-spin part of
segment 1 of the higher overtones, for example from about `a ≈ 0.81` for
`n = 26` and from about `a ≈ 0.91`–`0.97` for `n = 21`–`25`.

| `stop_reason` | meaning |
| :--- | :--- |
| `:negative_imaginary_axis_endpoint` | the requested segment ended on the NIA before `a` |
| `:negative_imaginary_axis_gap` | `a` lies between the end of `n₀` and the start of `n₁` |
| `:no_further_segment` | `segment=1` was requested where segment 0 exists |
| `:negative_imaginary_axis_tangent_accumulation` | segment 1 stalled on the axis after its birth (touches accumulate toward `a = 1`) |
| `:near_axis_precision_unverified`, `:near_axis_position_*` | a segment-1 root next to the axis could not be confirmed in BigFloat arithmetic |
| `:segment_attribution_unavailable` | `m ≠ 0`: re-emergence events are listed in `provenance.reemergence_events`, but their labels are not assigned |
| `:segment_unsupported_at_extremal` | `segment ≥ 1` at `abs(a) = 1` |
| `:initial_guess_off_branch` | a Newton start on or beyond the NIA for its branch |

For `m ≠ 0` sequences end near, not on, the axis, and Cook notes that the
re-emerged labels are heuristic, so segment 1 is not assigned. The mirror of
segment `k` is `-conj` of the ordinary segment `k`. `segment` is not available
for the `8₀`/`8₁` pair and is limited to `0`, `1`, and `:auto`.

### Leaver inversion index

With `convention=:leaver`, `n` (or an explicit `inversion_index`) is the
inversion index of the continued fraction, which does not identify a mode by
itself, so a `guess` is required:

```julia
leaver_root = qnm_frequency(
    0.68, -2, 2, 2, 0;
    convention=:leaver,
    inversion_index=root.inversion_index,
    guess=root.omega,
)
```

## Detailed radial solutions

The default call does not construct callable radial functions. Request them,
together with the excitation factor in the Teukolsky normalization, only when
they are needed:

```julia
detailed = qnm(0.68, -2, 2, 2, 0; detailed=true)

r = 10.0
rstar = rstar_from_r(0.68, r)
detailed.X(rstar)
detailed.Y(r)
detailed.R(r)
```

`X` uses the tortoise coordinate $r_*$; `Y` and `R` use Boyer-Lindquist $r$.

`detailed.teukolsky.excitation_factor` uses the convention

```math
B^{\mathrm{T}}_{\ell m\omega_q}
=\frac{(B^{\mathrm{ref}}/B^{\mathrm{trans}})_{\ell m\omega_q}}
{2i\omega_q\left.\partial_\omega
(B^{\mathrm{inc}}/B^{\mathrm{trans}})_{\ell m\omega}\right|_{\omega_q}}.
```

A reference defined with `2 omega` instead of `2 i omega` must therefore be
divided by `i` before comparison.

## Result types

- `QNMResult`: the root, radial, normalization and derivative checks all pass.
- `QNMEstimate`: a finite result with one check short of its target; the check
  is named in `limiting_gate`.
- `QNMFailure`: the stage that failed and the reason.
- `QNMEndpointResult`: an exact-extremal endpoint (`a = ±1`, `ω = m a/2`),
  which is not a simple pole.

At the Schwarzschild algebraically special frequency and at the synchronous
exact-extremal endpoint the excitation factor is not defined.

## Exact extremal spin

For `a = 1` or `a = -1` (exact equality) a mode in the packaged table starts
from its tabulated `a = 1` root, which must be reproduced by the recurrence
below to `1e-10`. Any other mode is first followed to the near-extremal spins `a = sqrt(1 - κ²)`, `κ = 0.04, …, 0.0025`. A branch whose
distance from `m a/2` halves with `κ` is a zero-damping branch and returns the
synchronous endpoint `ω = m a/2` as a `QNMEndpointResult` (displayed as
`ExtremalQuasiNormalMode` with `family = ZDM`, a zero-damping mode). Every other
branch is a damped mode (`family = DM`): its frequency is the root of the exact-extremal radial
recurrence of Richartz (Phys. Rev. D 93, 064062 (2016)), seeded by the
near-extremal limit. The root is accepted only when it is stable under doubling
of the recurrence depth and lies on the branch of the near-extremal ladder.
The `P_3` coefficient of that paper is used with `-2im(1 + 6iω)`; the printed
`-2im(1 - 6iω)` does not satisfy the Teukolsky equation for `m ≠ 0`.

The radial amplitudes at `a = ±1` are integrated along a complex path out of
the irregular singular horizon. If the incidence amplitude does not vanish to
its target there, `qnm` keeps the exact frequency and reports the failed check
instead of an excitation factor. At exact extremality the mirror root is the reflection `-conj(ω)` of the ordinary
`(s, l, -m, n)` root at the same spin.

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
