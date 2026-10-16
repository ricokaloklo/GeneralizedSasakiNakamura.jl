# Quasinormal modes

`qnm` finds the QNM frequency as a root of Leaver's continued fraction, then
evaluates the GSN scattering amplitudes and the excitation factor at that
frequency with a radial solver independent of the root search. A few
frequencies come from another equation: for some damped modes at
$0.999 < |a| < 1$ from a near-extremal series, at $|a| = 1$ from an
exact-extremal recurrence or, for a zero-damping mode, as $\omega = ma/2$ (see
below), and for the pair `8₀`/`8₁` that leaves $-2i$ (see Overtone labels) as a
zero of the GSN incidence amplitude at $0 < |a| < 0.015$ (for `8₁` at
$0 < |a| < 0.05$).

Units are $G = c = M = 1$. Perturbations are taken proportional to
$e^{-i\omega t + im\phi}$, so a damped QNM has $\operatorname{Im}\omega < 0$.
The spin satisfies $-1 \le a \le 1$; a negative `a` is evaluated through the
symmetry $(a, m) \to (-a, -m)$.

## Compact interface

```julia
using GeneralizedSasakiNakamura

ordinary_mode = qnm(0.68, -2, 2, 2, 0)
mirror_mode = qnm(0.68, -2, 2, 2, 0, mirror)
```

The five numerical inputs are `(a, s, l, m, n)`. `ordinary` (the default) is the
mode with positive real frequency and `mirror` the one with negative real
frequency. By default the mirror frequency is the exact reflection `-conj(ω)` of
the ordinary `(s, l, -m, n)` root at the same spin. Only when that root is not
accepted, and did not stop at the negative imaginary axis, is the mirror root
polished instead, starting from the reflection. Pass `mirror_guess` to polish
it from your own starting value. Its amplitudes and frequency derivative are
evaluated at the mirror frequency, not copied from the ordinary branch.

`qnm` also takes the keywords `convention`, `segment` and `inversion_index` of
`qnm_frequency` (see Root-only interface) and `detailed`; where `qnm_frequency`
takes `guess`, `qnm` takes `primary_guess` for the ordinary root and
`mirror_guess` for the mirror root. [`qnm_pair`](@ref)`(a, s, l, m, n)` computes
both branches in one [`QNMPairResult`](@ref).

A successful compact result exposes:

```julia
ordinary_mode.omega                  # ω (the same value as .frequency)
ordinary_mode.lambda                 # Teukolsky separation constant λ
ordinary_mode.incidence_amplitude    # A^{GSN,inc}
ordinary_mode.reflection_amplitude   # A^{GSN,ref}
ordinary_mode.incidence_derivative   # α^{GSN}
ordinary_mode.excitation_factor      # B^{GSN}
```

Check `status` before using these values (see Result types below). `lambda` is
the separation constant of the radial Teukolsky equation,
$\lambda = A + c^2 - 2mc$ with $c = a\omega$ and $A$ the angular eigenvalue in
Leaver's convention (`angular_A`). The amplitudes belong to the `IN` solution of
the GSN equation with unit transmission amplitude:

```math
X^{\rm in} \to
\begin{cases}
e^{-ipr_*}, & r_* \to -\infty, \\
A^{{\rm GSN},{\rm inc}}_{\ell m\omega} e^{-i\omega r_*}
+ A^{{\rm GSN},{\rm ref}}_{\ell m\omega} e^{+i\omega r_*}, & r_* \to \infty,
\end{cases}
```

with $p = \omega - m\Omega_H$ and $\Omega_H = a/(2r_+)$. At a QNM frequency
$\omega_q$ the incidence amplitude vanishes. Define

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

For `s = -2, …, 2`, `|s| ≤ l ≤ 15`, `n ≤ 10` and a damped root with
$\operatorname{Re}\omega \neq 0$, the amplitudes and $\alpha$ come from
integrating the GSN equation together with its frequency derivative along a
complex contour at the root frequency only, at two resolutions (the direct
route). It is accepted when the two resolutions agree to `1e-8`,
$|A^{\rm inc}|/\max(1, |A^{\rm ref}|) \le 10^{-8}$ and the Wronskian is
constant to `1e-8`. If any check fails, and for all other modes, $\alpha$ is
Richardson-extrapolated from central differences in real frequency of ISEM
amplitudes from `GSN_radial` (the Richardson route; for `a = ±1` see Exact
extremal spin). That route is accepted when the same incidence and Wronskian
checks pass and the step drift (the relative change between successive
Richardson estimates) is at most `1e-7`. When the route that is used misses a
target, or the root itself has status `:estimated`, `qnm` returns a
[`QNMEstimate`](@ref) in the cases listed in its docstring and a `QNMFailure`
otherwise.

`excitation_options` with `backend=:richardson`, an explicit stencil `h`, a
`direction_tolerance`, or a custom `radial_solver` or `teukolsky_solver` always
selects the Richardson route. `backend=:variational` forces the direct route
and returns its result without fallback. Outside the scope above, or with those
options, it raises an `ArgumentError`: `qnm` returns it as a `QNMFailure` with
`stage = :radial_evaluation` and the text in `error_message` (for a damped root
at `a = ±1`, `qnm` raises it instead), and `qnm_excitation_factor` raises it.
`result.excitation.metadata` has the entry `backend = :variational` only for a
direct-route result, and `variational_fallback_reason` after a fallback.

## Root-only interface

Use `qnm_frequency` when amplitudes are not needed:

```julia
root = qnm_frequency(0.68, -2, 2, 2, 0)                  # (a, s, l, m, n), as for qnm
mirror_root = qnm_frequency(0.68, -2, 2, 2, 0, mirror)

sequence = qnm_sequence([0.0, 0.3, 0.68], -2, 2, 2, 0)   # (spins, s, l, m, n)
```

Both functions take the same numbers as [`qnm`](@ref) (`qnm_sequence` takes a
list of spins in place of `a`) and the keywords `convention` (how `n` is
counted, below), `segment` (which piece of a sequence that ends on the negative
imaginary axis and re-emerges from it, below), `inversion_index` and `guess`.
For the ordinary branch, `qnm_frequency(QNMMode(s, l, m, n), a)` is equivalent;
a `QNMMode` with `:negative_real` needs a `guess`.

`qnm_frequency` returns a `LeaverResult` with the root only: `omega`
(`ComplexF64`, or `Complex{BigFloat}` when the root needed extended precision),
`lambda`, `angular_A`, `status` (`:accepted`, `:estimated` or `:failed`),
`stop_reason`, `overtone_index`, `inversion_index` and `provenance`.
`:estimated` marks a root that passed only relaxed checks, for example a
residual floor close to `|a| = 1`, an uncertain overtone label or an
unconfirmed root at `|a| = 1`; `stop_reason` gives the reason. A failed result can hold a
finite `omega`; it is a diagnostic value of the root search (for example the
last rejected candidate), not a QNM frequency. `qnm_sequence` returns a vector
of `LeaverResult`s; each accepted root seeds the next spin, and the vector
stops after the first result that is not `:accepted`, so it can be shorter than
the list of spins. The mirror branch follows the rule of `qnm`: without a
`guess`, `qnm_sequence` reflects the ordinary root at its first spin and
polishes each later spin from the previous mirror root.

### Overtone labels: `convention`

Overtones are numbered by damping, but two numberings are in use. They differ
only for `s = -2`, `l = 2`, `|m| = 2`. For co-rotating modes (`a*m > 0` on the
ordinary branch, `a*m < 0` on the mirror branch) the Schwarzschild
algebraically special frequency `-2i` (Schwarzschild `n = 8`) splits at
`a ≠ 0` into a pair of roots, which Cook calls `8₀` and `8₁`:

| root | `convention=:overtone` (default; Python `qnm`, Cook & Zalutskiy) | `convention=:complete_spectrum` (Forteza & Mourier, 2107.11829) |
| :--- | :--- | :--- |
| `n = 0, …, 7` | `n` | `n` |
| `8₀` (less damped of the pair) | `n = 8` | `n = 8` |
| `8₁` (more damped of the pair) | not numbered | `n = 9` |
| next roots | `n = 9, 10, …` | `n = 10, 11, …` |

The `:overtone` labels follow the Python [`qnm`](https://github.com/duetosymmetry/qnm)
package (Stein 2019) and the Cook & Zalutskiy catalogue: overtone `n` is the Kerr
continuation of the `n`-th Schwarzschild root, the same mode that
[`qnm_sequence`](@ref) follows from `a = 0`. `8₁` has no Schwarzschild root of
its own, so it is not numbered and every higher overtone keeps its
Schwarzschild number. `:complete_spectrum` numbers every root in order of
damping, so from `8₁` on its labels are one higher. For example, at `a = 0.68`:

```julia
qnm_frequency(0.68, -2, 2, 2, 8)                                    # 8₀
qnm_frequency(0.68, -2, 2, 2, 9; convention=:complete_spectrum)     # 8₁
qnm_frequency(0.68, -2, 2, 2, 9)                                    # overtone 9 = complete-spectrum 10
```

For counter-rotating `s = -2`, `l = 2`, `|m| = 2` modes `:overtone` keeps the
Schwarzschild number for every `n`, while `:complete_spectrum` still treats
`n = 8` and `n = 9` as the pair that leaves `-2i` and maps `n ≥ 10` to overtone
`n - 1`, as it also does at `a = 0`. For all other `(s, l, m)` the two
conventions give the same labels.

At `a = 0` the pair has merged into `-2i`: `:overtone` and `:complete_spectrum`
`n = 8` (for `|s| = 2`, `l = 2`, any `m`) and `:complete_spectrum` `n = 9` (for
`s = -2`, `l = 2`, `|m| = 2`) return exactly `-2i`. `qnm_frequency` reports it
with `status = :accepted` and
`stop_reason = :exact_algebraically_special_endpoint`; `qnm` returns a
`QNMFailure` with `stage = :pole_skipping_endpoint` and
`stop_reason = :matsubara_qnm_pole_skipping_no_simple_gsn_residue`, because the
incidence amplitude has no simple zero at `-2i` and the excitation factor is not
defined.

The `8₀`/`8₁` numbering is implemented for `s = -2` only. For `s = +2` the label
`n` goes to the root search unchanged: the packaged table has no accepted
co-rotating `s = +2`, `l = 2`, `n = 8` root at `a > 0`, and the `s = +2`,
`l = 2`, `n = 10` roots start from the Schwarzschild root `0.0826 - 2.8154i`,
which is overtone 11 for `s = -2`. QNM frequencies do not depend on the sign of
`s`, so for `l = 2` and `n ≥ 8` use `s = -2`.

A third convention, `:leaver`, makes `n` (or an explicit `inversion_index`,
which overrides it) the index at which Leaver's continued fraction is inverted.
All inversions have the same roots, so this index does not identify a mode by
itself: `:leaver` requires a `guess` and returns `overtone_index = missing`.

```julia
leaver_root = qnm_frequency(
    0.68, -2, 2, 2, 0;
    convention=:leaver,
    inversion_index=root.inversion_index,
    guess=root.omega,
)
```

### Sequences that end on the negative imaginary axis

Some sequences reach the negative imaginary axis (NIA) at finite spin. On the
axis the Leaver continued fraction has no minimal solution, and a mode there is
a polynomial solution of the confluent Heun equation (Cook & Zalutskiy, PRD 94,
104074). For `s = -2`, `l = 2`, `m = 0` and `9 ≤ n ≤ 26`, segment `0` of the
sequence (Cook's `n₀`, attached to the Schwarzschild root) ends on the axis.
After a spin interval without a mode, segment `1` (`n₁`) leaves the axis and
continues to larger spin, touching it again at points that are not modes.
Values are never interpolated across the gap. Most modes have only segment `0`,
and nothing needs to be set. `segment` selects the piece:

| `segment` | returns |
| :--- | :--- |
| `0` (default) | `n₀`; fails with `:negative_imaginary_axis_endpoint` beyond its end |
| `1` | `n₁` (`m = 0` only); fails with `:negative_imaginary_axis_gap` before it starts |
| `:auto` | for `m = 0`, segment `0` where it exists and otherwise segment `1`; in the gap, the segment-0 failure `:negative_imaginary_axis_endpoint` with `provenance.segment_search = :negative_imaginary_axis_gap`. For `m ≠ 0`, segment `0` |

```julia
qnm_frequency(0.35, -2, 2, 0, 10)                  # 10₀
qnm_frequency(0.5, -2, 2, 0, 10)                   # fails: 10₀ ended at a ≈ 0.391144
qnm_frequency(0.5, -2, 2, 0, 10; segment=1)        # 10₁
qnm_frequency(0.5, -2, 2, 0, 10; segment=:auto)    # first segment that exists at a
qnm(0.5, -2, 2, 0, 10; segment=1)                  # amplitudes of 10₁
```

A segment-0 failure on the axis records `endpoint_spin` and
`endpoint_frequency` (extrapolated to `Re ω = 0`) in `provenance`. With
`segment=1` the solver scans spins above that endpoint for roots that leave the
axis, and keeps the first whose axis frequency lies between the segment-0
endpoints of overtones `n` and `n + 1`. The result records `reemergence` (where
`n₁` leaves the axis) and a `branch_label` such as `:kerr_overtone_n10_1`.
The touches of `n₁` accumulate toward `a = 1`; where the continuation stalls
on the axis, it stops with `:negative_imaginary_axis_tangent_accumulation`.
For `s = -2`, `l = 2`, `m = 0` this happens, for example, from `a ≈ 0.81` for
`n = 26` and from `a ≈ 0.91`–`0.97` for `n = 21`–`25`.

| `stop_reason` | meaning |
| :--- | :--- |
| `:negative_imaginary_axis_endpoint` | segment 0 reached the NIA (for `m ≠ 0`, stalled next to it) before `a` |
| `:negative_imaginary_axis_gap` | `a` lies between the end of `n₀` and the start of `n₁` |
| `:no_further_segment` | `segment=1` was requested where segment 0 exists |
| `:negative_imaginary_axis_tangent_accumulation` | the re-emerged segment 1 stalled on the axis (touches accumulate toward `a = 1`) |
| `:near_axis_precision_unverified`, `:near_axis_position_*` | a segment-1 root next to the axis could not be confirmed in BigFloat arithmetic |
| `:segment_attribution_unavailable` | `m ≠ 0`: re-emergence events are listed in `provenance.reemergence_events`, but their labels are not assigned |
| `:segment_unsupported_at_extremal` | `segment ≥ 1` at `abs(a) = 1` |
| `:initial_guess_off_branch` | the root search started from a predicted frequency with `Re ω` zero or of the wrong sign for its branch (a `guess` like that raises an `ArgumentError`) |

For `m ≠ 0` sequences end near, not on, the axis, and Cook notes that the
re-emerged labels are heuristic, so segment 1 is not assigned; segment 0 still
stops with `:negative_imaginary_axis_endpoint` next to the axis, for example
for `(s, l, m, n) = (-2, 2, -2, 15)` at `a = 0.7`. The mirror of segment `k` is
the reflection `-conj(ω)` of segment `k` of the ordinary `(s, l, -m, n)`
sequence. `segment` is limited to `0`, `1` and `:auto`. For the `8₀`/`8₁` pair,
and for `-2i` itself at `a = 0`, any `segment` other than `0`, `:auto`
included, raises an `ArgumentError`.

### Packaged roots

The package ships the ordinary roots of the default continuation for
`s = -2, …, 2`, `|s| ≤ l ≤ 15`, every `m` and `n ≤ nmax(l)` on 60 spins
`a = 1 - u²` (`u` from `1` to `0.01`, so `a` from `0` to `0.9999`):

| `l` | 0–5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 | 13 | 14 | 15 |
| :-- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `nmax(l)` | 10 | 9 | 8 | 8 | 7 | 6 | 5 | 4 | 4 | 3 | 2 |

The table holds 8374 modes. 116 labels in this range are absent, because the
tabulated root disagreed with the default continuation (high overtones near
algebraically special frequencies, for example `s = -2`, `l = 2`, `m = 2`,
`:complete_spectrum` `n = 9` and `10`), and 554 modes lack a root at some
spins, mostly above `a = 0.95`. A second table holds the `a = 1` roots of 5493
modes (5393 damped, 100 zero-damping).

The tabulated roots are only starting points, used on the ordinary branch (and
so for mirror reflections) for a `Float64` or integer spin with the default
`convention`, `inversion_index` and `segment` and no `guess`; the returned root
is polished with the continued fraction. A mode or spin missing from the table
is continued from `a = 0`, which is slower and can fail.

### Near-extremal spins

For `0.999 < |a| < 1` the ordinary root of a damped mode is evaluated from the
near-extremal expansion `ω = c₀ + c₁x + … + c₈x⁸` in `x = 1 - a²`. `c₀` is the
packaged `a = 1` root, and `c₁, …, c₈` are solved order by order from the
Teukolsky equation at `a = 1` on first use of the mode (0.03–0.3 s) and then
cached. This root is not polished with the continued fraction; its
`stop_reason` is `:accepted_near_extremal_expansion`. The expansion is used
only when all of these hold:

- the call has a `Float64` spin and only default keywords (`convention=:overtone`,
  no `guess`, `inversion_index` or `root_options`, `segment = 0`), so
  `qnm_sequence` uses it at most for its first spin;
- the mode has a damped entry in the `a = 1` table and an accepted root at the
  largest tabulated spin `a ≤ 0.999`, and the expansion reproduces that root to
  `1e-8`;
- the term the expansion omits, of order `exp(-2π Re(2c₀ - m)/√x)`, is below
  `1e-16` already at that tabulated spin (`a ≈ 0.998`), which needs
  `Re(2c₀ - m) ≥ 0.3606`; the expansion then applies at every spin above
  `0.999`.

Otherwise the root is solved with the continued fraction; above `|a| = 0.9995`
its starting point is marched from the tabulated root nearest below
`a = 0.9995` to the requested spin, with the packaged `a = 1` root as the far
anchor of each prediction, when both tables have the mode.

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
`X` is the `IN` solution of the GSN equation with unit GSN transmission
amplitude, and `R` the `IN` solution of the Teukolsky equation with unit
Teukolsky transmission amplitude ($R \to \Delta^{-s} e^{-ipr_*}$ at the
horizon). `Y` is the `IN` solution for `s = -2` and the `UP` solution for
`s = +2`, and is `missing` for other `s`.

`detailed=true` also adds checks: the GSN and Teukolsky normalizations of the
excitation factor must agree to `1e-10`, and on the Richardson route the
derivative is also taken along imaginary frequency and must agree with the real
direction to `1e-7`. Its status can therefore be stricter than that of the
compact call.

`detailed.teukolsky` holds `reflection_amplitude`, `incidence_derivative` and
`excitation_factor` in the Teukolsky normalization. With $B^{\mathrm{inc}}$,
$B^{\mathrm{ref}}$ and $B^{\mathrm{trans}}$ the Teukolsky incidence, reflection
and transmission amplitudes of the `IN` solution, the excitation factor is

```math
B^{\mathrm{T}}_{\ell m\omega_q}
=\frac{(B^{\mathrm{ref}}/B^{\mathrm{trans}})_{\ell m\omega_q}}
{2i\omega_q\left.\partial_\omega
(B^{\mathrm{inc}}/B^{\mathrm{trans}})_{\ell m\omega}\right|_{\omega_q}}.
```

To compare it with a reference defined with $2\omega_q$ in place of
$2i\omega_q$, divide the reference value by $i$.

## Result types

`qnm` returns one of four types:

| type | `status` | REPL display | returned when |
| :--- | :--- | :--- | :--- |
| `QNMResult` | `:accepted` | `QuasiNormalMode(`; at `a = ±1` `ExtremalQuasiNormalMode(` with `family = DM` | every check passes |
| `QNMEstimate` | `:estimated` | `QuasiNormalModeEstimate(`; at `a = ±1` `ExtremalQuasiNormalModeEstimate(` with `family = DM` | the values are finite, but a target is missed or the root itself is `:estimated` (below) |
| `QNMFailure` | `:failed` | `QNMFailure(` | a step fails; `stage` and `stop_reason` say which and why |
| `QNMEndpointResult` | `:accepted` | `ExtremalQuasiNormalMode(` with `family = ZDM` and `simple_pole = false` | a zero-damping mode at `a = ±1` |

`DM` stands for damped mode and `ZDM` for zero-damping mode. The cases that give
a [`QNMEstimate`](@ref) (named in `limiting_gate`), the stages of a
[`QNMFailure`](@ref) and the meaning of its `omega` are given in their
docstrings under Reference below; `target_relative_precision` is `1e-7` on the Richardson
route and `1e-8` on the direct route. A [`QNMEndpointResult`](@ref) has
`missing` amplitude fields and `X`, `Y`, `R`, because the synchronous frequency
`ω = m a/2` is a branch point, not a simple pole. The `ordinary` and `mirror`
fields of a [`QNMPairResult`](@ref) each hold one of these types.

## Exact extremal spin

For `a = 1` or `a = -1` (exact equality; `a = -1` is evaluated as `a = 1` with
`m → -m`) the root is found as follows:

1. With default keywords, a mode in the packaged `a = 1` table starts from its
   tabulated root. A zero-damping entry returns the synchronous endpoint
   `ω = m a/2` directly; a damped entry is re-solved with the exact-extremal
   recurrence below and kept when it reproduces the table to `1e-10`.
2. Otherwise, with default keywords, the roots at
   `|a| = 0.9, 0.94, 0.96, 0.98, 0.99` are extrapolated in `x = 1 - a²` to seed
   the recurrence. The result is kept when all these roots are accepted, the
   recurrence root lies within `1e-4` (relative) of the extrapolation and
   `Re(2ω - m a) > 0.02`.
3. Otherwise the mode is followed to the near-extremal spins
   `|a| = sqrt(1 - κ²)`, `κ = 0.04, 0.02, 0.01, 0.005, 0.0025`. An
   ordinary-branch root with `a*m > 0` whose distance from `m a/2` drops to at
   most `0.7` of its previous value over each of the last two steps (it halves
   with `κ`) is a zero-damping mode, and `qnm` returns the synchronous
   endpoint `ω = m a/2` as a `QNMEndpointResult`. A root whose distance keeps
   at least `0.8` of its previous value over both steps, and every other root
   (a mirror root that is polished rather than reflected included), is a damped
   mode (`family = DM`); other ratios leave the family ambiguous. For a damped
   or ambiguous branch the frequency is the root of the exact-extremal radial
   recurrence of Richartz (Phys. Rev. D 93, 064062
   (2016)), seeded by the extrapolation of the ladder roots to `κ = 0`. The
   root is accepted when it is stable under doubling of the recurrence depth,
   lies on the branch of the ladder and the family is not ambiguous. Otherwise
   it is returned with `status = :estimated` and
   `stop_reason = :extremal_root_off_ladder_branch`,
   `:extremal_family_ambiguous` or the reason the recurrence gives (for example
   `:extremal_depth_not_stable`), and `qnm` reports it as a `QNMEstimate`, or as
   a `QNMFailure` if its checks fail. When no angular eigenvalue can be
   computed at the recurrence root, the root has status `:failed` and `qnm`
   returns a `QNMFailure` with `stage = :root`.

If a ladder spin fails for an `|m| = l` co-rotating mode, the synchronous
endpoint is returned without the distance test. `qnm` still returns a
`QNMEndpointResult` with `status = :accepted`; such an unverified endpoint is
recognized by `stop_reason = :extremal_family_unverified` (and
`root.status = :estimated`).

The coefficient $P_3$ of the Richartz paper is used with the term
$-2im(1 + 6i\omega)$; the printed $-2im(1 - 6i\omega)$ does not satisfy the
Teukolsky equation for $m \neq 0$.

At `a = ±1` the direct route integrates along a complex path out of the
irregular singular horizon. At `a = 1` a Richardson route chosen from the start
(a mode outside the direct-route scope, or `excitation_options` that select it,
without a custom `radial_solver`) does the same; in the compact call without
`h` that integration also gives $\alpha$, and the step drift is the implied
root shift $|A^{\rm inc}/\alpha|/\max(1, |\omega|)$. On the Richardson route at
`a = -1`, and after a failed direct route, the amplitudes come from the
exact-extremal solver of `GSN_radial`. A failed incidence check never moves the
exact-extremal frequency.

## Reference

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
