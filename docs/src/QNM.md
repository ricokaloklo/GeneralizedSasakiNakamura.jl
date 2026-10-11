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
mode = QNMMode(-2, 2, 2, 0)
root = qnm_frequency(mode, 0.68)

sequence = qnm_sequence(mode, [0.0, 0.3, 0.68])
```

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
qnm_frequency(QNMMode(-2, 2, 2, 8), 0.68)                                    # 8₀
qnm_frequency(QNMMode(-2, 2, 2, 9), 0.68; convention=:complete_spectrum)     # 8₁
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

The continued-fraction inversion convention is available separately and
requires both an inversion index and a root guess:

```julia
leaver_root = qnm_frequency(
    mode, 0.68;
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
