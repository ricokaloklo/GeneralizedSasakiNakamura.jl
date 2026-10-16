# API Reference

This page documents the public interface intended for direct user calls. The package uses units with $G = c = M = 1$: frequencies are in units of $1/M$, radii in units of $M$, and the spin `a` is the dimensionless $a/M$ with $-1 \le a \le 1$.

## QNM Functions

The [Quasinormal modes](@ref) page documents these functions in full, with the overtone labels and the lower-level functions.

```julia
result = qnm(a, s, l, m, n)                # ordinary branch
result = qnm(a, s, l, m, n, mirror)        # mirror branch
pair = qnm_pair(a, s, l, m, n)             # both branches
root = qnm_frequency(a, s, l, m, n)        # root only: a LeaverResult, no amplitudes
roots = qnm_sequence(spins, s, l, m, n)    # Vector{LeaverResult} along a list of spins
# keywords: convention = :overtone | :complete_spectrum | :leaver, segment = 0 | 1 | :auto,
#   inversion_index; qnm and qnm_pair also take primary_guess, mirror_guess,
#   root_options, excitation_options and detailed; qnm_frequency and qnm_sequence take guess
```

`ordinary` (the default) is the mode with positive real frequency and `mirror`
the one with negative real frequency. `convention = :leaver` needs an explicit
frequency guess. `qnm` returns a `QNMResult` when every check passes, a
`QNMEstimate` when one accuracy check misses its target, a `QNMFailure` when a
step fails, and a `QNMEndpointResult`, with `status = :accepted` but no
amplitudes, for a zero-damping mode at $|a| = 1$; check `status` and the result
type before using the numbers. The result contains:

| field | meaning |
| :--- | :--- |
| `omega` or `frequency` | QNM frequency ($\operatorname{Im}\omega < 0$), or $m a/2$ for a zero-damping mode at $\lvert a\rvert = 1$ |
| `lambda` | separation constant $\lambda$ of the radial Teukolsky equation, i.e. the spin-weighted spheroidal eigenvalue (the same quantity as `Mode.lambda`) |
| `incidence_amplitude` | GSN incidence amplitude $A^{\rm inc}$ of the `IN` solution with unit transmission amplitude |
| `reflection_amplitude` | GSN reflection amplitude $A^{\rm ref}$ of the same solution |
| `incidence_derivative` | $\alpha = dA^{\rm inc}/d\omega$ at the root |
| `excitation_factor` | GSN excitation factor $A^{\rm ref}/(2\omega\alpha)$ |
| `status`, `stop_reason` | `status` is `:accepted`, `:estimated` or `:failed`; `stop_reason` names the outcome |

Set `detailed=true` to also get `result.X(rstar)` and `result.R(r)` (the `IN`
solutions), `result.Y(r)` (`IN` for `s = -2`, `UP` for `s = 2`, `missing`
otherwise) and `result.teukolsky.excitation_factor`, the excitation factor in
the Teukolsky normalization.

## Boundary and Normalization Constants

Perturbations are taken proportional to $e^{-i\omega t + im\phi}$ in
Boyer–Lindquist coordinates, so a damped mode has $\operatorname{Im}\omega < 0$.
Below, $r_*$ is the tortoise coordinate (next section),
$p = \omega - m\Omega_H$, $\Omega_H = a/(2r_+)$ and $\Delta = r^2 - 2r + a^2$.

### `BoundaryCondition`

| value | meaning |
| :--- | :--- |
| `IN` | purely ingoing at the horizon, $\propto e^{-ipr_*}$ |
| `UP` | purely outgoing at infinity, $\propto e^{+i\omega r_*}$ |
| `OUT` | purely outgoing at the horizon, $\propto e^{+ipr_*}$ |
| `DOWN` | purely ingoing at infinity, $\propto e^{-i\omega r_*}$ |

`OUT` and `DOWN` are built as linear combinations of `IN` and `UP`.

### `NormalizationConvention`

| value | meaning |
| :--- | :--- |
| `UNIT_GSN_TRANS` | the stored GSN solution is normalized to unit GSN transmission amplitude |
| `UNIT_TEUKOLSKY_TRANS` | the stored Teukolsky solution is normalized to unit Teukolsky transmission amplitude |

With unit transmission amplitude, the GSN functions behave as

```math
X^{\rm in} \to
\begin{cases}
e^{-ipr_*}, & r \to r_+,\\
A^{\rm inc} e^{-i\omega r_*} + A^{\rm ref} e^{+i\omega r_*}, & r \to \infty,
\end{cases}
\qquad
X^{\rm up} \to
\begin{cases}
C^{\rm ref} e^{-ipr_*} + C^{\rm inc} e^{+ipr_*}, & r \to r_+,\\
e^{+i\omega r_*}, & r \to \infty,
\end{cases}
```

and the Teukolsky functions as

```math
R^{\rm in} \to
\begin{cases}
\Delta^{-s} e^{-ipr_*}, & r \to r_+,\\
B^{\rm inc} r^{-1} e^{-i\omega r_*} + B^{\rm ref} r^{-2s-1} e^{+i\omega r_*}, & r \to \infty,
\end{cases}
\qquad
R^{\rm up} \to
\begin{cases}
C^{\rm inc} e^{+ipr_*} + C^{\rm ref} \Delta^{-s} e^{-ipr_*}, & r \to r_+,\\
r^{-2s-1} e^{+i\omega r_*}, & r \to \infty.
\end{cases}
```

`transmission_amplitude` holds the coefficient set to 1 (for the exceptions
see `GSNRadialFunction` and `Y_radial`), and `incidence_amplitude` and
`reflection_amplitude` hold $A$ (GSN) or $B$ (Teukolsky) for `IN` and $C$ for
`UP`, in the formalism of the object; the GSN and Teukolsky $C$ values differ.
For `OUT` and `DOWN`, the coefficient set to 1 is that of the purely outgoing
wave at the horizon (`OUT`) or of the purely ingoing wave at infinity (`DOWN`).
`incidence_amplitude` and `reflection_amplitude` then hold the coefficients of
the ingoing and outgoing waves at the horizon for `DOWN`, and of the outgoing
and ingoing waves at infinity for `OUT`, with the wave forms shown above
(arXiv:2101.04592, Eq. (93)).

## Tortoise Coordinate

`rstar_from_r` computes the tortoise coordinate (the formula in its docstring
below fixes the additive constant); `r_from_rstar` inverts it. A
`GSNRadialFunction` is evaluated at $r_*$, a `TeukolskyRadialFunction` and a
`YRadialFunction` at the Boyer–Lindquist $r$: `X(rstar_from_r(a, 10.0))` is the
GSN function at $r = 10$, `X(10.0)` its value at $r_* = 10$. For
solutions built by the ODE methods (`"linear"`, `"Riccati"`) at complex
`omega` and $|a| < 1$, `X` takes the real parameter $\rho$ of the rotated
integration path instead (see `GSN_radial`).

```@docs
rstar_from_r
r_from_rstar
```

## Homogeneous Radial Solutions

### `Teukolsky_radial`

Main call forms:

```julia
Rin, Rup = Teukolsky_radial(s, l, m, a, omega; method = "auto")
R = Teukolsky_radial(s, l, m, a, omega, IN; method = "auto")       # also UP, OUT, DOWN
R = Teukolsky_radial(s, l, m, a, omega, IN, rsin, rsout; method = "linear")
```

The tuple form returns the `IN` and `UP` solutions, the boundary-specific form
one solution. The forms without `rsin, rsout` evaluate negative `a` through the
symmetry $(m, a) \to (-m, -a)$.

Keywords accepted by each call form (any other keyword raises a `MethodError`):

| call form | keywords |
| :--- | :--- |
| `(s, l, m, a, omega)` | `method`, `tol`, `tolerance`, `xm`, `N`, `sfe`, `lfe`, `data_type`, `ODE_algorithm` |
| `(s, l, m, a, omega, bc)` | `method`, `tol`, `xm`, `N`, `sfe`, `lfe` |
| `(s, l, m, a, omega, bc, rsin, rsout)` | the keywords of the tuple form, plus `horizon_expansion_order`, `infinity_expansion_order` (series orders of the ODE methods) and `rsmp` (crossing point of the rotated path) |

Meaning of the keywords:

| keyword | default | meaning |
| :--- | :--- | :--- |
| `method` | `"auto"` | `"auto"` tries ISEM (iterative series expansion matching of the GSN equation) first; if ISEM throws and `xm` was not given, it logs a warning and falls back to the ODE methods. `"ISEM"` never falls back. `"linear"` and `"Riccati"` select the ODE methods, which integrate the GSN equation numerically (for `OUT` and `DOWN`, see below) |
| `tol` | `nothing` | tolerance of ISEM (chosen automatically when `nothing`; `1e-14` at real frequency), or `abstol = reltol` of the ODE methods (`1e-12` when `nothing`) |
| `tolerance` | `nothing` | the same setting as `tol`; `tol` takes precedence when both are given |
| `xm` | `nothing` | ISEM matching point in $x = (r - r_+)/(r - r_-)$, which is `0` at the horizon and `1` at infinity; `"auto"` and `"ISEM"` require `0 < xm < 1` |
| `N` | `nothing` | ISEM expansion order; chosen automatically when omitted |
| `sfe`, `lfe` | `nothing` | switches for the small-frequency expansion (MST, the Mano–Suzuki–Takasugi series) and the large-frequency (eikonal) expansion of ISEM; `nothing` selects them automatically |
| `data_type` | `ComplexF64` | number type of the ODE methods, e.g. `Complex{BigFloat}` |
| `ODE_algorithm` | `AutoVern9(Rosenbrock23(autodiff=false))` | solver of the ODE methods; to pass another, add and load DifferentialEquations.jl in your own environment |

At $|a| < 1$, a keyword described for ISEM only, or for the ODE methods only,
is ignored by the other methods. In the tuple form, `data_type` and
`ODE_algorithm` take effect only with `"linear"` or `"Riccati"`; in the form
with `rsin, rsout` also when `"auto"` falls back. Leave `xm` unset: an explicit
`xm` can introduce a cancellation error in the incidence amplitude that the
matching tolerance does not bound. The v0.9.0 keywords `rhom`, `TSinInf`, `TSoutInf`, `TSinHor` and
`TSoutHor` were removed.

`rsin` and `rsout` are the tortoise-coordinate ends of the numerical
integration done by the ODE methods. ISEM ignores them at $|a| < 1$. At
$|a| = 1$ exactly (not at spins such as `0.9999`) every method uses the
exact-extremal solver, which integrates between them as given. The forms
without them use `rsin = -50` and `rsout = 1000`, except at $|a| = 1$ and in
the tuple form with an ODE method at complex `omega`, where the domain is
chosen automatically. The ODE methods move `rsin` inward for $|a| > 0.9$ and
`rsout` outward for $|\omega| < 0.03$ (the GSN function records the adjusted
values). At complex `omega` they integrate along a rotated path in the complex
$r_*$ plane, and `rsin`, `rsout` bound the real parameter $\rho$ of that path.

For `OUT` and `DOWN` at $|a| < 1$, `"linear"` and `"Riccati"` still build the
`IN` and `UP` parts with `"auto"` (ISEM first, the ODE methods only as its
fallback).

Only `omega == 0` exactly returns the closed-form static solutions (incidence
and reflection amplitudes `missing`); small nonzero frequencies keep the
finite-frequency boundary conditions. For $|a| < 1$ at exactly
$\omega = m\Omega_H$, the `UP` solution is built with ISEM for every `method`;
it is logarithmic at the horizon, so its `incidence_amplitude` and
`reflection_amplitude` are coefficients of that logarithmic basis, not
comparable with plane-wave amplitudes.

Returned objects are `TeukolskyRadialFunction`s, callable at the
Boyer–Lindquist $r$:

```julia
Rin, Rup = Teukolsky_radial(-2, 2, 2, 0.9, 0.5)
Rin(10.0)                       # R at r = 10
Rin.Teukolsky_solution(10.0)    # R and dR/dr at r = 10
```

`Rin.Teukolsky_solution(r)` returns the vector `[R, dR/dr]` with the default
`method = "auto"` and with the ODE methods, and the tuple
`(R, dR/dr, d²R/dr², err)` with `method = "ISEM"` (for `err` see Accuracy).
At $|a| = 1$ and at `omega = 0` every method returns `[R, dR/dr]`; for
$|a| < 1$ the `UP` solution at $\omega = m\Omega_H$ returns the four-entry
tuple with every method.

```@docs
Teukolsky_radial
```

### `GSN_radial`

Main call forms:

```julia
Xin, Xup = GSN_radial(s, l, m, a, omega; method = "auto")
X = GSN_radial(s, l, m, a, omega, IN; method = "auto")       # also UP, OUT, DOWN
X = GSN_radial(s, l, m, a, omega, IN, rsin, rsout; method = "linear")
```

`GSN_radial` has the same call forms, methods and keywords as
`Teukolsky_radial`, plus `tolerance` in its boundary-specific form and `lambda`
(a precomputed spin-weighted spheroidal eigenvalue) in its form with
`rsin, rsout`. `lambda` is used only by `"linear"` and `"Riccati"` for `IN` and
`UP` at $|a| < 1$ and by the exact-extremal solver at $|a| = 1$; otherwise the
eigenvalue is computed from `s, l, m, a*omega`.

Returned objects are `GSNRadialFunction`s, callable at the tortoise coordinate
$r_*$:

```julia
Xin, Xup = GSN_radial(-2, 2, 2, 0.9, 0.5)
Xin(rstar_from_r(0.9, 10.0))                # X at r = 10
Xin.GSN_solution(rstar_from_r(0.9, 10.0))   # (X, dX/dr_*, err) at r = 10
```

`Xin.GSN_solution(rstar)` returns `(X, dX/dr_*, err)` when ISEM built the
solution or its `IN` and `UP` parts at $|a| < 1$ (it then takes the real $r_*$
also at complex `omega`; for `err` see Accuracy), and `(X, dX/dr_*)`
otherwise.

```@docs
GSN_radial
```

### `Y_radial`

```julia
Y = Y_radial(-2, l, m, a, omega, IN; method = "auto")
Y = Y_radial(2, l, m, a, omega, UP; method = "auto")
```

Constructs the radial function $Y$ that the point-particle convolutions
integrate against the source. Up to a constant factor,
$Y = e^{-i\phi}R^{\rm in}/r^2$ for `s = -2` and
$Y = e^{+i\phi}\Delta^2 R^{\rm up}/r^2$ for `s = +2`, where
$d\phi/dr = K/\Delta$ and $K = (r^2 + a^2)\omega - am$. Only these two branches
exist, `omega` must be nonzero, and there is no tuple form. The keywords are
`method`, `tol`, `xm`, `N`, `sfe` and `lfe`, with the meanings given for
`Teukolsky_radial`. With `method = "auto"` or `"ISEM"`, $Y$ is built with ISEM
(with the exact-extremal solver at $|a| = 1$), and a failure throws an error
without falling back; `"linear"` and `"Riccati"` build it from a numerically
integrated GSN function and accept complex `omega` only at $|a| = 1$.

The returned `YRadialFunction` is callable at the Boyer–Lindquist $r$:

| call | returns |
| :--- | :--- |
| `Y(r)` | $Y(r)$ |
| `Y.Y_solution(r)` | `(Y, dY/dr, X, err)`, where `X` is the GSN function at $r_*(r)$ and `err` is always `0` |
| `Y.X_solution(r)` | the GSN function at $r_*(r)$ |

Its amplitudes are the Teukolsky amplitudes of the solution obtained from the
GSN function: `transmission_amplitude` is the GSN-to-Teukolsky conversion
factor, not 1, although `normalization_convention` reads
`UNIT_TEUKOLSKY_TRANS`. For a real-frequency ISEM solution, the GSN function
used has the transmission amplitude `numerical_GSN_solution.solution_scale` of
the `GSN_radial` result for the same mode and keywords; this scale is 1 unless
the amplitudes approach the `Float64` range, and `transmission_amplitude` is
the conversion factor times it.

```@docs
Y_radial
YRadialFunction
```

### Accuracy

- `X.GSN_solution(rstar)` of an ISEM solution returns `err = 0.0` at real
  frequency, where no error estimate is computed (`0.0` does not mean exact),
  and an estimate of the local error, not a bound, at complex frequency.
- With `method = "ISEM"`, the `err` entry of `R.Teukolsky_solution(r)` is
  likewise an estimate at complex frequency; at real frequency it is `0.0`,
  except at large `r` on an `UP` solution and on the `OUT` and `DOWN`
  solutions built from it, where it is the truncation estimate of the
  asymptotic series used there.
- No error is returned for the amplitudes.

These estimates do not include the error of the angular eigenvalue or of the
equation coefficients near the horizon (up to about `1e-11`). Reflection
amplitudes that are small compared with the dominant amplitude, and solutions
at strongly damped frequencies, can lose relative accuracy. At a QNM
frequency, where $A^{\rm inc} = 0$, judge the computed $A^{\rm inc}$ by its
size relative to the largest of 1, $|A^{\rm ref}|$ and the transmission
amplitude, as `qnm` does.

Below `|omega| = 1e-8` the ODE methods lose the `UP` solution;
`method = "auto"` warns when it falls back to them there.

## Point-Particle Single-Mode Fluxes

### `Teukolsky_pointparticle_mode`

```julia
mode = Teukolsky_pointparticle_mode(s, l, m, n, k, a, p, e, x; method = "auto", N = -1, K = -1, Nmax = 2^12, Kmax = 2^9, tol = 1e-8)
```

Computes one inhomogeneous Teukolsky mode for a bound Kerr geodesic:

| parameter | meaning |
| :--- | :--- |
| `s` | spin weight; use `-2` for flux at infinity and `+2` for horizon flux |
| `l`, `m` | spheroidal harmonic indices |
| `n` | radial harmonic index |
| `k` | polar harmonic index |
| `a` | Kerr spin parameter |
| `p` | semi-latus rectum |
| `e` | eccentricity |
| `x` | inclination parameter, $x = \cos\theta_\mathrm{inc}$ |

With the default `method`, the radial and polar grids start from `N = 256` and `K = 64` and are doubled, up to `Nmax` and `Kmax`, until the relative change of the mode's energy flux between successive grids is at most `min(sample_tol, 3*sqrt(tol))` (`3e-4` by default) or the flux falls below `truncation_floor`. `N, K ≤ 0` select the default starting sizes; `N`, `K`, `Nmax` and `Kmax` must be powers of 2.

Method choices:

| method | meaning |
| :--- | :--- |
| `"auto"` | selects `"isem_trapezoidal"` |
| `"isem_trapezoidal"` | ISEM radial solution with adaptive trapezoidal convolution |
| `"isem_levin"` | ISEM radial solution with Levin convolution (circular equatorial orbits need no convolution) |
| `"trapezoidal"` | `Y` obtained by numerically integrating its radial equation, with trapezoidal convolution on one fixed `N × K` grid; `tol`, `Nmax`, `Kmax` and `sample_tol` do not apply |
| `"levin"` | as `"trapezoidal"`, with Levin convolution |

It returns a `TeukolskyPointParticleMode` with the fields:

| field | meaning |
| :--- | :--- |
| `mode` | `PointParticleMode` containing `s,l,m,n,k,a,omega,lambda` |
| `amplitude` | Teukolsky amplitude $Z$ at infinity for `s = -2`, or at the horizon for `s = +2` |
| `energy_flux` | energy flux of the mode |
| `angular_momentum_flux` | angular-momentum flux of the mode |
| `Carter_const_flux` | Carter-constant flux of the mode |
| `trajectory` | geodesic data used by the convolution |
| `Y_solution` | the radial function `Y` used by the convolution: a `YRadialFunction` for the ISEM methods, except that `"isem_trapezoidal"` (the default) stores a `Dict` of `Y` values sampled on its grid for eccentric and generic orbits; a numerically integrated solution for `"trapezoidal"` and `"levin"` |
| `SWSH` | spin-weighted spheroidal harmonic data |
| `method` | named tuple of the settings used: `method`, `radial_method` (always `"ISEM"`, also for `"trapezoidal"` and `"levin"`), `radial_sfe` (the small-frequency-expansion setting of the radial solve), the starting grid sizes `N` and `K`, and `truncation_floor` |

The particle mass $\mu$ is not an input: amplitudes are per unit $\mu$ and
fluxes per unit $\mu^2$. The mode frequency
$\omega = (m\Upsilon_\phi + n\Upsilon_r + k\Upsilon_\theta)/\Gamma$, built from
the Mino-time frequencies of the orbit, refers to the Boyer–Lindquist time $t$,
and the fluxes are $t$-rates of the energy, angular momentum and Carter
constant carried by the mode to infinity (`s = -2`) or into the horizon
(`s = +2`). At infinity the energy flux is
$|Z|^2/(4\pi\omega^2)$ and the angular-momentum flux $m|Z|^2/(4\pi\omega^3)$.
A horizon energy flux is negative for superradiant modes,
$\omega(\omega - m\Omega_H) < 0$.

Harmonics that the orbit does not have (`n ≠ 0` on circular orbits, `k ≠ 0` on
equatorial orbits) return zero amplitude and fluxes, with
`mode.mode.omega === nothing`, `mode.mode.lambda` `missing`,
`mode.method.radial_sfe` `nothing`, and `Y_solution` and `SWSH` set to
`nothing`. The mode `m = n = k = 0` also returns zero amplitude and fluxes and
`nothing` for `Y_solution` and `SWSH`, with `mode.mode.omega = 0`. The REPL
display prints the `mode` field as `Mode(...)` and adds `_inf` (`s = -2`) or
`_hor` (`s = +2`) to the amplitude and flux names (e.g. `amplitude_inf`); the
fields themselves have no suffix.

```@docs
Teukolsky_pointparticle_mode
```

### `GSN_pointparticle_mode`

```julia
mode = GSN_pointparticle_mode(s, l, m, n, k, a, p, e, x; method = "auto", N = -1, K = -1)
```

Computes the corresponding inhomogeneous GSN mode: the same fluxes, with the amplitude in the GSN normalization.

```@docs
GSN_pointparticle_mode
```

## Point-Particle Total Fluxes

### `Teukolsky_pointparticle_flux`

```julia
flux = Teukolsky_pointparticle_flux(a, p, e, x; tol = 1e-8, lmax = 30, nmax = 500, kmax = 20, truncation_floor = 1e-16)
```

Computes the total fluxes with a mode summation chosen by orbit type:

| orbit type | `orbit_type` | condition |
| :--- | :--- | :--- |
| circular equatorial | `:circular` | `e == 0` and `x == ±1` |
| eccentric equatorial | `:eccentric` | `e != 0` and `x == ±1` |
| circular inclined | `:inclined` | `e == 0` and `abs(x) != 1` |
| generic | `:generic` | otherwise |

Keywords:

| keyword | default | meaning |
| :--- | :--- | :--- |
| `tol` | `1e-8` | relative truncation tolerance of the mode sums, which add modes in shells of one ℓ (circular orbits), one k (inclined) or one n (eccentric, generic) |
| `lmax` | `30` | maximum $\ell$ index |
| `nmax` | `500` | maximum radial shell index |
| `kmax` | `20` | maximum polar shell index |
| `minimum_consecutive` | `nothing` (5 for circular equatorial orbits, 2 otherwise) | number of consecutive negligible terms that stops a sum over ℓ, m or (generic orbits) k; the outer `n` sum (eccentric, generic) and `k` sum (inclined) stop after `minimum_consecutive + 1` negligible shells |
| `N`, `N0` | `64` | initial radial grid interval count of the adaptive sampling (`N`, `K`, `Nmax` and `Kmax` must be powers of 2) |
| `K`, `K0` | `16` | initial polar grid interval count of the adaptive sampling |
| `Nmax` | `2^14` | maximum radial grid interval count |
| `Kmax` | `2^12` | maximum polar grid interval count |
| `sample_tol` | `1e-3` | adaptive single-mode sampling tolerance |
| `truncation_floor` | `1e-16` | absolute single-mode refinement floor; the total-flux summation tracks separate infinity and horizon floor values |
| `generic_beta` | `1e-4` | generic orbits, shells with `n ≥ 1`: raises the cutoff below which a mode counts as negligible to at least `generic_beta * tol` times the total flux summed so far (at infinity or at the horizon); `0` disables this |
| `generic_deep` | `1e-5` | generic orbits: an upward scan in ℓ stops as soon as one mode's energy flux is at most `generic_deep` times the cutoff; `0` disables this |
| `tail_levin` | `nothing` | eccentric and generic orbits: `true` evaluates the radial shells from `levin_nmin` on with adaptive Levin quadrature instead of a refined trapezoidal grid, `false` never does, `nothing` does so when it is much cheaper |
| `levin_nmin` | `50` | first radial shell that `tail_levin` applies to |
| `levin_local_n` | `16` | Chebyshev nodes per adaptive Levin segment (for generic orbits, intervals per direction) |
| `levin_max_depth` | `8` | maximum adaptive Levin subdivision depth |
| `info` | `false` | no effect; kept for compatibility |
| `record` | `false` | `true` writes every evaluated mode to an HDF5 file (see [HDF5 Mode Records](@ref)) |
| `record_path` | `nothing` | optional HDF5 output path |
| `fast` | `true` | `false` turns `tail_levin` off and computes generic-orbit modes without reusing cached orbit samples |

It accepts stable bound orbits and marginally stable (critical) bound circular orbits, as classified by `KerrGeodesics`; for other parameters it warns and returns `nothing`.

It returns a `TeukolskyPointParticleFlux`:

| field | meaning |
| :--- | :--- |
| `a`, `p`, `e`, `x` | input orbital parameters |
| `orbit_type` | `:circular`, `:eccentric`, `:inclined` or `:generic` |
| `infinity_energy_flux` | total energy flux at infinity |
| `infinity_angular_momentum_flux` | total angular-momentum flux at infinity |
| `infinity_carter_constant_flux` | total Carter-constant flux at infinity |
| `horizon_energy_flux` | total energy flux at the horizon |
| `horizon_angular_momentum_flux` | total angular-momentum flux at the horizon |
| `horizon_carter_constant_flux` | total Carter-constant flux at the horizon |
| `total_modes` | number of single modes evaluated |
| `tolerance` | the `tol` used; when a sum stopped at `nmax`, the display adds the last shell's relative energy and a larger `nmax` to try |
| `truncation_floor` | final absolute truncation floors, as `(infinity = …, horizon = …)` |
| `reached` | highest shell indices reached: `l_reached_inf`, `l_reached_hor` and their maximum `l_reached` (circular), `k_reached_inf`, `k_reached_hor` (inclined), or `n_reached_inf`, `n_reached_hor` (eccentric, generic); `nothing` in a `k_reached_*` or `n_reached_*` entry means the sum reached `kmax` or `nmax` before converging. The display groups the `l` and `n` entries as `l_reached = (infinity = …, horizon = …)` and `n_reached = (…)` |
| `reached.convolution_integral` | named tuple describing whether convolution integrals were needed and which tail strategy was used |
| `cost` | wall-clock runtime in seconds |
| `result` | raw mode-summation result |

Example:

```julia
flux = Teukolsky_pointparticle_flux(0.9, 6.0, 0.0, 1.0)   # circular equatorial, about 0.3 s after warm-up
flux.infinity_energy_flux
flux.horizon_energy_flux
flux.reached.l_reached_inf
```

A generic orbit such as `(0.9, 6.0, 0.7, 0.5)` sums about 1.3 million modes;
the run shown on the Home page took about 22 minutes at `tol = 1e-8`.

```@docs
Teukolsky_pointparticle_flux
```

## HDF5 Mode Records

With `record = true`, `Teukolsky_pointparticle_flux` writes every evaluated mode
of an eccentric, inclined or generic orbit to an HDF5 file; circular equatorial
orbits are not recorded. Without `record_path` the file is created in the
current directory with a name built from the orbit, e.g.
`eccentric_mode_data_a_..._p_..._e_....h5`. An existing file at the output path
is overwritten.

Each mode is stored in its own group:

| orbit | group of one mode, e.g. |
| :--- | :--- |
| eccentric | `n_000003/m_+00002/l_0002/I` |
| inclined | `k_000001/m_+00002/l_0002/I` |
| generic | `n_000003/k_+00001/m_+00002/l_0002/I` |

The last level is `I` for the solution at infinity (`s = -2`) and `H` for the
horizon (`s = +2`). Each such group holds the mode indices, the amplitude split
into `Amplitude_re` and `Amplitude_im`, the fluxes `EnergyFlux`,
`AngularMomentumFlux` and `CarterConstantFlux`, and the grid size used
(`N_sample` for eccentric orbits, `K_sample` for inclined orbits, both for
generic orbits). The group `meta` holds the inputs and the totals.

To read the file, add `HDF5.jl` to your own environment first
(`using Pkg; Pkg.add("HDF5")`):

```julia
using HDF5

flux = Teukolsky_pointparticle_flux(0.9, 6.0, 0.3, 1.0; record = true,
    record_path = "eccentric_mode_data.h5")
h5open("eccentric_mode_data.h5", "r") do f
    keys(f)
end
```

## Returned Types

`TeukolskyPointParticleMode`, `GSNPointParticleMode`, `PointParticleMode` and
`TeukolskyPointParticleFlux` are not exported; to name them, qualify them, e.g.
`GeneralizedSasakiNakamura.TeukolskyPointParticleFlux`.

### `Mode`

Stores homogeneous mode metadata:

| field | meaning |
| :--- | :--- |
| `s` | spin weight |
| `l` | harmonic index $\ell$ |
| `m` | azimuthal index $m$ |
| `a` | Kerr spin parameter |
| `omega` | frequency |
| `lambda` | spin-weighted spheroidal eigenvalue $\lambda$, the separation constant of the radial Teukolsky equation |

### `PointParticleMode`

Stores point-particle mode metadata:

| field | meaning |
| :--- | :--- |
| `s`, `l`, `m` | spin and spheroidal harmonic indices |
| `n`, `k` | radial and polar harmonic indices |
| `a` | Kerr spin parameter |
| `omega` | mode frequency; `nothing` for a harmonic that the orbit does not have |
| `lambda` | spin-weighted spheroidal eigenvalue; `missing` for such a harmonic |

### `TeukolskyRadialFunction`

Stores a homogeneous Teukolsky radial solution.

| field | meaning |
| :--- | :--- |
| `mode` | `Mode` metadata |
| `boundary_condition` | one of `IN`, `UP`, `OUT`, `DOWN` |
| `transmission_amplitude` | Teukolsky transmission amplitude |
| `incidence_amplitude` | Teukolsky incidence amplitude |
| `reflection_amplitude` | Teukolsky reflection amplitude |
| `P_solution` | not for direct use; `missing` except for ISEM solutions at complex frequency |
| `GSN_solution` | the `GSNRadialFunction` the solution was built from; `missing` for the static solution and, at $\lvert a\rvert < 1$, for `OUT` and `DOWN` built with `method = "ISEM"` |
| `Teukolsky_solution` | callable radial solution |
| `normalization_convention` | `UNIT_TEUKOLSKY_TRANS` |

### `GSNRadialFunction`

Stores a homogeneous GSN radial solution.

| field | meaning |
| :--- | :--- |
| `mode` | `Mode` metadata |
| `boundary_condition` | one of `IN`, `UP`, `OUT`, `DOWN` |
| `rsin`, `rsout` | integration ends of the ODE methods after their widening (bounds of $\rho$ at complex frequency) or of the exact-extremal solver; `missing` for ISEM solutions at $\lvert a\rvert < 1$ and for the static solution |
| `rsmp` | for ISEM `IN` and `UP` solutions at $\lvert a\rvert < 1$, $r_*$ of the matching point (it can be `missing` at complex frequency); for the ODE methods at complex frequency, the $r_*$ where the rotated path crosses the real axis; otherwise `nothing` or `missing` |
| `horizon_expansion_order`, `infinity_expansion_order` | orders of the series at the horizon and at infinity used by the method; for ISEM at complex frequency the order `N` when given, otherwise `missing`; `missing` for the static solution |
| `transmission_amplitude` | GSN transmission amplitude: 1, except for a real-frequency ISEM solution whose rescaling to unit transmission would overflow, which keeps its scale here (then `numerical_GSN_solution.normalization == :scaled_gsn_transmission`); `missing` for the static solution |
| `incidence_amplitude` | GSN incidence amplitude |
| `reflection_amplitude` | GSN reflection amplitude |
| `numerical_GSN_solution` | ODE solution for the ODE methods; matching metadata for ISEM |
| `numerical_Riccati_solution` | Riccati ODE solution when applicable |
| `GSN_solution` | callable GSN radial solution |
| `normalization_convention` | `UNIT_GSN_TRANS`; `UNIT_TEUKOLSKY_TRANS` for the static solution at `omega = 0` |
| `method` | the method that built the solution: `"ISEM"` (also for the exact-extremal solver under `"auto"` and `"ISEM"`), `"linear"`, `"Riccati"` or `"static"`; at $\lvert a\rvert < 1$, `OUT` and `DOWN` record `"linear"` or `"Riccati"` when requested and `"auto"` after a fallback |

### `TeukolskyPointParticleMode` and `GSNPointParticleMode`

Both have the fields listed under `Teukolsky_pointparticle_mode`.

### `TeukolskyPointParticleFlux`

Stores total flux results from `Teukolsky_pointparticle_flux`; see [Point-Particle Total Fluxes](@ref).

## Parallel Calls

Start Julia with `julia --threads=auto`, and avoid nested BLAS parallelism:

```julia
using GeneralizedSasakiNakamura, LinearAlgebra
BLAS.set_num_threads(1)

parameters = [(-2, 2, 2, 0.7, 0.3), (-2, 2, 2, 0.9, 0.5)]
solutions = Vector{Any}(undef, length(parameters))
Threads.@threads for i in eachindex(parameters)
    solutions[i] = GSN_radial(parameters[i]..., IN)
end
```

Radial solutions (built, or one shared solution evaluated, from several
threads), QNM frequencies and excitation factors are thread-safe and give
identical results in serial and threaded runs; point-particle modes and fluxes
are thread-safe as well. If you construct a Leaver-method harmonic yourself with
SpinWeightedSpheroidalHarmonics 1.4 and share it between threads, do not
evaluate its first polar derivative from several threads at once.

### Point-particle submissions

Several single-mode or total-flux calls on the same orbit can share one
read-only trajectory and sampling geometry:

```julia
modes = with_pointparticle_submission(0.9, 10.0, 0.2, 0.3) do submission
    tasks = [spawn_pointparticle_submission!(submission) do
                 Teukolsky_pointparticle_mode(-2, 2, 2, n, 0, 0.9, 10.0, 0.2, 0.3; submission)
             end for n in 0:2]
    fetch.(tasks)
end
```

| function | meaning |
| :--- | :--- |
| `with_pointparticle_submission(f, a, p, e, x; Nmax = 2^14, Kmax = 2^12, stop_requested = nothing)` | create a submission for one orbit and sampling grid, call `f(submission)`, then wait for every admitted mode and task before releasing it |
| `spawn_pointparticle_submission!(f, submission)` | start a task registered with the submission; it is joined before the submission is released |
| `cancel_pointparticle_submission!(submission; reason)` | stop admitting new modes; modes already running finish normally |

```@docs
with_pointparticle_submission
spawn_pointparticle_submission!
cancel_pointparticle_submission!
```

`Teukolsky_pointparticle_mode`, `GSN_pointparticle_mode` and
`Teukolsky_pointparticle_flux` accept `submission` and `stop_requested`
(a zero-argument function returning `true` to request cancellation). A
submission belongs to one orbit `(a, p, e, x)` and one `(Nmax, Kmax)`; passing it
to a call with other values is an error. Single-mode calls that pass a
submission therefore use its `Nmax` and `Kmax` (by default `2^14` and `2^12`)
instead of their own defaults (`2^12` and `2^9`); create the submission with
`Nmax = 2^12, Kmax = 2^9` to reproduce standalone results. Returned mode and
flux objects remain valid after the submission is released. After
cancellation, a call that asks the submission to admit a new mode or task
throws `GeneralizedSasakiNakamura.ModeSummation.PublicSubmissionCancelled`.

For process-based batches, start the workers and load the package on each of
them:

```julia
using Distributed
addprocs(4)                                   # or start workers with the cluster launcher
@everywhere using GeneralizedSasakiNakamura
parameters = [(-2, 2, 2, 0.7, 0.3), (-2, 2, 2, 0.9, 0.5)]
solutions = pmap(p -> GSN_radial(p..., IN), parameters)
```

Each worker has its own package state and compiles the package on its first
call, so time warm calls separately from package loading.
