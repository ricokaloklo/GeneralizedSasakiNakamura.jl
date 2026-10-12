# API Reference

This page documents the public interface intended for direct user calls. The package uses units with $G = c = M = 1$.

## Exported Functions

The native QNM frequency, independent ISEM validation, and excitation-factor
interfaces are documented on the [Quasinormal modes](@ref) page.

```@docs
Teukolsky_radial
```

```@docs
GSN_radial
```

```@docs
Teukolsky_pointparticle_mode
```

```@docs
GSN_pointparticle_mode
```

```@docs
Teukolsky_pointparticle_flux
```

```@docs
rstar_from_r
```

```@docs
r_from_rstar
```

## Quasinormal Modes

```julia
result = qnm(a, s, l, m, n)
result = qnm(a, s, l, m, n, mirror)
pair = qnm_pair(a, s, l, m, n)
root = qnm_frequency(a, s, l, m, n)                 # frequency only
roots = qnm_sequence(spins, s, l, m, n)
# keywords: convention=:overtone | :complete_spectrum, segment=0 | 1 | :auto
```

`ordinary` is the default branch. `mirror` requests the negative-real root
and its amplitudes. The compact result contains:

| field | meaning |
| :--- | :--- |
| `omega` or `frequency` | polished Leaver QNM frequency |
| `lambda` | radial Teukolsky separation constant |
| `incidence_amplitude` | GSN incoming amplitude |
| `reflection_amplitude` | GSN reflected amplitude |
| `incidence_derivative` | Richardson-extrapolated $dA_{\rm in}/d\omega$ |
| `excitation_factor` | GSN excitation factor |
| `status`, `stop_reason` | combined root, radial, and derivative status |

Set `detailed=true` to construct callable `X`, `Y`, and `R` solutions. Root
label conventions, result types, and lower-level interfaces are documented on
the [Quasinormal modes](@ref) page.

## Boundary and Normalization Constants

### `BoundaryCondition`

Boundary conditions are represented by the constants below.

| value | meaning |
| :--- | :--- |
| `IN` | purely ingoing at the horizon |
| `UP` | purely outgoing at infinity |
| `OUT` | purely outgoing at the horizon |
| `DOWN` | purely ingoing at infinity |

### `NormalizationConvention`

| value | meaning |
| :--- | :--- |
| `UNIT_GSN_TRANS` | the stored GSN solution is normalized to unit GSN transmission amplitude |
| `UNIT_TEUKOLSKY_TRANS` | the stored Teukolsky solution is normalized to unit Teukolsky transmission amplitude |

## Homogeneous Radial Solutions

### `Teukolsky_radial`

Main call forms:

```julia
Rin, Rup = Teukolsky_radial(s, l, m, a, omega; method = "auto")
R = Teukolsky_radial(s, l, m, a, omega, IN; method = "auto")
R = Teukolsky_radial(s, l, m, a, omega, UP; method = "auto")
```

The tuple form returns the `IN` and `UP` solutions. The boundary-specific form also supports `OUT` and `DOWN` where available.

Keyword summary:

| keyword | default | meaning |
| :--- | :--- | :--- |
| `method` | `"auto"` | `"auto"` tries GSN-ISEM first and falls back to the legacy automatic route only on direct failure with no explicit `xm`; `"GSN-ISEM"` is strict; `"ISEM"` selects the previous implementation; `"linear"` and `"Riccati"` select the legacy GSN ODE solvers |
| `tol` | internal default | ISEM matching tolerance or alias for the ODE tolerance in high-level calls |
| `tolerance` | internal default | legacy ODE tolerance; `tol` takes precedence when both are supplied |
| `xm` | `nothing` | matching point: `0 < xm < 1` for GSN-ISEM; a negative matching coordinate for the original ISEM |
| `rhom` | `nothing` | original ISEM matching control; use with `method = "ISEM"` |
| `N` | `nothing` | optional ISEM expansion order override; when omitted, ISEM uses selector/default controls plus local-`N` rescue |
| `sfe`, `lfe` | `nothing` | optional small/large-frequency expansion switches |
| `TSinInf`, `TSoutInf`, `TSinHor`, `TSoutHor` | `nothing` | original ISEM identity switches; use with `method = "ISEM"` |

Omit optional matching controls to use automatic selection. Do not transfer an
`xm` value between `"GSN-ISEM"` and `"ISEM"`; their coordinates differ.
If construction at an explicit `xm` fails, `"auto"` reports the failure rather
than falling back to a legacy solver that ignores that junction.
Automatic matching is recommended. For the tested UP mode
`(-2, 2, 2, 0.7, 0.4198385710926635 - 2.09355961067422im)`,
explicit `xm = 0.4` gives an incidence error of `1.8151e-5` in reflection-amplitude
units (`5.048e-10` in transmission-amplitude units); matching tolerance alone
does not bound this cancellation error.
Original-ISEM-only controls require `method = "ISEM"`; supplying them to
`"auto"` or `"GSN-ISEM"` raises an `ArgumentError`.

Only `omega == 0` uses the analytic static branch. Small nonzero frequencies
retain their finite-frequency boundary conditions. The dedicated UP
horizon-threshold branch requires exact equality with `omega = m a / (2 r_+)`;
nearby frequencies are not substituted by that value. Likewise, exact-extremal
dispatch requires `a == 1` or `a == -1`.

Returned objects are `TeukolskyRadialFunction`. They are callable:

```julia
Rin, Rup = Teukolsky_radial(-2, 2, 2, 0.9, 0.5)
Rin(10.0)
Rin.Teukolsky_solution(10.0)
```

`TeukolskyRadialFunction(r)` returns only the radial function value. `Teukolsky_solution(r)` returns the stored radial data from the selected backend.

### `GSN_radial`

Main call forms:

```julia
Xin, Xup = GSN_radial(s, l, m, a, omega; method = "auto")
X = GSN_radial(s, l, m, a, omega, IN; method = "auto")
X = GSN_radial(s, l, m, a, omega, UP; method = "auto")
```

`method = "auto"` tries the direct GSN-ISEM route first and falls back to the legacy automatic radial route only after a direct-construction failure with no explicit `xm`. `method = "GSN-ISEM"` is strict and does not fall back. The former `"direct_ISEM"` spelling remains a strict compatibility alias. The previous implementation remains available as `"ISEM"`, and `"linear"` and `"Riccati"` select the legacy GSN ODE evolution. Matching controls, the static branch, and the superradiance-threshold branch follow the same rules as `Teukolsky_radial`.

Returned objects are `GSNRadialFunction`. They are callable:

```julia
Xin, Xup = GSN_radial(-2, 2, 2, 0.9, 0.5)
Xin(rstar_from_r(0.9, 10.0))
Xin.GSN_solution(rstar_from_r(0.9, 10.0))
```

`GSNRadialFunction(rstar)` returns only the GSN function value. `GSN_solution(rstar)` returns the stored GSN data from the selected backend.

Keyword summary (in addition to the matching controls shared with `Teukolsky_radial`):

| keyword | default | meaning |
| :--- | :--- | :--- |
| `tol` | `nothing` | tolerance; alias for `tolerance` that takes precedence when both are supplied |
| `tolerance` | `nothing` (`1e-12` for the legacy ODE paths) | legacy ODE tolerance |
| `lambda` | `nothing` | precomputed spin-weighted spheroidal eigenvalue (method with `rsin`, `rsout`); computed from `s, l, m, a*omega` when omitted |
| `ODE_algorithm` | `AutoVern9(Rosenbrock23(autodiff=false))` | ODE solver of the legacy `"linear"` and `"Riccati"` paths |
| `data_type` | `ComplexF64` | number type of the legacy ODE paths |

### `Y_radial`

```julia
Y = Y_radial(s, l, m, a, omega, IN; method = "auto")
```

Constructs the auxiliary ISEM radial function used by the point-particle convolutions, with the same matching controls as `Teukolsky_radial`. The returned `YRadialFunction` is callable.

```@docs
Y_radial
YRadialFunction
```

## Point-Particle Single-Mode Fluxes

### `Teukolsky_pointparticle_mode`

```julia
mode = Teukolsky_pointparticle_mode(s, l, m, n, k, a, p, e, x; method = "auto", N = -1, K = -1, Nmax = 2^12, Kmax = 2^9, tol = 1e-8)
```

With the default `method`, the radial and polar grids start from `N = 256` and `K = 64` and are doubled until the mode converges to `tol`, up to `Nmax` and `Kmax`.

Computes one inhomogeneous Teukolsky mode for a bound Kerr geodesic specified by:

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

Method choices:

| method | meaning |
| :--- | :--- |
| `"auto"` | selects `"isem_trapezoidal"` |
| `"isem_trapezoidal"` | ISEM radial solution with adaptive trapezoidal convolution |
| `"isem_levin"` | ISEM radial solution with Levin convolution where supported |
| `"trapezoidal"` | legacy radial solution with trapezoidal convolution |
| `"levin"` | legacy radial solution with Levin convolution |

The return type is `TeukolskyPointParticleMode`. Its main fields are:

| field | meaning |
| :--- | :--- |
| `mode` | `PointParticleMode` containing `s,l,m,n,k,a,omega,lambda` |
| `amplitude` | Teukolsky amplitude at infinity for `s = -2`, or at the horizon for `s = +2` |
| `energy_flux` | single-mode energy flux |
| `angular_momentum_flux` | single-mode angular-momentum flux |
| `Carter_const_flux` | single-mode Carter-constant flux |
| `trajectory` | geodesic data used by the convolution |
| `Y_solution` | auxiliary ISEM radial solution used by the convolution |
| `SWSH` | spin-weighted spheroidal harmonic data |
| `method` | named tuple recording the method and grid sizes |

### `GSN_pointparticle_mode`

```julia
mode = GSN_pointparticle_mode(s, l, m, n, k, a, p, e, x; method = "auto", N = -1, K = -1)
```

Computes the corresponding inhomogeneous GSN mode. Fluxes are formalism-independent; the amplitude is expressed in the GSN normalization.

## Point-Particle Total Fluxes

### `Teukolsky_pointparticle_flux`

```julia
flux = Teukolsky_pointparticle_flux(a, p, e, x; tol = 1e-8, lmax = 30, nmax = 500, kmax = 20, truncation_floor = 1e-16)
```

Computes total point-particle fluxes by selecting the appropriate mode-summation strategy from the orbit type:

| orbit type | condition | summation |
| :--- | :--- | :--- |
| circular equatorial | `e == 0` and `x == ±1` | circular mode summation |
| eccentric equatorial | `e != 0` and `x == ±1` | eccentric mode summation |
| circular inclined | `e == 0` and `abs(x) != 1` | inclined mode summation |
| generic | otherwise | generic mode summation |

For `x == -1`, the equatorial retrograde case is internally mapped to the corresponding positive-inclination convention by flipping the sign of `a`.

Important keywords:

| keyword | default | meaning |
| :--- | :--- | :--- |
| `tol` | `1e-8` | global shell truncation tolerance |
| `lmax` | `30` | maximum $\ell$ index |
| `nmax` | `500` | maximum radial shell index |
| `kmax` | `20` | maximum polar shell index |
| `minimum_consecutive` | `nothing` (5 for circular equatorial orbits, 2 otherwise) | number of consecutive small shells required for truncation |
| `N`, `N0` | `64` | initial radial grid interval count for adaptive sampling |
| `K`, `K0` | `16` | initial polar grid interval count for adaptive sampling |
| `Nmax` | `2^14` | maximum radial grid interval count |
| `Kmax` | `2^12` | maximum polar grid interval count |
| `sample_tol` | `1e-3` | adaptive single-mode sampling tolerance |
| `truncation_floor` | `1e-16` | absolute single-mode refinement floor; the total-flux summation tracks separate infinity and horizon floor values |
| `generic_beta` | `1e-4` | generic orbits: per-mode budget floor `generic_beta * tol * abs(branch_total)` for outward ℓ tails; `0` disables it |
| `generic_deep` | `1e-5` | generic orbits: an outward ℓ tail stops once a mode is at most this fraction of its cutoff; `0` disables it. Peak-search walks keep the original cutoff and confirmation count |
| `tail_levin` | `nothing` (automatic) | adaptive Levin quadrature for high-index eccentric and generic radial tails: `true`, `false`, or automatic selection |
| `levin_nmin` | `50` | radial shell at which tail Levin quadrature becomes eligible |
| `levin_local_n` | `ConvolutionIntegrals.DEFAULT_ADAPTIVE_LEVIN_LOCAL_N` | number of local nodes per adaptive Levin subinterval |
| `levin_max_depth` | `8` | maximum adaptive Levin subdivision depth |
| `info` | `false` | print solver fallback diagnostics from `Y_radial` when `method = "auto"` changes radial backend |
| `record` | `false` | when set to `true`, records mode details to HDF5 where supported |
| `record_path` | `nothing` | optional HDF5 output path |
| `fast` | `true` | use cached/presampled fast summation path where available |

If the supplied orbital parameters do not define a bound orbit according to `KerrGeodesics`, the function emits a warning and returns `nothing`.

The public tolerance keyword is `tol`. The returned object displays this value as `tolerance`.

For high-index tail modes, eccentric and generic summations can use adaptive Levin quadrature instead of uniformly increasing a trapezoidal grid. The adaptive rule refines radial phase intervals only where the oscillatory integral has not stabilized. In generic two-dimensional convolutions, this radial adaptive Levin rule is paired with fixed Clenshaw-Curtis sampling in the polar direction, so the smooth polar dependence is resolved with a compact cosine-spaced grid while radial oscillations receive targeted refinement.

The return type is `TeukolskyPointParticleFlux`:

| field | meaning |
| :--- | :--- |
| `a`, `p`, `e`, `x` | input orbital parameters |
| `orbit_type` | detected orbit type as a `Symbol` |
| `infinity_energy_flux` | total energy flux at infinity |
| `infinity_angular_momentum_flux` | total angular-momentum flux at infinity |
| `infinity_carter_constant_flux` | total Carter-constant flux at infinity |
| `horizon_energy_flux` | total energy flux at the horizon |
| `horizon_angular_momentum_flux` | total angular-momentum flux at the horizon |
| `horizon_carter_constant_flux` | total Carter-constant flux at the horizon |
| `total_modes` | number of single modes evaluated |
| `tolerance` | effective summation tolerance |
| `truncation_floor` | named tuple of final infinity and horizon absolute truncation floors |
| `reached` | named tuple of reached shell indices, e.g. `l_reached`, `n_reached`, or `k_reached` with separate infinity and horizon entries where applicable |
| `reached.convolution_integral` | named tuple describing whether convolution integrals were needed and which tail strategy was used |
| `cost` | wall-clock runtime in seconds |
| `result` | raw mode-summation result |

Example:

```julia
flux = Teukolsky_pointparticle_flux(0.9, 6.0, 0.7, cos(pi / 4))
flux.infinity_energy_flux
flux.horizon_energy_flux
```

## HDF5 Mode Records

When `record = true` is used in the mode-summation interface, mode details are written to an HDF5 file. Records are stored hierarchically by shell indices and mode indices, with datasets for amplitudes, fluxes, and grid sizes. Complex quantities are stored by real and imaginary parts for portable HDF5 access.

Use `HDF5.jl` to inspect the file:

```julia
using HDF5

h5open("eccentric_mode_data.h5", "r") do f
    keys(f)
end
```

## Returned Types

### `Mode`

Stores homogeneous mode metadata:

| field | meaning |
| :--- | :--- |
| `s` | spin weight |
| `l` | harmonic index $\ell$ |
| `m` | azimuthal index $m$ |
| `a` | Kerr spin parameter |
| `omega` | frequency |
| `lambda` | spin-weighted spheroidal eigenvalue |

### `PointParticleMode`

Stores point-particle mode metadata:

| field | meaning |
| :--- | :--- |
| `s`, `l`, `m` | spin and spheroidal harmonic indices |
| `n`, `k` | radial and polar harmonic indices |
| `a` | Kerr spin parameter |
| `omega` | mode frequency |
| `lambda` | spin-weighted spheroidal eigenvalue |

### `TeukolskyRadialFunction`

Stores a homogeneous Teukolsky radial solution.

| field | meaning |
| :--- | :--- |
| `mode` | `Mode` metadata |
| `boundary_condition` | one of `IN`, `UP`, `OUT`, `DOWN` |
| `transmission_amplitude` | Teukolsky transmission amplitude |
| `incidence_amplitude` | Teukolsky incidence amplitude |
| `reflection_amplitude` | Teukolsky reflection amplitude |
| `P_solution` | internal ISEM `P` solution when available; this is an internal callable, not the public rescue interface |
| `GSN_solution` | associated `GSNRadialFunction` when available |
| `Teukolsky_solution` | callable radial solution |
| `normalization_convention` | normally `UNIT_TEUKOLSKY_TRANS` |

### `GSNRadialFunction`

Stores a homogeneous GSN radial solution.

| field | meaning |
| :--- | :--- |
| `mode` | `Mode` metadata |
| `boundary_condition` | one of `IN`, `UP`, `OUT`, `DOWN` |
| `rsin`, `rsout`, `rsmp` | legacy ODE integration and matching coordinates, or `missing` for ISEM |
| `horizon_expansion_order` | legacy ODE horizon expansion order, or `missing` for ISEM |
| `infinity_expansion_order` | legacy ODE infinity expansion order, or `missing` for ISEM |
| `transmission_amplitude` | GSN transmission amplitude |
| `incidence_amplitude` | GSN incidence amplitude |
| `reflection_amplitude` | GSN reflection amplitude |
| `numerical_GSN_solution` | ODE solution for legacy methods; matching metadata for `method == "ISEM"` or `method == "GSN-ISEM"` |
| `numerical_Riccati_solution` | Riccati ODE solution when applicable |
| `GSN_solution` | callable GSN radial solution |
| `normalization_convention` | normally `UNIT_GSN_TRANS` |
| `method` | solver method string |

Reported errors are estimates, not complete error bounds; they exclude equation-coefficient evaluation error, which can be of order `1e-11` near the horizon.

They also exclude uncertainty in the input angular eigenvalue. In tested real-frequency MST amplitudes, relative errors reach about `3e-10`, with estimates low by up to a factor of `2.6`.

Complex-frequency amplitude checks show residual underestimation by up to about `1.33` in the tested cases. At an incidence zero, use absolute or transmission-normalized error rather than relative error with respect to the zero. Real-frequency estimates can miss larger relative errors in small reflection amplitudes.

Amplitude estimates do not bound errors in `X`, `R`, or `Y` evaluations. Even at real or weakly damped frequencies, transmission-normalized `R(10)` can have relative errors above `1e-8`; strong damping can cause order-unity errors on the real radius axis despite successful construction. Small reflection amplitudes may also lose relative precision. No total radial-value error estimate is provided.

For example, tested real-frequency cases with `s = -2`, `l = m = 2:4`, and `abs(omega) = 3` have reflection amplitudes with no reliable relative digits. These are not high-accuracy reflection results even when the dominant amplitude is accurate.

Explicit legacy methods have separate low-frequency limits. At
`(s,l,m,a,omega) = (2,2,-1,0.47302630436742865,-1e-4)`, `method="ISEM"`
returns an unreliable `R(6)` (reference-relative error about `674`; the
error slot is about `1.08` in returned-value units). `method="Riccati"`
agrees with the default near-zone value to about `5e-14`, but its incidence
amplitude is wrong with the default outer boundary. An accurate local
radial value does not certify the legacy far-field amplitudes.

### `TeukolskyPointParticleMode` and `GSNPointParticleMode`

These store single-mode inhomogeneous results. `TeukolskyPointParticleMode` stores the amplitude in the Teukolsky formalism; `GSNPointParticleMode` stores the amplitude in the GSN formalism. Flux fields are shared in meaning.

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

Serial/thread bitwise checks cover radial construction, complex-frequency
evaluation, QNM frequencies and excitation factors, and concurrent evaluation
of a shared radial solution. Point-particle modes and fluxes also depend on the
angular and geodesic packages. KerrGeodesics 0.5 evaluates its Jacobi elliptic
functions without a shared workspace. With SpinWeightedSpheroidalHarmonics 1.4.0,
do not take the first polar derivative of one shared Leaver-method harmonic
object from several threads at once; harmonics constructed inside each call are
not shared.

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
| `with_pointparticle_submission(f, a, p, e, x; Nmax = 2^14, Kmax = 2^12, stop_requested = nothing)` | create a submission for one orbit and sampling grid, call `f(submission)`, then wait for every admitted mode and task before releasing internal caches |
| `spawn_pointparticle_submission!(f, submission)` | start a task registered with the submission; it is joined before the submission is released |
| `cancel_pointparticle_submission!(submission; reason)` | stop admitting new modes; modes already running finish normally |

`Teukolsky_pointparticle_mode`, `GSN_pointparticle_mode` and
`Teukolsky_pointparticle_flux` accept `submission` and `stop_requested`
(a zero-argument function returning `true` to request cancellation). A
submission belongs to one orbit `(a, p, e, x)` and one `(Nmax, Kmax)`; passing it
to a call with other values is an error. Each mode keeps its own mutable
workspace; returned mode and flux objects remain valid after the submission
is released.

For process-based batches, load the package once on each worker:

```julia
using Distributed
@everywhere using GeneralizedSasakiNakamura
solutions = pmap(p -> GSN_radial(p..., IN), parameters)
```

Configure workers through the cluster launcher or `addprocs`; each worker has
its own package state. Package loading and first-call compilation are separate
from warmed solve time.
