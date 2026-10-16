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
# keywords: convention=:overtone | :complete_spectrum | :leaver, segment=0 | 1 | :auto
```

`ordinary` is the default branch. `mirror` requests the mode with negative real
frequency and its amplitudes. The compact result contains:

| field | meaning |
| :--- | :--- |
| `omega` or `frequency` | QNM frequency (root of Leaver's continued fraction) |
| `lambda` | radial Teukolsky separation constant |
| `incidence_amplitude` | GSN incoming amplitude |
| `reflection_amplitude` | GSN reflected amplitude |
| `incidence_derivative` | $dA_{\rm in}/d\omega$ at the root |
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
| `method` | `"auto"` | `"auto"` tries ISEM first and falls back to the legacy automatic route only on an ISEM failure with no explicit `xm`; `"ISEM"` is strict; `"linear"` and `"Riccati"` select the legacy GSN ODE solvers |
| `tol` | internal default | ISEM matching tolerance, or the tolerance of the ODE methods |
| `tolerance` | internal default | tolerance of the ODE methods; `tol` takes precedence when both are supplied. Tuple form and the form with explicit `rsin, rsout`; not the boundary-specific form |
| `xm` | `nothing` | ISEM matching point in the compactified radial coordinate `x ∈ (0, 1)`, which maps the horizon to `0` and infinity to `1` |
| `N` | `nothing` | optional ISEM expansion order override; chosen automatically when omitted |
| `sfe`, `lfe` | `nothing` | switches for the small-frequency expansion (MST, the Mano–Suzuki–Takasugi series) and the large-frequency (eikonal) expansion; `nothing` selects them automatically |

Omit optional matching controls to use automatic selection.
If construction at an explicit `xm` fails, `"auto"` reports the failure rather
than falling back to a legacy solver that ignores that junction.
Automatic matching is recommended: an explicit `xm` can introduce a
cancellation error in the incidence amplitude that the matching tolerance does
not bound.
The v0.9.0 ISEM keywords `rhom`, `TSinInf`, `TSoutInf`, `TSinHor`, `TSoutHor`
and the negative `xm` coordinate no longer exist; passing them raises an error.

Only `omega == 0` uses the analytic static branch. Small nonzero frequencies
retain their finite-frequency boundary conditions. At the superradiant threshold
`omega = m a / (2 r_+)` the UP solution is always built with ISEM, using a
logarithmic horizon solution; its incidence and reflection amplitudes are then the
coefficients of that logarithmic horizon basis. Nearby frequencies are not
substituted by the threshold value. Likewise, exact-extremal
dispatch requires `a == 1` or `a == -1`.

Returned objects are `TeukolskyRadialFunction`. They are callable:

```julia
Rin, Rup = Teukolsky_radial(-2, 2, 2, 0.9, 0.5)
Rin(10.0)
Rin.Teukolsky_solution(10.0)
```

`TeukolskyRadialFunction(r)` returns the value of `R(r)`. `Teukolsky_solution(r)` returns `R` and `dR/dr` as its first two entries; with `method = "ISEM"` it returns the tuple `(R, dR/dr, d²R/dr², err)`.

### `GSN_radial`

Main call forms:

```julia
Xin, Xup = GSN_radial(s, l, m, a, omega; method = "auto")
X = GSN_radial(s, l, m, a, omega, IN; method = "auto")
X = GSN_radial(s, l, m, a, omega, UP; method = "auto")
```

`method = "auto"` tries ISEM first and falls back to the legacy automatic radial route only after an ISEM failure with no explicit `xm`. `method = "ISEM"` is strict and does not fall back, and `"linear"` and `"Riccati"` select the legacy GSN ODE evolution. Matching controls, the static branch, and the superradiance-threshold branch follow the same rules as `Teukolsky_radial`.

Returned objects are `GSNRadialFunction`. They are callable:

```julia
Xin, Xup = GSN_radial(-2, 2, 2, 0.9, 0.5)
Xin(rstar_from_r(0.9, 10.0))
Xin.GSN_solution(rstar_from_r(0.9, 10.0))
```

`GSNRadialFunction(rstar)` returns the value of `X(r_*)`. With ISEM, `GSN_solution(rstar)` returns the tuple `(X, dX/dr_*, err)`, where `err` is an estimate of the local error, not a bound.

Keyword summary (in addition to the matching controls shared with `Teukolsky_radial`):

| keyword | default | meaning |
| :--- | :--- | :--- |
| `tol` | `nothing` | tolerance; alias for `tolerance` that takes precedence when both are supplied |
| `tolerance` | `nothing` (`1e-12` for the ODE methods) | tolerance of the ODE methods |
| `lambda` | `nothing` | precomputed spin-weighted spheroidal eigenvalue; computed from `s, l, m, a*omega` when omitted. Only the form with explicit `rsin, rsout` |
| `ODE_algorithm` | `AutoVern9(Rosenbrock23(autodiff=false))` | ODE solver of the `"linear"` and `"Riccati"` methods. Tuple form and the form with explicit `rsin, rsout` |
| `data_type` | `ComplexF64` | number type of the ODE methods. Tuple form and the form with explicit `rsin, rsout` |

### `Y_radial`

```julia
Y = Y_radial(s, l, m, a, omega, IN; method = "auto")
```

Constructs the auxiliary radial function `Y` that the point-particle convolutions integrate against the source, with the same matching controls as `Teukolsky_radial`. It is available for `s = -2` with `IN` and `s = +2` with `UP`, for nonzero `omega`; there is no tuple form. The returned `YRadialFunction` is callable.

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
| `levin_local_n` | `16` | number of local nodes per adaptive Levin subinterval |
| `levin_max_depth` | `8` | maximum adaptive Levin subdivision depth |
| `info` | `false` | no effect; kept for compatibility |
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
| `reached` | named tuple of the highest shell indices reached: `l_reached`, `l_reached_inf`, `l_reached_hor` (circular), `k_reached_inf`, `k_reached_hor` (inclined) or `n_reached_inf`, `n_reached_hor` (eccentric and generic); the display groups the infinity and horizon values |
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

A generic orbit such as `(0.9, 6.0, 0.7, 0.5)` sums over about a million modes
and takes about 20 minutes at `tol = 1e-8`.

## HDF5 Mode Records

When `record = true` is used in the mode-summation interface, mode details are written to an HDF5 file. Records are stored hierarchically by shell indices and mode indices, with datasets for amplitudes, fluxes, and grid sizes. Complex quantities are stored by real and imaginary parts for portable HDF5 access.

Without `record_path` the file name is built from the orbit, e.g.
`eccentric_mode_data_a_..._p_..._e_....h5`. Use `HDF5.jl` to inspect the file:

```julia
using HDF5

flux = Teukolsky_pointparticle_flux(0.9, 6.0, 0.3, 1.0; record = true,
    record_path = "eccentric_mode_data.h5")
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
| `P_solution` | internal `P`-equation solution when available (complex frequencies); an internal callable |
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
| `numerical_GSN_solution` | ODE solution for legacy methods; matching metadata for `method == "ISEM"` |
| `numerical_Riccati_solution` | Riccati ODE solution when applicable |
| `GSN_solution` | callable GSN radial solution |
| `normalization_convention` | normally `UNIT_GSN_TRANS` |
| `method` | solver method string |

#### Accuracy

The error entry returned by `GSN_solution` is an estimate, not a bound. It does not include the error of the input angular
eigenvalue or of the equation coefficients near the horizon (up to about
`1e-11`). Reflection amplitudes that are small compared with the dominant
amplitude, and solutions at strongly damped frequencies, can lose relative
accuracy; at an incidence zero, judge the incidence amplitude by its size
relative to the transmission amplitude.

The ODE methods (`"linear"`, `"Riccati"`) have their own low-frequency limits:
with the default outer boundary their incidence amplitudes can be wrong at low
frequency, and below `|omega| = 1e-8` they lose the `UP` solution;
`method = "auto"` warns when it falls back to them there.

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

Radial construction, complex-frequency evaluation, QNM frequencies and
excitation factors, and concurrent evaluation of one shared radial solution are
thread-safe and give identical results in serial and threaded runs.
Point-particle modes and fluxes also depend on the angular and geodesic
packages. KerrGeodesics 0.5 evaluates its Jacobi elliptic
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

```@docs
with_pointparticle_submission
spawn_pointparticle_submission!
cancel_pointparticle_submission!
```

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
