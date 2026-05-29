# API Reference

This page documents the public interface intended for direct user calls. The package uses units with $G = c = M = 1$.

## Exported Functions

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
| `method` | `"auto"` | `"auto"` and `"ISEM"` use the ISEM homogeneous solver; `"linear"` and `"Riccati"` use the legacy GSN ODE solvers |
| `tolerance` | internal default | ODE or ISEM matching tolerance, depending on method |
| `xm`, `rhom` | `nothing` | optional ISEM matching controls |
| `N` | `nothing` | optional ISEM expansion order override |
| `sfe`, `lfe` | `nothing` | optional small/large-frequency expansion switches |
| `TSinInf`, `TSoutInf`, `TSinHor`, `TSoutHor` | `nothing` | optional Teukolsky-Starobinsky identity switches |

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

`method = "auto"` currently routes real-frequency homogeneous solutions through ISEM. Legacy `"linear"` and `"Riccati"` methods remain available for direct GSN ODE evolution.

Returned objects are `GSNRadialFunction`. They are callable:

```julia
Xin, Xup = GSN_radial(-2, 2, 2, 0.9, 0.5)
Xin(rstar_from_r(10.0, 0.9))
Xin.GSN_solution(rstar_from_r(10.0, 0.9))
```

`GSNRadialFunction(rstar)` returns only the GSN function value. `GSN_solution(rstar)` returns the stored GSN data from the selected backend.

## Point-Particle Single-Mode Fluxes

### `Teukolsky_pointparticle_mode`

```julia
mode = Teukolsky_pointparticle_mode(s, l, m, n, k, a, p, e, x; method = "auto", N = -1, K = -1)
```

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
| `"auto"` | currently selects `"isem_trapezoidal"` |
| `"isem_trapezoidal"` | ISEM radial solution with trapezoidal convolution |
| `"isem_levin"` | ISEM radial solution with Levin convolution |
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
| `Y_solution` | auxiliary radial solution used by the convolution |
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
flux = Teukolsky_pointparticle_flux(a, p, e, x; tol = 1e-8, lmax = 30, nmax = 500, kmax = 50)
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
| `kmax` | `50` | maximum polar shell index |
| `minimum_consecutive` | `2` | number of consecutive small shells required for truncation |
| `N`, `N0` | `64` | initial radial grid interval count for adaptive sampling |
| `K`, `K0` | `16` | initial polar grid interval count for adaptive sampling |
| `Nmax` | `2^14` | maximum radial grid interval count |
| `Kmax` | `2^12` | maximum polar grid interval count |
| `sample_tol` | `1e-3` | adaptive single-mode sampling tolerance |
| `record` | `0` | when set to `1`, records mode details to HDF5 where supported |
| `record_path` | `nothing` | optional HDF5 output path |
| `fast` | `true` | use cached/presampled fast summation path where available |

If the supplied orbital parameters do not define a bound orbit according to `KerrGeodesics`, the function emits a warning and returns `nothing`.

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
| `reached` | named tuple of reached shell indices, e.g. `n_reached_inf` or `k_reached_hor` |
| `cost` | wall-clock runtime in seconds |
| `result` | raw mode-summation result |

Example:

```julia
flux = Teukolsky_pointparticle_flux(0.9, 6.0, 0.7, cos(pi / 4))
flux.infinity_energy_flux
flux.horizon_energy_flux
```

## HDF5 Mode Records

When `record = 1` is used in the mode-summation interface, mode details are written to an HDF5 file. Records are stored hierarchically by shell indices and mode indices, with datasets for amplitudes, fluxes, and grid sizes. Complex quantities are stored by real and imaginary parts for portable HDF5 access.

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
| `P_solution` | internal ISEM `P` solution when available |
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
| `numerical_GSN_solution` | ODE solution for legacy methods; ISEM metadata for `method == "ISEM"` |
| `numerical_Riccati_solution` | Riccati ODE solution when applicable |
| `GSN_solution` | callable GSN radial solution |
| `normalization_convention` | normally `UNIT_GSN_TRANS` |
| `method` | solver method string |

### `TeukolskyPointParticleMode` and `GSNPointParticleMode`

These store single-mode inhomogeneous results. `TeukolskyPointParticleMode` stores the amplitude in the Teukolsky formalism; `GSNPointParticleMode` stores the amplitude in the GSN formalism. Flux fields are shared in meaning.

### `TeukolskyPointParticleFlux`

Stores total flux results from `Teukolsky_pointparticle_flux`; see [Point-Particle Total Fluxes](@ref).
