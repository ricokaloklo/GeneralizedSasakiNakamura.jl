# APIs

There are 6 functions that are exported, namely

- [`Teukolsky_radial`](@ref)
- [`GSN_radial`](@ref)
- [`Teukolsky_pointparticle_mode`](@ref)
- [`GSN_pointparticle_mode`] (@ref)
- [`rstar_from_r`](@ref)
- [`r_from_rstar`](@ref)

and there are 7 custom types that are exported, i.e.

- [BoundaryCondition](@ref)
- [NormalizationConvention](@ref)
- [Mode](@ref)
- [GSNRadialFunction](@ref)
- [TeukolskyRadialFunction](@ref)
- [GSNPointParticleMode](@ref)
- [TeukolskyPointParticleMode](@ref)

Currently, only the exported functions and types are documented below. Documentations for private (i.e. unexported) functions will be added at a later stage.

## Functions
```@docs
Teukolsky_radial
```

```@docs
GSN_radial
```

```@docs
rstar_from_r
```

```@docs
r_from_rstar
```

```@docs
Teukolsky_pointparticle_mode
```

```@docs
GSN_radial
```

## Types
#### BoundaryCondition
This is an enum type that can take either one of the four values

| value |  |
| :--- | :--- |
| `IN` | purely ingoing at the horizon | 
| `UP` | purely outgoing at infinity |
| `OUT` | purely outgoing at the horizon |
| `DOWN`| purely ingoing at infinity |

#### NormalizationConvention
This is an enum type that can take either one of the two values

| value |   |
| :--- | :--- |
| `UNIT_GSN_TRANS` | normalized to have a unit transmission amplitude for the GSN function | 
| `UNIT_TEUKOLSKY_TRANS` | normalized to have a unit transmission amplitude for the Teukolsky function |

#### Mode
This is a composite struct type that stores information about a mode

| field |   |
| :--- | :--- |
| `s` | spin weight $s$ |
| `l` | harmonic index $\ell$ |
| `m` | azimuthal index $m$ |
| `a` | Kerr spin parameter $a/M$ |
| `omega` | frequency $M\omega$ |
| `lambda` | spin-weighted spheroidal eigenvalue $\lambda$ |

#### GSNRadialFunction
This is a composite struct type that stores the output from [`GSN_radial`](@ref)

!!! tip

    `GSNRadialFunction(rstar)` is equivalent to `GSNRadialFunction.GSN_solution(rstar)[1]`, 
    returning only the value of the GSN function evaluated at the *tortoise coordinate* `rstar`

| field |    |
| :--- | :--- |
| `mode` | a [Mode](@ref) object storing information about the mode |
| `boundary_condition` | a [BoundaryCondition](@ref) object storing which boundary condition this function satisfies |
| `rsin` | numerical inner boundary $r_{*}^{\mathrm{in}}/M$ where the GSN equation is numerically evolved ($r_{*}$ is a tortoise coordinate) |
| `rsout` | numerical outer boundary $r_{*}^{\mathrm{out}}/M$ where the GSN equation is numerically evolved ($r_{*}$ is a tortoise coordinate) |
| `rsmp` | The matching point in tortoise coordinate $r_{*}^{\mathrm{mp}}$ if used  |
| `horizon_expansion_order` | order of the asymptotic expansion at the horizon |
| `infinity_expansion_order` | order of the asymptotic expansion at infinity |
| `transmission_amplitude` | transmission amplitude in the GSN formalism of this function |
| `incidence_amplitude` | incidence amplitude in the GSN formalism of this function |
| `reflection_amplitude` | reflection amplitude in the GSN formalism of this function |
| `numerical_GSN_solution` | numerical solution ([`ODESolution`](https://docs.sciml.ai/DiffEqDocs/stable/types/ode_types/#SciMLBase.ODESolution) object from `DifferentialEquations.jl`) to the GSN equation in [`rsin`, `rsout`] if applicable; output is a vector $[ \hat{X}(r_{*}), d\hat{X}(r_{*})/dr_{*} ]$ |
| `numerical_Riccati_solution` | numerical solution ([`ODESolution`](https://docs.sciml.ai/DiffEqDocs/stable/types/ode_types/#SciMLBase.ODESolution) object from `DifferentialEquations.jl`) to the GSN equation in the Riccati form if applicable; output is a vector $[ \hat{\Phi}(r_{*}), d\hat{\Phi}(r_{*})/dr_{*} ]$ |
| `GSN_solution` | full GSN solution where asymptotic solutions are smoothly attached; output is a vector $[ \hat{X}(r_{*}), d\hat{X}(r_{*})/dr_{*} ]$ |
| `normalization_convention` | a [NormalizationConvention](@ref) object storing which normalization convention this function adheres to |

#### TeukolskyRadialFunction
This is a composite struct type that stores the output from [`Teukolsky_radial`](@ref)

!!! tip

    `TeukolskyRadialFunction(r)` is equivalent to `TeukolskyRadialFunction.Teukolsky_solution(r)[1]`, 
    returning only the value of the Teukolsky function evaluated at the *Boyer-Lindquist coordinate* `r`

| field |    |
| :--- | :--- |
| `mode` | a [Mode](@ref) object storing information about the mode |
| `boundary_condition` | a [BoundaryCondition](@ref) object storing which boundary condition this function satisfies |
| `transmission_amplitude` | transmission amplitude in the Teukolsky formalism of this function |
| `incidence_amplitude` | incidence amplitude in the Teukolsky formalism of this function |
| `reflection_amplitude` | reflection amplitude in the Teukolsky formalism of this function |
| `GSN_solution` | a [GSNRadialFunction](@ref) object storing the corresponding GSN function
| `Teukolsky_solution` | Teukolsky solution where asymptotic solutions are smoothly attached; output is a vector $[ \hat{R}(r), d\hat{R}(r)/dr ]$ |
| `normalization_convention` | a [NormalizationConvention](@ref) object storing which normalization convention this function adheres to |

#### GSNPointParticleMode
This is a composite struct type that stores the output from [`GSN_pointparticle_mode`](@ref)

| field |    |
| :--- | :--- |
| `mode` | a [Mode](@ref) object storing information about the mode (including the frequency) |
| `amplitude_inf` | the amplitude of the inhomogeneous GSN solution at infinity such that $X^{\infty}_{\ell m \omega} = $ `amplitude_inf` $ \times e^{i \omega r_*}$ |
| `energy_flux_inf` | the energy flux emitted towards infinity |
| `angular_momentum_flux_inf ` | the angular momentum flux towards infinity |
| `Carter_const_flux_inf` | the Carter constant flux towards infinity |
| `trajectory` | the trajectory/geodesic that the particle follows |
| `Y_solution` | the auxiliary function $Y$ used in the calculation |
| `SWSH` | the spin-weighted spheroidal harmonic` used in the calculation |
| `method` | a named tuple storing the method and grid size used in computing the convolution integral |

#### TeukolskyPointParticleMode
This is a composite struct type that stores the output from [`Teukolsky_pointparticle_mode`](@ref)

| field |    |
| :--- | :--- |
| `mode` | a [Mode](@ref) object storing information about the mode (including the frequency) |
| `amplitude_inf` | the amplitude of the inhomogeneous Teukolsky solution at infinity such that $R^{\infty}_{\ell m \omega} = $ `amplitude_inf` $ \times r^{-(2s+1)} e^{i \omega r_*}$ |
| `energy_flux_inf` | the energy flux emitted towards infinity |
| `angular_momentum_flux_inf ` | the angular momentum flux towards infinity |
| `Carter_const_flux_inf` | the Carter constant flux towards infinity |
| `trajectory` | the trajectory/geodesic that the particle follows |
| `Y_solution` | the auxiliary function $Y$ used in the calculation |
| `SWSH` | the spin-weighted spheroidal harmonic` used in the calculation |
| `method` | a named tuple storing the method and grid size used in computing the convolution integral |
