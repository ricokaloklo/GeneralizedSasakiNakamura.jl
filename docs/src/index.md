# Home

GeneralizedSasakiNakamura.jl computes solutions to the frequency-domain radial Teukolsky equation with the Generalized Sasaki-Nakamura (GSN) formalism.

The code is capable of handling *both ingoing and outgoing* radiation of scalar, electromagnetic, and gravitational type (corresponding to spin weight of $s = 0, \pm 1, \pm 2$ respectively).

The angular Teukolsky equation is solved with an accompanying julia package [SpinWeightedSpheroidalHarmonics.jl](https://github.com/ricokaloklo/SpinWeightedSpheroidalHarmonics.jl) using a spectral decomposition method.

Both codes are capable of handling *complex* frequencies, and we use $M = 1$ convention throughout.

The paper describing both the GSN formalism and the implementation can be found in [2306.16469](https://arxiv.org/abs/2306.16469). A set of Mathematica notebooks deriving all the equations used in the code can be found in [10.5281/zenodo.8080241](https://zenodo.org/records/8080242).

Starting from v0.8.0, the code is also capable of computing the gravitational waveform amplitude and fluxes at infinity and at the horizon due to a test particle orbiting around a Kerr black hole in a _generic (eccentric, inclined) timelike bound orbit_ by solving the inhomogeneous SN equation using integration by parts.

Starting from v0.9.0, the package includes the ISEM solver, short for
_iterative series expansion matching_, and a high-level total-flux interface,
[`Teukolsky_pointparticle_flux`](@ref), which automatically selects the
circular, eccentric, inclined, or generic mode-summation strategy.

Version 0.10.0 replaces it with a GSN-based ISEM solver that directly solves the generalized Sasaki–Nakamura equation using iterative series expansion matching (ISEM), along with a Kerr quasinormal-mode solver. To support large-scale parallel computations, the point-particle mode and flux solvers now support thread-safe concurrent calls on the default route, sharing synchronized caches within each submission. It requires Julia 1.12, KerrGeodesics v0.5 and SpinWeightedSpheroidalHarmonics v1.4.

For high-index tail modes in eccentric and generic flux summations, large radial indices n make the convolution integrands highly oscillatory. These integrals can now use adaptive Levin quadrature in the radial direction, combined with Clenshaw–Curtis quadrature in the polar direction for generic orbits, instead of trapezoidal quadrature. Method selection and refinement are handled automatically: users only need to specify the desired accuracy, and the solver adapts the computation to meet that target.

## Installation
Julia 1.12 or later in the 1.x series is required.

To install the package using the Julia package manager, simply type the following in the Julia REPL:
```julia
using Pkg
Pkg.add("GeneralizedSasakiNakamura")
```

*Note: [SpinWeightedSpheroidalHarmonics.jl](https://github.com/ricokaloklo/SpinWeightedSpheroidalHarmonics.jl) v1.4 or later is installed automatically from the General registry.*

## Highlights
### Numerical solver path: linear/Riccati integration
The numerical path directly evolves the radial equation with the `linear` or `Riccati` methods and stitches analytical boundary ansatzes near the horizon and infinity:
```@raw html
<table>
  <tr>
    <th>GeneralizedSasakiNakamura.jl</th>
    <th><a href="https://github.com/BlackHolePerturbationToolkit/Teukolsky">Teukolsky</a> Mathematica package using the MST method </th>
  </tr>
  <tr>
    <td><p align="center"><img width="100%" src="https://github-production-user-asset-6210df.s3.amazonaws.com/55488840/248965077-7d216deb-5bae-433f-a699-d40a35f0e35d.gif"></p></td>
    <td><p align="center"><img width="100%" src="https://github-production-user-asset-6210df.s3.amazonaws.com/55488840/248966033-9e7d8027-81ee-4762-98d9-0ad0a1c030ad.gif"></p></td>
  </tr>
</table>
```

*(There was no caching in this benchmark; the equation was solved on the fly. The notebook generating this animation can be found [here](https://github.com/ricokaloklo/GeneralizedSasakiNakamura.jl/blob/main/examples/realtime-demo.ipynb).)*

### Semi-analytical solver path: ISEM
ISEM constructs the GSN function $X$ directly, by matching series expansions of
the generalized Sasaki–Nakamura equation, and derives the Teukolsky function $R$
and the source-adapted function $Y$ from the same matched solution.
`method = "auto"` tries ISEM first and falls back to the legacy automatic radial
solver only after an explicit construction failure; `method = "ISEM"` selects it
strictly with no fallback, while `linear` and `Riccati` select the legacy ODE solvers.
In v0.10.0 this GSN-based solver replaces the Teukolsky-based ISEM of v0.9.0,
whose `rhom` and `TSinInf`/`TSoutInf`/`TSinHor`/`TSoutHor` controls and negative
`xm` coordinate are gone.

```@raw html
<p align="center">
  <img width="80%" src="isem_matching_original_30fps.gif">
</p>
```

```@raw html
<p align="center">
  <img width="80%" src="https://raw.githubusercontent.com/CuberYyc808/CuberYyc808.github.io/3b847bfb92587eecb078bf7284f1cbc86e476ddb/images/teukolsky_gsn_solver_demo_60s.gif">
</p>
```

Static/zero-frequency solutions are solved analytically with Gauss hypergeometric functions.

Superradiance-threshold solutions with $\omega = m a / (2 r_+)$ are built with a logarithmic horizon solution.

### QNM interface

The native QNM interface computes a Leaver root and validates its scattering
data with the ISEM radial route:

```julia
ordinary_mode = qnm(0.9, -2, 2, 2, 10)
mirror_mode = qnm(0.9, -2, 2, 2, 10, mirror)
```

The five numerical inputs are `(a, s, l, m, n)`. `ordinary` is the default
branch; `mirror` requests a directly polished negative-real root and directly
evaluated amplitudes. The compact result exposes `omega`, `lambda`,
`incidence_amplitude`, `reflection_amplitude`, `incidence_derivative`, and
`excitation_factor`. Use `detailed=true` only when the Teukolsky bridge and
callable `X`, `Y`, and `R` objects are needed. See [Quasinormal modes](@ref).

### Solutions that are accurate everywhere
Numerical-path solutions are *smoothly stitched* to analytical ansatzes near the horizon and infinity at user-specified locations `rsin` and `rsout` respectively:

```@raw html
<p align="center">
  <img width="50%" src="https://github-production-user-asset-6210df.s3.amazonaws.com/55488840/248724944-9707332b-1238-4b3b-b1c0-ac426a1b3dc6.gif">
</p>
```

### Easy to use
The following code snippet lets you solve the (source-free) Teukolsky function (in frequency domain) for the mode $s=-2, \ell=2, m=2, a/M=0.7, M\omega=0.5$ that satisfies the purely-ingoing boundary condition at the horizon, $R^{\textrm{in}}$, and the purely-outgoing boundary condition at spatial infinity, $R^{\textrm{up}}$, respectively:
```julia
using GeneralizedSasakiNakamura # This is going to take some time to pre-compile, mostly due to DifferentialEquations.jl

# Specify which mode to solve
s=-2; l=2; m=2; a=0.7; omega=0.5;

# NOTE: julia uses 'just-ahead-of-time' compilation. Calling this the first time in each session will take some time
Rin, Rup = Teukolsky_radial(s, l, m, a, omega)
```
That's it! If you run this on Julia REPL, it should give you something like this
```
(
TeukolskyRadialFunction(
    mode = Mode(s = -2, l = 2, m = 2, a = 0.7, omega = 0.5, lambda = 1.6966094016353415),
    boundary_condition = IN,
    transmission_amplitude = 1.0 + 0.0im,
    incidence_amplitude = 6.536587661184236 - 4.941203897065118im,
    reflection_amplitude = -0.1282466191290958 - 0.4404813349644446im,
    normalization_convention = UNIT_TEUKOLSKY_TRANS
),
TeukolskyRadialFunction(
    mode = Mode(s = -2, l = 2, m = 2, a = 0.7, omega = 0.5, lambda = 1.6966094016353415),
    boundary_condition = UP,
    transmission_amplitude = 1.0 + 0.0im,
    incidence_amplitude = -1.1698840333865992 - 2.5455723340434875im,
    reflection_amplitude = 2.5169908586317717 - 8.644964686261376im,
    normalization_convention = UNIT_TEUKOLSKY_TRANS
))
```
In Julia REPL, you can check out all the asymptotic amplitudes at a glimpse using something like
```julia
julia> Rin
TeukolskyRadialFunction(
    mode = Mode(s = -2, l = 2, m = 2, a = 0.7, omega = 0.5, lambda = 1.6966094016353415),
    boundary_condition = IN,
    transmission_amplitude = 1.0 + 0.0im,
    incidence_amplitude = 6.536587661184236 - 4.941203897065118im,
    reflection_amplitude = -0.1282466191290958 - 0.4404813349644446im,
    normalization_convention = UNIT_TEUKOLSKY_TRANS
)
```

For example, if we want to evaluate the Teukolsky function $R^{\textrm{in}}$ at the location $r = 10M$, simply do
```julia
Rin(10)
```
This should give
```
77.57508416830544 - 429.4029095224015im
```

#### Solving radial equations at complex frequencies
The radial interfaces also accept a supplied complex frequency. For example,
the $s=-2, \ell=2, m=2, a/M=0.68$ fundamental QNM frequency can be checked
independently with
```julia
Rin, Rup = Teukolsky_radial(-2, 2, 2, 0.68, 0.5239751-0.0815126im)
```
We can check out the $R^{\textrm{up}}$ solution using
```julia
julia> Rup
TeukolskyRadialFunction(
    mode = Mode(s = -2, l = 2, m = 2, a = 0.68, omega = 0.5239751 - 0.0815126im, lambda = 1.6550030805786913 + 0.36026765638858793im),
    boundary_condition = UP,
    transmission_amplitude = 1.0 + 0.0im,
    incidence_amplitude = -5.8212734659171046e-8 - 3.805329485211967e-7im,
    reflection_amplitude = 1.1011632105029665 + 2.1300597368117073im,
    normalization_convention = UNIT_TEUKOLSKY_TRANS
)
```
We see that the incidence amplitude is indeed very small numerically as a QNM solution should. This can be accessed using
```julia
Rup.incidence_amplitude
```

This should give
```julia
-5.8212734659171046e-8 - 3.805329485211967e-7im
```

#### Solving the inhomogeneous radial Teukolsky/SN equation with a point-particle source on a generic timelike bound orbit
This can now be done easily with this code.
Suppose we want to compute the inhomogeneous solution to the radial Teukolsky equation at infinity for the $s = -2$, $\ell = m = 2$ mode driven by a test particle on a bound geodesic with $a/M = 0.9, p = 6M, e = 0.7, x = \cos(\pi/4)$, one can simply do
```julia
mode_info = Teukolsky_pointparticle_mode(-2, 2, 2, 0, 0, 0.9, 6, 0.7, cos(π/4))
```
where $n = 0$ and $k = 0$ label the radial and polar modes, respectively.
To have a glimpse of the output, one can do so with
```julia
julia> mode_info
TeukolskyPointParticleMode(
    mode = Mode(s = -2, l = 2, m = 2, n = 0, k = 0, a = 0.9, omega = 0.06568724726732753, lambda = 3.606789012119982),
    amplitude_inf = 0.00023429507956768906 - 6.558414441160938e-5im,
    energy_flux_inf = 1.0917330113048546e-6,
    angular_momentum_flux_inf = 3.324033375494718e-5,
    Carter_const_flux_inf = 5.890504443058139e-5,
    method = (method = "isem_trapezoidal", radial_method = "ISEM", radial_sfe = false, N = 256, K = 64, truncation_floor = 1.0e-16),
)
```
To access for example the amplitude at infinity
```julia
julia> mode_info.amplitude
0.00023429507956768906 - 6.558414441160938e-5im
```
which is the value for $Z^{\infty}_{\ell m n k}$, the amplitude of the inhomogeneous radial Teukolsky solution near infinity for that particular frequency.

If we want to compute the inhomogeneous solution to the radial Teukolsky equation at the event horizon for the 
same set of parameters, we can simply change the sign of $s$ to $2$
```julia
mode_info = Teukolsky_pointparticle_mode(2, 2, 2, 0, 0, 0.9, 6, 0.7, cos(π/4))
```
The output should be
```julia
julia> mode_info
TeukolskyPointParticleMode(
    mode = Mode(s = 2, l = 2, m = 2, n = 0, k = 0, a = 0.9, omega = 0.06568724726732753, lambda = -0.39321098788001757),
    amplitude_hor = 0.006089946888790668 - 0.0014130019665083518im,
    energy_flux_hor = -2.843814878429853e-9,
    angular_momentum_flux_hor = -8.658651402627282e-8,
    Carter_const_flux_hor = -1.5343956812851395e-7,
    method = (method = "isem_trapezoidal", radial_method = "ISEM", radial_sfe = false, N = 256, K = 64, truncation_floor = 1.0e-16),
)
```
To access for example the amplitude at the horizon
```julia
julia> mode_info.amplitude
0.006089946888790668 - 0.0014130019665083518im
```
which is the value for $Z^{\mathrm{H}}_{\ell m n k}$, the amplitude of the inhomogeneous radial Teukolsky solution near the horizon for that particular frequency.

Total fluxes can be computed with the orbit-aware high-level interface. A generic-orbit total-flux run can be substantially slower than a single-mode calculation because it performs a full mode summation.
```julia
flux = Teukolsky_pointparticle_flux(0.9, 6.0, 0.7, 0.5)
flux.infinity_energy_flux
flux.horizon_energy_flux
```

The function automatically dispatches to circular, eccentric, inclined, or generic mode summation according to the supplied orbital parameters.

For modes with large `n`, eccentric and generic summations can switch from uniform-grid trapezoidal sampling to adaptive Levin quadrature. For generic two-dimensional convolutions, the default accelerated tail path uses adaptive Levin in the radial direction and Clenshaw-Curtis sampling in the polar direction, reducing the need for a uniformly dense two-dimensional grid.

!!! details "Circular equatorial orbit"
    ```julia
    julia> flux = Teukolsky_pointparticle_flux(0.9, 6.0, 0.0, 1.0; tol=1e-8)
    TeukolskyPointParticleFlux(
        orbital_parameters(a = 0.9, p = 6.0, e = 0.0, x = 1.0),
        orbit_type = circular,
        infinity_energy_flux = 0.0005658659548571375,
        infinity_angular_momentum_flux = 0.008825776472648025,
        infinity_carter_constant_flux = 0.0,
        horizon_energy_flux = -4.177363290666302e-6,
        horizon_angular_momentum_flux = -6.515407815579624e-5,
        horizon_carter_constant_flux = 0.0,
        total_modes = 378,
        l_reached = (infinity = 19, horizon = 12),
        convolution_integral = (strategy = "no convolution integral is needed",),
        tolerance = 1.0e-8,
        truncation_floor = (infinity = 1.0e-16, horizon = 1.0e-16),
        cost = 0.6137290000915527 seconds,
    )
    ```

    With `e = 0` and `x = 1` the call dispatches to the circular equatorial summation, which needs no convolution integral. In this example the warm run averaged about `1.79 ms` per computed mode.

!!! details "Spherical (inclined circular) orbit"
    ```julia
    julia> flux = Teukolsky_pointparticle_flux(0.9, 6.0, 0.0, 0.5; tol=1e-8)
    TeukolskyPointParticleFlux(
        orbital_parameters(a = 0.9, p = 6.0, e = 0.0, x = 0.5),
        orbit_type = inclined,
        infinity_energy_flux = 0.0006637621566020543,
        infinity_angular_momentum_flux = 0.0060183149433979885,
        infinity_carter_constant_flux = 0.04294455189603736,
        horizon_energy_flux = -3.0865749224058106e-6,
        horizon_angular_momentum_flux = -0.00010495330048599927,
        horizon_carter_constant_flux = 8.871154473182051e-5,
        total_modes = 11466,
        k_reached_inf = 15,
        k_reached_hor = 13,
        convolution_integral = (strategy = "ISEM adaptive trapezoidal for all modes",),
        tolerance = 1.0e-8,
        truncation_floor = (infinity = 1.0e-16, horizon = 1.0e-16),
        cost = 10.54723310470581 seconds,
    )
    ```

    With `e = 0` and `|x| < 1` the call dispatches to the inclined (spherical-orbit) summation. In this example the warm run averaged about `0.66 ms` per computed mode.

After an API warmup call,
```julia
Teukolsky_pointparticle_flux(0.9, 6.0, 0.01, 0.99)
```
the example above produced the following result on the test machine:

```julia
julia> flux = Teukolsky_pointparticle_flux(0.9, 6.0, 0.7, 0.5; tol = 1e-8)
TeukolskyPointParticleFlux(
    orbital_parameters(a = 0.9, p = 6.0, e = 0.7, x = 0.5),
    orbit_type = generic,
    infinity_energy_flux = 0.001146913075693756,
    infinity_angular_momentum_flux = 0.006380279786756462,
    infinity_carter_constant_flux = 0.041463294475529405,
    horizon_energy_flux = -9.20041236363546e-6,
    horizon_angular_momentum_flux = -0.0003627225703136479,
    horizon_carter_constant_flux = 0.0009653818407531111,
    total_modes = 1146074,
    n_reached = (infinity = 145, horizon = 84),
    convolution_integral = (strategy = "ISEM adaptive trapezoidal for all computed n; tail ISEM adaptive Levin enabled but not triggered", tail_levin_nmin = 50, tail_levin_local_n = 16, tail_levin_max_depth = 8),
    tolerance = 1.0e-8,
    truncation_floor = (infinity = 1.0e-16, horizon = 1.0e-16),
    cost = 1269.3770909309387 seconds,
)
```


## How to cite
If you have used this code in your research that leads to a publication, please cite the following article:
```
@article{Lo:2023fvv,
    author = "Lo, Rico K. L.",
    title = "{Recipes for computing radiation from a Kerr black hole using a generalized Sasaki-Nakamura formalism: Homogeneous solutions}",
    eprint = "2306.16469",
    archivePrefix = "arXiv",
    primaryClass = "gr-qc",
    doi = "10.1103/PhysRevD.110.124070",
    journal = "Phys. Rev. D",
    volume = "110",
    number = "12",
    pages = "124070",
    year = "2024"
}
```

Additionally, if you have used this code's capability to solve for solutions with complex frequencies, please also cite the following article:
```
@article{Lo:2025njp,
    author = "Lo, Rico K. L. and Sabani, Leart and Cardoso, Vitor",
    title = "{Quasinormal modes and excitation factors of Kerr black holes}",
    eprint = "2504.00084",
    archivePrefix = "arXiv",
    primaryClass = "gr-qc",
    doi = "10.1103/PhysRevD.111.124002",
    journal = "Phys. Rev. D",
    volume = "111",
    number = "12",
    pages = "124002",
    year = "2025"
}
```

If you have used this code's capability to solve for the gravitational waveform amplitudes and fluxes at infinity and at horizon with a test particle orbiting a Kerr black hole in a generic timelike bound and stable orbit (e.g., for extreme mass ratio inspiral waveforms), please cite the following articles:
```
@article{Yin:2025kls,
    author = "Yin, Yucheng and Lo, Rico K. L. and Chen, Xian",
    title = "{Gravitational radiation from Kerr black holes using the Sasaki-Nakamura formalism: waveforms and fluxes at infinity}",
    eprint = "2511.08673",
    archivePrefix = "arXiv",
    primaryClass = "gr-qc",
    month = "11",
    year = "2025"
}

@article{Lo:2025lpo,
    author = "Lo, Rico K. L. and Yin, Yucheng",
    title = "{Near-horizon gravitational perturbations of rotating black holes}",
    eprint = "2512.07937",
    archivePrefix = "arXiv",
    primaryClass = "gr-qc",
    doi = "10.1103/bljh-l413",
    journal = "Phys. Rev. D",
    volume = "113",
    number = "6",
    pages = "L061505",
    year = "2026"
}
```

## License
The package is licensed under the MIT License.
