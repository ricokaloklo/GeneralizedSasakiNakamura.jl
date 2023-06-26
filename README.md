# GeneralizedSasakiNakamura.jl

![license](https://img.shields.io/github/license/ricokaloklo/GeneralizedSasakiNakamura.jl)
[![GitHub release](https://img.shields.io/github/v/release/ricokaloklo/GeneralizedSasakiNakamura.jl.svg)](https://github.com/ricokaloklo/GeneralizedSasakiNakamura.jl/releases)
[![Documentation](https://img.shields.io/badge/Documentation-ready)](http://ricokaloklo.github.io/GeneralizedSasakiNakamura.jl)

GeneralizedSasakiNakamura.jl computes solutions to the frequency-domain radial Teukolsky equation with the Generalized Sasaki-Nakamura (GSN) formalism.

The angular Teukolsky equation is solved with another julia package [SpinWeightedSpheroidalHarmonics.jl](https://github.com/ricokaloklo/SpinWeightedSpheroidalHarmonics.jl) using a spectral decomposition method.

## Installation
To install the package using the Julia package manager, simply type the following in the Julia REPL:
```julia
using Pkg
Pkg.add("GeneralizedSasakiNakamura")
```

*Note: There is no need to install [SpinWeightedSpheroidalHarmonics.jl](https://github.com/ricokaloklo/SpinWeightedSpheroidalHarmonics.jl) manually as it should be automatically installed by the package manager.*

## Highlights
### Performant frequency-domain Teukolsky solver
Takes on average only a few tens of milliseconds:
<p align="center">
  <img width="75%" src="https://github-production-user-asset-6210df.s3.amazonaws.com/55488840/248602615-eb611e3f-00b9-49c3-905c-0440b0607683.gif">
</p>

*(There was no caching! We solved the equation on-the-fly! The notebook generating this animation can be found [here](https://github.com/ricokaloklo/GeneralizedSasakiNakamura.jl/blob/main/examples/realtime-demo.ipynb))*

### Solutions that are accurate everywhere
Numerical solutions are *smoothly stitched* to analytical ansatzes near the horizon and infinity at user-specified locations `rsin` and `rsout` respectively:

<p align="center">
  <img width="75%" src="https://github-production-user-asset-6210df.s3.amazonaws.com/55488840/248601145-8ce9b236-d197-4321-83b9-fcb9692a51e8.gif">
</p>

### Easy to use
The following code snippet lets you solve the (source-free) Teukolsky function (in frequency domain) for the mode $s=-2, \ell=2, m=2, a=0.7, \omega=0.5$ that satisfies the purely-ingoing boundary condition at the horizon:
```julia
using GeneralizedSasakiNakamura # This is going to take some time to pre-compile, mostly due to DifferentialEquations.jl

# Specify which mode and what boundary condition
s=-2; l=2; m=2; a=0.7; omega=0.5; bc=IN;
# Specify where to match to ansatzes
rsin=-20; rsout=250;

# NOTE: julia uses 'just-ahead-of-time' compilation. Calling this the first time in each session will take some time
R = Teukolsky_radial(s, l, m, a, omega, bc, rsin, rsout) 
```
That's it! If you run this on Julia REPL, it should give you something like this
```
TeukolskyRadialFunction(
    mode=Mode(s=-2, l=2, m=2, a=0.7, omega=0.5, lambda=1.6966094016353415),
    boundary_condition=IN,
    transmission_amplitude=1.0 + 0.0im,
    incidence_amplitude=6.536587661197995 - 4.941203897068852im,
    reflection_amplitude=-0.128246619129379 - 0.44048133496664404im,
    normalization_convention=UNIT_TEUKOLSKY_TRANS
)
```

## License
The package is licensed under the MIT License.
