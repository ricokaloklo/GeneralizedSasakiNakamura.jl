# GeneralizedSasakiNakamura.jl

![license](https://img.shields.io/github/license/ricokaloklo/GeneralizedSasakiNakamura.jl)
[![GitHub release](https://img.shields.io/github/v/release/ricokaloklo/GeneralizedSasakiNakamura.jl.svg)](https://github.com/ricokaloklo/GeneralizedSasakiNakamura.jl/releases)
[![Documentation](https://img.shields.io/badge/Documentation-ready)](http://ricokaloklo.github.io/GeneralizedSasakiNakamura.jl)

GeneralizedSasakiNakamura.jl computes solutions to the frequency-domain radial Teukolsky equation with the Generalized Sasaki-Nakamura (GSN) formalism.

The code is capable of handling *both in-going and out-going* radiation of scalar, electromagnetic, and gravitational type (corresponding to spin weight of $s = 0, \pm 1, \pm 2$ respectively).

The angular Teukolsky equation is solved with an accompanying julia package [SpinWeightedSpheroidalHarmonics.jl](https://github.com/ricokaloklo/SpinWeightedSpheroidalHarmonics.jl) using a spectral decomposition method.

## Installation
To install the package using the Julia package manager, simply type the following in the Julia REPL:
```julia
using Pkg
Pkg.add("GeneralizedSasakiNakamura")
```

*Note: There is no need to install [SpinWeightedSpheroidalHarmonics.jl](https://github.com/ricokaloklo/SpinWeightedSpheroidalHarmonics.jl) separately as it should be automatically installed by the package manager.*

## Highlights
### Performant frequency-domain Teukolsky solver
Works well at *both low and high frequencies*, and takes only a few tens of milliseconds on average:
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

*(There was no caching! We solved the equation on-the-fly! The notebook generating this animation can be found [here](https://github.com/ricokaloklo/GeneralizedSasakiNakamura.jl/blob/main/examples/realtime-demo.ipynb))*

Static/zero-frequency solutions are solved analytically with Gauss hypergeometric functions.

### Solutions that are accurate everywhere
Numerical solutions are *smoothly stitched* to analytical ansatzes near the horizon and infinity at user-specified locations `rsin` and `rsout` respectively:

<p align="center">
  <img width="50%" src="https://github-production-user-asset-6210df.s3.amazonaws.com/55488840/248724944-9707332b-1238-4b3b-b1c0-ac426a1b3dc6.gif">
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
For example, if we want to evaluate the Teukolsky function at the location $r = 10M$, simply do
```julia
R(10)
```
This should give
```
77.57508416832009 - 429.40290952257226im
```

## How to cite
If you have used this code in your research that leads to a publication, please cite the following article:
```
@article{Lo:2023fvv,
    author = "Lo, Rico K. L.",
    title = "{Recipes for computing radiation from a Kerr black hole using Generalized Sasaki-Nakamura formalism: I. Homogeneous solutions}",
    eprint = "2306.16469",
    archivePrefix = "arXiv",
    primaryClass = "gr-qc",
    month = "6",
    year = "2023"
}
```

## License
The package is licensed under the MIT License.
