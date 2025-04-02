# GeneralizedSasakiNakamura.jl

![license](https://img.shields.io/github/license/ricokaloklo/GeneralizedSasakiNakamura.jl)
[![GitHub release](https://img.shields.io/github/v/release/ricokaloklo/GeneralizedSasakiNakamura.jl.svg)](https://github.com/ricokaloklo/GeneralizedSasakiNakamura.jl/releases)
[![Documentation](https://img.shields.io/badge/Documentation-ready)](http://ricokaloklo.github.io/GeneralizedSasakiNakamura.jl)

GeneralizedSasakiNakamura.jl computes solutions to the frequency-domain radial Teukolsky equation with the Generalized Sasaki-Nakamura (GSN) formalism.

The code is capable of handling *both ingoing and outgoing* radiation of scalar, electromagnetic, and gravitational type (corresponding to spin weight of $s = 0, \pm 1, \pm 2$ respectively).

The angular Teukolsky equation is solved with an accompanying julia package [SpinWeightedSpheroidalHarmonics.jl](https://github.com/ricokaloklo/SpinWeightedSpheroidalHarmonics.jl) using a spectral decomposition method.

Both codes are capable of handling *complex* frequencies, and we use $M = 1$ convention throughout.

The paper describing both the GSN formalism and the implementation can be found in [2306.16469](https://arxiv.org/abs/2306.16469). A set of Mathematica notebooks deriving all the equations used in the code can be found in [10.5281/zenodo.8080241](https://zenodo.org/records/8080242).

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
(TeukolskyRadialFunction(mode=Mode(s=-2, l=2, m=2, a=0.7, omega=0.5, lambda=1.696609401635342), boundary_condition=IN), TeukolskyRadialFunction(mode=Mode(s=-2, l=2, m=2, a=0.7, omega=0.5, lambda=1.696609401635342), boundary_condition=UP))
```
In Julia REPL, you can check out all the asymptotic amplitudes at a glimpse using something like
```julia
julia> Rin
TeukolskyRadialFunction(
    mode=Mode(s=-2, l=2, m=2, a=0.7, omega=0.5, lambda=1.696609401635342),
    boundary_condition=IN,
    transmission_amplitude=1.0 + 0.0im,
    incidence_amplitude=6.5365876612287765 - 4.9412038970871555im,
    reflection_amplitude=-0.1282466191307726 - 0.440481334972911im,
    normalization_convention=UNIT_TEUKOLSKY_TRANS
)
```

For example, if we want to evaluate the Teukolsky function $R^{\textrm{in}}$ at the location $r = 10M$, simply do
```julia
Rin(10)
```
This should give
```
77.57508416835319 - 429.40290952262677im
```

#### Solving for complex frequencies
One can use the same interface to compute solutions with complex frequencies. For example, the QNM solution of the $s=-2, \ell=2, m=2, a/M=0.68$ fundamental tone can be obtained using
```julia
Rin, Rup = Teukolsky_radial(-2, 2, 2, 0.68, 0.5239751-0.0815126im)
```
We can check out the $R^{\textrm{up}}$ solution using
```julia
julia> Rup
TeukolskyRadialFunction(
    mode=Mode(s=-2, l=2, m=2, a=0.68, omega=0.5239751 - 0.0815126im, lambda=1.655003080578682 + 0.3602676563885877im),
    boundary_condition=UP,
    transmission_amplitude=1.0 + 0.0im,
    incidence_amplitude=-5.850900444651249e-8 - 3.80716581300155e-7im,
    reflection_amplitude=1.1011632133920028 + 2.1300597377432497im,
    normalization_convention=UNIT_TEUKOLSKY_TRANS
)
```
We see that the incidence amplitude is indeed very small numerically as a QNM solution should. This can be accessed using
```julia
Rup.incidence_amplitude
```

This should give
```julia
-5.850900444651249e-8 - 3.80716581300155e-7im
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
    month = "3",
    year = "2025"
}
```

## License
The package is licensed under the MIT License.
