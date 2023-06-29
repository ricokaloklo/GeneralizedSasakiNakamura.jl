var documenterSearchIndex = {"docs":
[{"location":"APIs/#APIs","page":"APIs","title":"APIs","text":"","category":"section"},{"location":"APIs/","page":"APIs","title":"APIs","text":"There are 4 functions that are exported, namely","category":"page"},{"location":"APIs/","page":"APIs","title":"APIs","text":"Teukolsky_radial\nGSN_radial\nrstar_from_r\nr_from_rstar","category":"page"},{"location":"APIs/","page":"APIs","title":"APIs","text":"and there are 5 custom types that are exported, i.e.","category":"page"},{"location":"APIs/","page":"APIs","title":"APIs","text":"BoundaryCondition\nNormalizationConvention\nMode\nGSNRadialFunction\nTeukolskyRadialFunction","category":"page"},{"location":"APIs/","page":"APIs","title":"APIs","text":"Currently, only the exported functions and types are documented below. Documentations for private (i.e. unexported) functions will be added at a later stage.","category":"page"},{"location":"APIs/#Functions","page":"APIs","title":"Functions","text":"","category":"section"},{"location":"APIs/","page":"APIs","title":"APIs","text":"Teukolsky_radial","category":"page"},{"location":"APIs/#GeneralizedSasakiNakamura.Teukolsky_radial","page":"APIs","title":"GeneralizedSasakiNakamura.Teukolsky_radial","text":"Teukolsky_radial(s::Int, l::Int, m::Int, a, omega, boundary_condition, rsin, rsout; horizon_expansion_order::Int=3, infinity_expansion_order::Int=6, data_type=Solutions._DEFAULTDATATYPE,  ODE_algorithm=Solutions._DEFAULTSOLVER, tolerance=Solutions._DEFAULTTOLERANCE)\n\nCompute the Teukolsky function for a given mode (specified by s the spin weight, l the harmonic index, m the azimuthal index, a the Kerr spin parameter, and omega the frequency)  and boundary condition (specified by boundary_condition which can be either IN for purely-ingoing at the horizon or UP for purely-outgoing at infinity) using the GSN formalism with GSN_radial.\n\nThe full GSN solution is converted to the corresponding Teukolsky solution (R(r) dRdr) and  the incidence, reflection and transmission amplitude are converted from the GSN formalism to the Teukolsky formalism  with the normalization convention that the transmission amplitude is normalized to 1 (i.e. normalization_convention=UNIT_TEUKOLSKY_TRANS). \n\n\n\n\n\n","category":"function"},{"location":"APIs/","page":"APIs","title":"APIs","text":"GSN_radial","category":"page"},{"location":"APIs/#GeneralizedSasakiNakamura.GSN_radial","page":"APIs","title":"GeneralizedSasakiNakamura.GSN_radial","text":"GSN_radial(s::Int, l::Int, m::Int, a, omega, boundary_condition, rsin, rsout; horizon_expansion_order::Int=3, infinity_expansion_order::Int=6, data_type=Solutions._DEFAULTDATATYPE,  ODE_algorithm=Solutions._DEFAULTSOLVER, tolerance=Solutions._DEFAULTTOLERANCE)\n\nCompute the GSN function for a given mode (specified by s the spin weight, l the harmonic index, m the azimuthal index, a the Kerr spin parameter, and omega the frequency)  and boundary condition (specified by boundary_condition which can be either IN for purely-ingoing at the horizon or UP for purely-outgoing at infinity).\n\nThe GSN function is numerically solved in the interval of tortoise coordinates r_* in [rsin, rsout] using the ODE solver (from DifferentialEquations.jl) specified by ODE_algorithm (default: Vern9())  with tolerance specified by tolerance (default: 1e-12). By default the data type used is ComplexF64 (i.e. double-precision floating-point number) but it can be changed by  specifying data_type (e.g. Complex{BigFloat} for complex arbitrary precision number).\n\nWhile the numerical GSN solution is only accurate in the range [rsin, rsout],  the full GSN solution is constructed by smoothly attaching the asymptotic solutions near horizon (up to horizon_expansion_order-th order)  and infinity (up to infinity_expansion_order-th order). Therefore, the now-semi-analytical GSN solution is accurate everywhere.\n\nReturn a GSNRadialFunction object which contains all the information about the GSN solution.\n\n\n\n\n\n","category":"function"},{"location":"APIs/","page":"APIs","title":"APIs","text":"rstar_from_r","category":"page"},{"location":"APIs/#GeneralizedSasakiNakamura.Coordinates.rstar_from_r","page":"APIs","title":"GeneralizedSasakiNakamura.Coordinates.rstar_from_r","text":"rstar_from_r(a, r)\n\nConvert a Boyer-Lindquist coordinate r to the corresponding tortoise coordinate rstar.\n\n\n\n\n\n","category":"function"},{"location":"APIs/","page":"APIs","title":"APIs","text":"r_from_rstar","category":"page"},{"location":"APIs/#GeneralizedSasakiNakamura.Coordinates.r_from_rstar","page":"APIs","title":"GeneralizedSasakiNakamura.Coordinates.r_from_rstar","text":"r_from_rstar(a, rstar)\n\nConvert a tortoise coordinate rstar to the corresponding Boyer-Lindquist coordiante r.  It uses a bisection method when rstar <= 0, and Newton method otherwise.\n\nThe function assumes that r geq r_+ where r_+ is the outer event horizon.\n\n\n\n\n\n","category":"function"},{"location":"APIs/#Types","page":"APIs","title":"Types","text":"","category":"section"},{"location":"APIs/#BoundaryCondition","page":"APIs","title":"BoundaryCondition","text":"","category":"section"},{"location":"APIs/","page":"APIs","title":"APIs","text":"This is an enum type that can take either one of the two values","category":"page"},{"location":"APIs/","page":"APIs","title":"APIs","text":"value \nIN purely ingoing at the horizon\nUP purely outgoing at infinity","category":"page"},{"location":"APIs/#NormalizationConvention","page":"APIs","title":"NormalizationConvention","text":"","category":"section"},{"location":"APIs/","page":"APIs","title":"APIs","text":"This is an enum type that can take either one of the two values","category":"page"},{"location":"APIs/","page":"APIs","title":"APIs","text":"value \nUNIT_GSN_TRANS normalized to have a unit transmission amplitude for the GSN function\nUNIT_TEUKOLSKY_TRANS normalized to have a unit transmission amplitude for the Teukolsky function","category":"page"},{"location":"APIs/#Mode","page":"APIs","title":"Mode","text":"","category":"section"},{"location":"APIs/","page":"APIs","title":"APIs","text":"This is a composite struct type that stores information about a mode","category":"page"},{"location":"APIs/","page":"APIs","title":"APIs","text":"field \ns spin weight s\nl harmonic index ell\nm azimuthal index m\na Kerr spin parameter aM\nomega frequency Momega\nlambda spin-weighted spheroidal eigenvalue lambda","category":"page"},{"location":"APIs/#GSNRadialFunction","page":"APIs","title":"GSNRadialFunction","text":"","category":"section"},{"location":"APIs/","page":"APIs","title":"APIs","text":"This is a composite struct type that stores the output from GSN_radial","category":"page"},{"location":"APIs/","page":"APIs","title":"APIs","text":"note: Note\nGSNRadialFunction(rstar) is equivalent to GSNRadialFunction.GSN_solution(rstar)[1],  returning only the value of the GSN function evaluated at the tortoise coordinate rstar","category":"page"},{"location":"APIs/","page":"APIs","title":"APIs","text":"field \nmode a Mode object storing information about the mode\nboundary_condition a BoundaryCondition object storing which boundary condition this function satisfies\nrsin numerical inner boundary r_*^mathrminM where the GSN equation is numerically evolved (r_* is a tortoise coordinate)\nrsout numerical outer boundary r_*^mathrmoutM where the GSN equation is numerically evolved (r_* is a tortoise coordinate)\nhorizon_expansion_order order of the asymptotic expansion at the horizon\ninfinity_expansion_order order of the asymptotic expansion at infinity\ntransmission_amplitude transmission amplitude in the GSN formalism of this function\nincidence_amplitude incidence amplitude in the GSN formalism of this function\nreflection_amplitude reflection amplitude in the GSN formalism of this function\nnumerical_GSN_solution numerical solution (ODESolution object from DifferentialEquations.jl) to the GSN equation in [rsin, rsout] if applicable; output is a vector  hatX(r_*) dhatX(r_*)dr_* \nnumerical_Riccati_solution numerical solution (ODESolution object from DifferentialEquations.jl) to the GSN equation in the Riccati form if applicable; output is a vector  hatPhi(r_*) dhatPhi(r_*)dr_* \nGSN_solution full GSN solution where asymptotic solutions are smoothly attached; output is a vector  hatX(r_*) dhatX(r_*)dr_* \nnormalization_convention a NormalizationConvention object storing which normalization convention this function adheres to","category":"page"},{"location":"APIs/#TeukolskyRadialFunction","page":"APIs","title":"TeukolskyRadialFunction","text":"","category":"section"},{"location":"APIs/","page":"APIs","title":"APIs","text":"This is a composite struct type that stores the output from Teukolsky_radial","category":"page"},{"location":"APIs/","page":"APIs","title":"APIs","text":"note: Note\nTeukolskyRadialFunction(r) is equivalent to TeukolskyRadialFunction.Teukolsky_solution(r)[1],  returning only the value of the Teukolsky function evaluated at the Boyer-Lindquist coordinate r","category":"page"},{"location":"APIs/","page":"APIs","title":"APIs","text":"field \nmode a Mode object storing information about the mode\nboundary_condition a BoundaryCondition object storing which boundary condition this function satisfies\ntransmission_amplitude transmission amplitude in the Teukolsky formalism of this function\nincidence_amplitude incidence amplitude in the Teukolsky formalism of this function\nreflection_amplitude reflection amplitude in the Teukolsky formalism of this function\nGSN_solution a GSNRadialFunction object storing the corresponding GSN function\nTeukolsky_solution Teukolsky solution where asymptotic solutions are smoothly attached; output is a vector  hatR(r) dhatR(r)dr \nnormalization_convention a NormalizationConvention object storing which normalization convention this function adheres to","category":"page"},{"location":"examples/#Examples","page":"Examples","title":"Examples","text":"","category":"section"},{"location":"#Home","page":"Home","title":"Home","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"GeneralizedSasakiNakamura.jl computes solutions to the frequency-domain radial Teukolsky equation with the Generalized Sasaki-Nakamura (GSN) formalism.","category":"page"},{"location":"","page":"Home","title":"Home","text":"The code is capable of handling both in-going and out-going radiation of scalar, electromagnetic, and gravitational type (corresponding to spin weight of s = 0 pm 1 pm 2 respectively).","category":"page"},{"location":"","page":"Home","title":"Home","text":"The angular Teukolsky equation is solved with an accompanying julia package SpinWeightedSpheroidalHarmonics.jl using a spectral decomposition method.","category":"page"},{"location":"#Installation","page":"Home","title":"Installation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"To install the package using the Julia package manager, simply type the following in the Julia REPL:","category":"page"},{"location":"","page":"Home","title":"Home","text":"using Pkg\nPkg.add(\"GeneralizedSasakiNakamura\")","category":"page"},{"location":"","page":"Home","title":"Home","text":"Note: There is no need to install SpinWeightedSpheroidalHarmonics.jl separately as it should be automatically installed by the package manager.","category":"page"},{"location":"#Highlights","page":"Home","title":"Highlights","text":"","category":"section"},{"location":"#Performant-frequency-domain-Teukolsky-solver","page":"Home","title":"Performant frequency-domain Teukolsky solver","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Takes on average only a few tens of milliseconds:","category":"page"},{"location":"","page":"Home","title":"Home","text":"<table>\n  <tr>\n    <th>GeneralizedSasakiNakamura.jl</th>\n    <th><a href=\"https://github.com/BlackHolePerturbationToolkit/Teukolsky\">Teukolsky</a> Mathematica package using the MST method </th>\n  </tr>\n  <tr>\n    <td><p align=\"center\"><img width=\"100%\" src=\"https://github-production-user-asset-6210df.s3.amazonaws.com/55488840/248965077-7d216deb-5bae-433f-a699-d40a35f0e35d.gif\"></p></td>\n    <td><p align=\"center\"><img width=\"100%\" src=\"https://github-production-user-asset-6210df.s3.amazonaws.com/55488840/248966033-9e7d8027-81ee-4762-98d9-0ad0a1c030ad.gif\"></p></td>\n  </tr>\n</table>","category":"page"},{"location":"","page":"Home","title":"Home","text":"(There was no caching! We solved the equation on-the-fly! The notebook generating this animation can be found here)","category":"page"},{"location":"#Solutions-that-are-accurate-everywhere","page":"Home","title":"Solutions that are accurate everywhere","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Numerical solutions are smoothly stitched to analytical ansatzes near the horizon and infinity at user-specified locations rsin and rsout respectively:","category":"page"},{"location":"","page":"Home","title":"Home","text":"<p align=\"center\">\n  <img width=\"50%\" src=\"https://github-production-user-asset-6210df.s3.amazonaws.com/55488840/248724944-9707332b-1238-4b3b-b1c0-ac426a1b3dc6.gif\">\n</p>","category":"page"},{"location":"#Easy-to-use","page":"Home","title":"Easy to use","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"The following code snippet lets you solve the (source-free) Teukolsky function (in frequency domain) for the mode s=-2 ell=2 m=2 a=07 omega=05 that satisfies the purely-ingoing boundary condition at the horizon:","category":"page"},{"location":"","page":"Home","title":"Home","text":"using GeneralizedSasakiNakamura # This is going to take some time to pre-compile, mostly due to DifferentialEquations.jl\n\n# Specify which mode and what boundary condition\ns=-2; l=2; m=2; a=0.7; omega=0.5; bc=IN;\n# Specify where to match to ansatzes\nrsin=-20; rsout=250;\n\n# NOTE: julia uses 'just-ahead-of-time' compilation. Calling this the first time in each session will take some time\nR = Teukolsky_radial(s, l, m, a, omega, bc, rsin, rsout) ","category":"page"},{"location":"","page":"Home","title":"Home","text":"That's it! If you run this on Julia REPL, it should give you something like this","category":"page"},{"location":"","page":"Home","title":"Home","text":"TeukolskyRadialFunction(\n    mode=Mode(s=-2, l=2, m=2, a=0.7, omega=0.5, lambda=1.6966094016353415),\n    boundary_condition=IN,\n    transmission_amplitude=1.0 + 0.0im,\n    incidence_amplitude=6.536587661197995 - 4.941203897068852im,\n    reflection_amplitude=-0.128246619129379 - 0.44048133496664404im,\n    normalization_convention=UNIT_TEUKOLSKY_TRANS\n)","category":"page"},{"location":"","page":"Home","title":"Home","text":"For example, if we want to evaluate the Teukolsky function at the location r = 10M, simply do","category":"page"},{"location":"","page":"Home","title":"Home","text":"R(10)","category":"page"},{"location":"","page":"Home","title":"Home","text":"This should give","category":"page"},{"location":"","page":"Home","title":"Home","text":"77.57508416832009 - 429.40290952257226im","category":"page"},{"location":"#License","page":"Home","title":"License","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"The package is licensed under the MIT License.","category":"page"}]
}
