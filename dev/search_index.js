var documenterSearchIndex = {"docs":
[{"location":"APIs/#APIs","page":"APIs","title":"APIs","text":"","category":"section"},{"location":"APIs/#Main-functions","page":"APIs","title":"Main functions","text":"","category":"section"},{"location":"APIs/","page":"APIs","title":"APIs","text":"Teukolsky_radial","category":"page"},{"location":"APIs/#GeneralizedSasakiNakamura.Teukolsky_radial","page":"APIs","title":"GeneralizedSasakiNakamura.Teukolsky_radial","text":"Teukolsky_radial(s::Int, l::Int, m::Int, a, omega, boundary_condition, rsin, rsout; horizon_expansion_order::Int=3, infinity_expansion_order::Int=6, data_type=Solutions._DEFAULTDATATYPE,  ODE_algorithm=Solutions._DEFAULTSOLVER, tolerance=Solutions._DEFAULTTOLERANCE)\n\nCompute the Teukolsky function for a given mode (specified by s the spin weight, l the harmonic index, m the azimuthal index, a the Kerr spin parameter, and omega the frequency)  and boundary condition (specified by boundary_condition which can be either IN for purely-ingoing at the horizon or UP for purely-outgoing at infinity) using the GSN formalism with GSN_radial.\n\nThe full GSN solution is converted to the corresponding Teukolsky solution (R(r) dRdr) and  the incidence, reflection and transmission amplitude are converted from the GSN formalism to the Teukolsky formalism  with the normalization convention that the transmission amplitude is normalized to 1 (i.e. normalization_convention=UNIT_TEUKOLSKY_TRANS). \n\n\n\n\n\n","category":"function"},{"location":"APIs/","page":"APIs","title":"APIs","text":"GSN_radial","category":"page"},{"location":"APIs/#GeneralizedSasakiNakamura.GSN_radial","page":"APIs","title":"GeneralizedSasakiNakamura.GSN_radial","text":"GSN_radial(s::Int, l::Int, m::Int, a, omega, boundary_condition, rsin, rsout; horizon_expansion_order::Int=3, infinity_expansion_order::Int=6, data_type=Solutions._DEFAULTDATATYPE,  ODE_algorithm=Solutions._DEFAULTSOLVER, tolerance=Solutions._DEFAULTTOLERANCE)\n\nCompute the GSN function for a given mode (specified by s the spin weight, l the harmonic index, m the azimuthal index, a the Kerr spin parameter, and omega the frequency)  and boundary condition (specified by boundary_condition which can be either IN for purely-ingoing at the horizon or UP for purely-outgoing at infinity).\n\nThe GSN function is numerically solved in the interval of tortoise coordinates r_* in [rsin, rsout] using the ODE solver (from DifferentialEquations.jl) specified by ODE_algorithm (default: Vern9())  with tolerance specified by tolerance (default: 1e-12). By default the data type used is ComplexF64 (i.e. double-precision floating-point number) but it can be changed by  specifying data_type (e.g. Complex{BigFloat} for complex arbitrary precision number).\n\nWhile the numerical GSN solution is only accurate in the range [rsin, rsout],  the full GSN solution is constructed by smoothly attaching the asymptotic solutions near horizon (up to horizon_expansion_order-th order)  and infinity (up to infinity_expansion_order-th order). Therefore, the now-semi-analytical GSN solution is accurate everywhere.\n\nReturn a GSNRadialFunction object which contains all the information about the GSN solution.\n\n\n\n\n\n","category":"function"},{"location":"APIs/","page":"APIs","title":"APIs","text":"rstar_from_r","category":"page"},{"location":"APIs/#GeneralizedSasakiNakamura.Coordinates.rstar_from_r","page":"APIs","title":"GeneralizedSasakiNakamura.Coordinates.rstar_from_r","text":"rstar_from_r(a, r)\n\nConvert a Boyer-Lindquist coordinate r to the corresponding tortoise coordinate rstar.\n\n\n\n\n\n","category":"function"},{"location":"APIs/","page":"APIs","title":"APIs","text":"r_from_rstar","category":"page"},{"location":"APIs/#GeneralizedSasakiNakamura.Coordinates.r_from_rstar","page":"APIs","title":"GeneralizedSasakiNakamura.Coordinates.r_from_rstar","text":"r_from_rstar(a, rstar)\n\nConvert a tortoise coordinate rstar to the corresponding Boyer-Lindquist coordiante r.  It uses a bisection method when rstar <= 0, and Newton method otherwise.\n\nThe function assumes that r geq r_+ where r_+ is the outer event horizon.\n\n\n\n\n\n","category":"function"},{"location":"examples/#Examples","page":"Examples","title":"Examples","text":"","category":"section"},{"location":"#Home","page":"Home","title":"Home","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"GeneralizedSasakiNakamura.jl computes solutions to the frequency-domain radial Teukolsky equation with the Generalized Sasaki-Nakamura (GSN) formalism.","category":"page"},{"location":"","page":"Home","title":"Home","text":"The angular Teukolsky equation is solved with another julia package SpinWeightedSpheroidalHarmonics.jl using a spectral decomposition method.","category":"page"},{"location":"#Installation","page":"Home","title":"Installation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"To install the package using the Julia package manager, simply type the following in the Julia REPL:","category":"page"},{"location":"","page":"Home","title":"Home","text":"using Pkg\nPkg.add(\"GeneralizedSasakiNakamura\")","category":"page"},{"location":"","page":"Home","title":"Home","text":"Note: There is no need to install SpinWeightedSpheroidalHarmonics.jl manually as it should be automatically installed by the package manager.","category":"page"},{"location":"#Highlights","page":"Home","title":"Highlights","text":"","category":"section"},{"location":"#Performant-frequency-domain-Teukolsky-solver","page":"Home","title":"Performant frequency-domain Teukolsky solver","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Takes on average only a few tens of milliseconds: <p align=\"center\">   <img width=\"50%\" src=\"https://github-production-user-asset-6210df.s3.amazonaws.com/55488840/248602615-eb611e3f-00b9-49c3-905c-0440b0607683.gif\"> </p>","category":"page"},{"location":"","page":"Home","title":"Home","text":"(There was no caching! We solved the equation on-the-fly! The notebook generating this animation can be found here)","category":"page"},{"location":"#Solutions-that-are-accurate-everywhere","page":"Home","title":"Solutions that are accurate everywhere","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Numerical solutions are smoothly stitched to analytical ansatzes near the horizon and infinity at user-specified locations rsin and rsout respectively:","category":"page"},{"location":"","page":"Home","title":"Home","text":"<p align=\"center\">   <img width=\"50%\" src=\"https://github-production-user-asset-6210df.s3.amazonaws.com/55488840/248601145-8ce9b236-d197-4321-83b9-fcb9692a51e8.gif\"> </p>","category":"page"},{"location":"#Easy-to-use","page":"Home","title":"Easy to use","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"The following code snippet lets you solve the (source-free) Teukolsky function (in frequency domain) for the mode s=-2 ell=2 m=2 a=07 omega=05 that satisfies the purely-ingoing boundary condition at the horizon:","category":"page"},{"location":"","page":"Home","title":"Home","text":"using GeneralizedSasakiNakamura # This is going to take some time to pre-compile, mostly due to DifferentialEquations.jl\n\n# Specify which mode and what boundary condition\ns=-2; l=2; m=2; a=0.7; omega=0.5; bc=IN;\n# Specify where to match to ansatzes\nrsin=-20; rsout=250;\n\n# NOTE: julia uses 'just-ahead-of-time' compilation. Calling this the first time in each session will take some time\nR = Teukolsky_radial(s, l, m, a, omega, bc, rsin, rsout) ","category":"page"},{"location":"","page":"Home","title":"Home","text":"That's it!","category":"page"},{"location":"#License","page":"Home","title":"License","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"The package is licensed under the MIT License.","category":"page"}]
}
