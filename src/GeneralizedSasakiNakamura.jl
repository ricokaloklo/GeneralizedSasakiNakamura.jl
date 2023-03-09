module GeneralizedSasakiNakamura

include("Kerr.jl")
include("Coordinates.jl")
include("AsymptoticExpansionCoefficients.jl")
include("InitialConditions.jl")
include("ConversionFactors.jl")
include("Potentials.jl")
include("Transformation.jl")
include("Solutions.jl")

using .Coordinates
export r_from_rstar, rstar_from_r

using .Potentials

using .Solutions
export solve_Xup, solve_Xin, Teukolsky_radial_function_from_Sasaki_Nakamura_function, Teukolsky_radial_function_from_Sasaki_Nakamura_function_conversion_matrix
solve_Xup = Solutions.solve_Xup
solve_Xin = Solutions.solve_Xin
Teukolsky_radial_function_from_Sasaki_Nakamura_function = Solutions.Teukolsky_radial_function_from_Sasaki_Nakamura_function
Teukolsky_radial_function_from_Sasaki_Nakamura_function_conversion_matrix = Solutions.Teukolsky_radial_function_from_Sasaki_Nakamura_function_conversion_matrix
export d2Rdr2_from_Rsoln
d2Rdr2_from_Rsoln = Solutions.d2Rdr2_from_Rsoln
export scaled_Wronskian
scaled_Wronskian = Solutions.scaled_Wronskian
export CincCref_SN_from_Xup, BrefBinc_SN_from_Xin
CincCref_SN_from_Xup = Solutions.CincCref_SN_from_Xup
BrefBinc_SN_from_Xin = Solutions.BrefBinc_SN_from_Xin

end
