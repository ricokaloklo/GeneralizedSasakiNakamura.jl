module QNM

using LinearAlgebra
using Printf
using SpinWeightedSpheroidalHarmonics
using ..GeneralizedSasakiNakamura:
    GSN_radial, Teukolsky_radial, Y_radial, IN, UP,
    ConversionFactors, Transformation, r_from_rstar, rstar_from_r, ISEM,
    UNIT_GSN_TRANS, UNIT_TEUKOLSKY_TRANS, ExtremalRadial, _swsh_eigenvalue

# Portions of ContinuedFractions.jl are translated from qnm/radial.py at
# commit f3abd18e59828e7e7d75d07f20c7cbc87925edfa.
#
# MIT License
# Copyright (c) 2019 Leo C. Stein
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

include("Types.jl")
include("ContinuedFractions.jl")
include("AngularAdapter.jl")
include("RadialLeaver.jl")
include("ExtremalLeaver.jl")
include("ExtremalRootTable.jl")
include("SpinRootTable.jl")
include("NearExtremalExpansion.jl")
include("RootPolishing.jl")
include("NegativeImaginaryAxis.jl")
include("ModeSequences.jl")
include("ISEMValidation.jl")
include("ExcitationFactors.jl")
include("VariationalExcitation.jl")
include("VariationalObservables.jl")
include("ModeObservables.jl")

export QNMMode, QNMBranch, QNMResult, QNMEstimate, QNMFailure,
    QNMEndpointResult,
    QNMPairResult
export ordinary, mirror, Ordinary, Mirror
export LeaverResult, ISEMValidationResult, ExcitationFactorResult
export qnm_frequency, qnm_sequence, validate_qnm_with_isem
export qnm_excitation_factor
export qnm, qnm_pair
export angular_A_to_lambda, lambda_to_angular_A

end
