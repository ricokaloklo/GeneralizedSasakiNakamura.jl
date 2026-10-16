# Leaver radial residual at (a, omega): solve the angular problem for A on
# `sheet_id`, then evaluate the continued fraction with inversion index
# `inversion_index`. Its zero is a QNM; returns the value, CF error and angular data.
function leaver_radial_residual(mode::QNMMode, a, omega;
        angular_order::Int=32,
        sheet_id::Symbol=:qnm_straight_from_spherical,
        inversion_index::Int=mode.n,
        cf_tolerance=nothing,
        cf_minimum_iterations::Int=300,
        cf_maximum_iterations::Int=4000,
        require_cf_convergence::Bool=true)
    angular = angular_branch(mode, a, omega;
        truncation_order=angular_order, sheet_id)
    angular.status in (
        :continued, :predictor_corrected, :spherical_anchor,
        :high_precision_refined) ||
        throw(ErrorException(
            "Angular continuation failed with status $(angular.status)."))
    fraction = leaver_cf_inversion(
        omega, a, mode.s, mode.m, angular.angular_A, inversion_index;
        tolerance=cf_tolerance,
        minimum_iterations=cf_minimum_iterations,
        maximum_iterations=cf_maximum_iterations,
        require_convergence=require_cf_convergence,
    )
    return (
        value=fraction.value,
        error=fraction.error,
        iterations=fraction.iterations,
        cf_converged=fraction.converged,
        inversion_index=inversion_index,
        angular=angular,
    )
end
