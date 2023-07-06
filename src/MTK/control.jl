"""
    avoid_shore(water_interpolant; tol)

Construct a `DiscreteCallback` which terminates integration when both the `u` and `v` components of `water_interpolant` are less than `tol`.
"""
function avoid_shore(
    water_interpolant::VectorField2DInterpolantEQR; 
    tol::Real = 0.1, 
    save_positions::NTuple{2, Bool} = (true, false))
    condition(u, t, integrator) = (abs(water_interpolant.u(u..., t)) < tol) && (abs(water_interpolant.v(u..., t)) < tol)
    affect!(integrator) = terminate!(integrator)
    return DiscreteCallback(condition, affect!, save_positions = save_positions)
end