"""
    simulate(rp, tspan; rhs, alg, abstol, reltol, return_raw)

Simulate a Sargassum raft with [`RaftParameters`](@ref) `rp` for a time interval `tspan` and return a [`RaftTrajectory`](@ref).

### Arguments 

- `rp`: A [`RaftParameters`](@ref) defining the raft.
- `tspan`: A `Tuple{Real, Real}` such that the integration is performed for `tspan[1] ≤ t ≤ tspan[2]` where `t` is in days.

### Optional Arguments

- `rhs`: The function to integrate, default [`Raft!`](@ref). Use [`WaterWind!`](@ref) to integrate particles with no springs or inertia.
- `alg`: The integration algorithm to use, default [`Tsit5()`](https://docs.sciml.ai/DiffEqDocs/stable/solvers/ode_solve/).
- `abstol`: The absolute tolerance of integration; default `1e-6`.
- `reltol`: The relative tolerance of integration; default `1e-6`.
- `showprogress`: If `true`, print a status indicator of the progress of the integration. Default `true`.
`return_raw`: If true, return the result of `OrdinaryDiffEq.solve`, rather than a [`RaftTrajectory`](@ref). Use this if you would
    like to manipulate the solution directly. Default `false`.
"""
function simulate(
    rp::RaftParameters, 
    tspan::Tuple{Real, Real}; 
    rhs::Function = Raft!,
    alg::OrdinaryDiffEq.OrdinaryDiffEqAlgorithm = Tsit5(),
    abstol::Real = 1e-6,
    reltol::Real = 1e-6,
    showprogress = true,
    return_raw::Bool = false)

    @assert tspan[1] < tspan[2] "initial time must be less than final time"

    prob_raft = ODEProblem(rhs, rp.ics.ics, tspan, rp)

    sol = solve(
        prob_raft, 
        alg, abstol = abstol, reltol = reltol,
        callback = CallbackSet(
            cb_update(showprogress = showprogress), 
            cb_land(rp.land), 
            cb_growth_death(rp.gd_model), 
            cb_connections())
        )

    if return_raw
        return sol
    else
        return RaftTrajectory(sol, rp, EQR_DEFAULT, dt = 0.1)
    end
end