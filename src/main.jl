"""
    simulate(rp; rhs, alg, abstol, reltol, showprogress, dt, return_raw)

Simulate a Sargassum raft with [`RaftParameters`](@ref) `rp` and return a [`RaftTrajectory`](@ref).

This function modifies `rp` significantly.

### Arguments 

- `rp`: A [`RaftParameters`](@ref) defining the raft.

### Optional Arguments

- `rhs`: The function to integrate, default [`Raft!`](@ref). Use [`Leeway!`](@ref) to integrate particles with no springs or inertia.
- `alg`: The integration algorithm to use, default [`Tsit5()`](https://docs.sciml.ai/DiffEqDocs/stable/solvers/ode_solve/).
- `abstol`: The absolute tolerance of integration; default `nothing`.
- `reltol`: The relative tolerance of integration; default `nothing`.
- `showprogress`: If `true`, print a status indicator of the progress of the integration. Default `true`.
- `dt`: The solution trajectories are uniformized (after integration) to be spaced in time by increments of `dt`. Note that the units of \
this quantity are implicity `UNITS["time"]`. Default `0.1`.
- `return_raw`: If true, return the result of `OrdinaryDiffEq.solve`, rather than a [`RaftTrajectory`](@ref). \
Use this if you would like to manipulate the solution directly. Default `false`.
"""
function simulate(
    rp::RaftParameters; 
    rhs::Function = Raft!,
    alg::OrdinaryDiffEq.OrdinaryDiffEqAlgorithm = Tsit5(),
    abstol::Union{Real, Nothing} = nothing,
    reltol::Union{Real, Nothing} = nothing,
    showprogress::Bool = true,
    dt::Real = 0.1,
    return_raw::Bool = false)

    prob_raft = ODEProblem(rhs, rp.ics.ics, rp.ics.tspan, rp)

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

        return RaftTrajectory(sol, rp, dt = dt)
    end
end

