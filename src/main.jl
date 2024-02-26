"""
    simulate(raft_parameters; rhs, alg, abstol, reltol, return_raw)

Simulate a Sargassum raft with [`RaftParameters`](@ref) `raft_parameters` for a time interval `tspan` and return a [`RaftTrajectory`](@ref).

This function modifies `raft_parameters` significantly.

### Arguments 

- `raft_parameters`: A [`RaftParameters`](@ref) defining the raft.

### Optional Arguments

- `rhs`: The function to integrate, default [`Raft!`](@ref). Use [`Leeway!`](@ref) to integrate particles with no springs or inertia.
- `alg`: The integration algorithm to use, default [`Tsit5()`](https://docs.sciml.ai/DiffEqDocs/stable/solvers/ode_solve/).
- `abstol`: The absolute tolerance of integration; default `1e-6`. Can be `nothing` to use `OrdinaryDiffEq.solve` defaults.
- `reltol`: The relative tolerance of integration; default `1e-6`. Can be `nothing` to use `OrdinaryDiffEq.solve` defaults.
- `showprogress`: If `true`, print a status indicator of the progress of the integration. Default `true`.
- `return_raw`: If true, return the result of `OrdinaryDiffEq.solve`, rather than a [`RaftTrajectory`](@ref). \
Use this if you would like to manipulate the solution directly. Default `false`.
"""
function simulate(
    rp::RaftParameters; 
    rhs::Function = Raft!,
    alg::OrdinaryDiffEq.OrdinaryDiffEqAlgorithm = Tsit5(),
    abstol::Union{Real, Nothing} = 1e-6,
    reltol::Union{Real, Nothing} = 1e-6,
    showprogress::Bool = true,
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
        return RaftTrajectory(sol, rp, EQR_DEFAULT, dt = 0.1)
    end
end

