"""
    simulate(raft_parameters; rhs, alg, abstol, reltol, return_raw)

Simulate a Sargassum raft with [`RaftParameters`](@ref) `raft_parameters` for a time interval `tspan` and return a [`RaftTrajectory`](@ref).

This function modifies `raft_parameters` significantly.

### Arguments 

- `raft_parameters`: A [`RaftParameters`](@ref) defining the raft.

### Optional Arguments

- `rhs`: The function to integrate, default [`Raft!`](@ref). Use [`WaterWind!`](@ref) to integrate particles with no springs or inertia.
- `alg`: The integration algorithm to use, default [`Tsit5()`](https://docs.sciml.ai/DiffEqDocs/stable/solvers/ode_solve/).
- `abstol`: The absolute tolerance of integration; default `1e-6`.
- `reltol`: The relative tolerance of integration; default `1e-6`.
- `showprogress`: If `true`, print a status indicator of the progress of the integration. Default `true`.
- `return_raw`: If true, return the result of `OrdinaryDiffEq.solve`, rather than a [`RaftTrajectory`](@ref). \
Use this if you would like to manipulate the solution directly. Default `false`.
"""
function simulate(
    rp::RaftParameters; 
    rhs::Function = Raft!,
    alg::OrdinaryDiffEq.OrdinaryDiffEqAlgorithm = Tsit5(),
    abstol::Real = 1e-6,
    reltol::Real = 1e-6,
    showprogress::Bool = true,
    return_raw::Bool = false)

    prob_raft = ODEProblem(rhs, rp.ics.ics, rp.tspan, rp)

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

"""
    yearmonth2tspan(ym_initial, ym_final; t_extra)

Convert two "yearmonths" (tuples of the form `(year, month)`) into a single tuple `(t_initial, t_final)` suitable for 
use in [`simulate`](@ref).

The times are meaured with respect to `WATER_ITP.x.time_start`.

### Optional Arguments 

- `t_extra`: An extra number of days added to the initial and final times. Default `(0, 0)`.
"""
function yearmonth2tspan(
    ym_initial::Tuple{Integer, Integer}, 
    ym_final::Tuple{Integer, Integer}; 
    t_extra::Tuple{Integer, Integer} = (0, 0))
    t1 = Day(DateTime(ym_initial...) - WATER_ITP.x.time_start).value |> float
    t2 = t1 + Day(DateTime(ym_final...) - DateTime(ym_initial...)).value |> float

    return (t1 + t_extra[1], t2 + t_extra[2])
end