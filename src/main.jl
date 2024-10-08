"""
    simulate(rp; rhs, alg, abstol, reltol, showprogress, dt, return_raw)

Simulate a Sargassum raft with [`RaftParameters`](@ref) `rp` and return a [`RaftTrajectory`](@ref).

This function modifies the fields of `rp` significantly.

### Arguments 

- `rp`: A [`RaftParameters`](@ref) defining the raft.

### Optional Arguments

- `leeway`: Use [`Leeway!`](@ref) to integrate particles with no springs or inertia. Default `false`.
- `alg`: The integration algorithm to use, default [`Tsit5()`](https://docs.sciml.ai/DiffEqDocs/stable/solvers/ode_solve/).
- `abstol`: The absolute tolerance of integration; default `nothing`.
- `reltol`: The relative tolerance of integration; default `nothing`.
- `showprogress`: If `true`, print a status indicator of the progress of the integration. Default `false`.
- `dt`: The solution trajectories are uniformized to be spaced in time by increments of `dt`. Note that the units of \
this quantity are implicity `UNITS["time"]`. Default `0.1`.
- `return_raw`: If true, return the result of `OrdinaryDiffEq.solve`, rather than a [`RaftTrajectory`](@ref). \
Use this if you would like to manipulate the solution directly. Default `false`.
"""
function simulate(
    rp::RaftParameters; 
    leeway::Bool = false,
    alg::OrdinaryDiffEq.OrdinaryDiffEqAlgorithm = Tsit5(),
    abstol::Union{Real, Nothing} = nothing,
    reltol::Union{Real, Nothing} = nothing,
    showprogress::Bool = false,
    dt::Real = 0.1,
    return_raw::Bool = false)

    tspan = rp.ics.tspan

    if leeway
        prob_raft = ODEProblem(Leeway!, rp.ics.ics, tspan, rp)
    elseif rp.dx_MR !== nothing && rp.dx_MR !== nothing
        prob_raft = ODEProblem(FastRaft!, rp.ics.ics, tspan, rp)
    else
        prob_raft = ODEProblem(Raft!, rp.ics.ics, tspan, rp)
    end

    callback = CallbackSet(
                DiscreteCallback(rp.land, rp.land),                             # LAND
                DiscreteCallback(rp.gd_model, rp.gd_model),                     # BIOLOGY 
                DiscreteCallback((u, t, integrator) -> true, rp.connections))   # CONNECTIONS
    integrator = init(prob_raft, alg, abstol = abstol, reltol = reltol, callback = callback)

    ts = range(tspan..., step = dt)         # "binned" times with width dt
    n_clumps = zeros(Int64, length(ts))     # the number of clumps that are alice at each time
    raft_com = zeros(length(ts), 2)         # the center-of-mass coordinates at each time
    trajs = [reshape(rp.ics.ics[:,i], 2, 1) for i = 1:rp.n_clumps_max] # trajs[i] is a matrix of coordinate locations of clump i
    lifespans = zeros(Bool, length(ts), rp.n_clumps_max) # lifespans[i, j] = true if clump j is alive at time ts[i]
    li = 0 # keeps track of how many steps have been taken such at at least one clump is alive
    
    for (u, t) in TimeChoiceIterator(integrator, ts)
        sum(rp.living) == 0 && continue
        li += 1
        n_clumps[li] = sum(rp.living)
        raft_com[li,:] .= com(u[:,rp.living])

        for i in (1:rp.n_clumps_max)[rp.living]
            trajs[i] = hcat(trajs[i], u[:,i])
            lifespans[li, i] = true
        end

        if showprogress
            t0, tend = integrator.sol.prob.tspan
            val = round(100*(t - t0)/(tend - t0), sigdigits = 3)
            print(WHITE_BG("Integrating: $(val)%   \r"))
            flush(stdout)
        end
    end

    return_raw && return integrator.sol

    l_idx = 1:li
    ts = collect(ts)[l_idx]
    trajs = [Trajectory(permutedims(xy2sph(trajs[i][:,2:end])), ts[lifespans[l_idx,i]]) for i = 1:rp.n_clumps_max if size(trajs[i], 2) > 1]
    trajs = Dict(1:length(trajs) .=> trajs)
    raft_com = raft_com[l_idx,:] |> permutedims |> xy2sph |> permutedims |> r -> Trajectory(r, ts)

    return RaftTrajectory(trajectories = trajs, n_clumps = n_clumps[l_idx], com = raft_com)
end

"""
    rk4(rp; rhs, alg, showprogress, dt, return_raw)

Similar to [`simulate`](@ref) but uses a stock RK4 algorithm with time step `dt` (default `0.1` days).

May be faster for large spring constants at the expense of accuracy.
"""
function rk4(
    rp::RaftParameters; 
    leeway::Bool = false,
    showprogress::Bool = false,
    dt::Real = 0.1,
    return_raw::Bool = false)

    tspan = rp.ics.tspan

    if leeway
        prob_raft = ODEProblem(Leeway!, rp.ics.ics, tspan, rp)
    elseif rp.dx_MR !== nothing && rp.dx_MR !== nothing
        prob_raft = ODEProblem(FastRaft!, rp.ics.ics, tspan, rp)
    else
        prob_raft = ODEProblem(Raft!, rp.ics.ics, tspan, rp)
    end

    integrator = init(prob_raft, SimpleRK4(), dt = dt)

    ts = range(tspan..., step = dt)
    n_clumps = zeros(Int64, length(ts))
    raft_com = zeros(length(ts), 2)
    trajs = [reshape(rp.ics.ics[:,i], 2, 1) for i = 1:rp.n_clumps_max]
    lifespans = zeros(Bool, length(ts), rp.n_clumps_max)
    li = 0
    
    while integrator.t <= last(tspan)
        step!(integrator)
        rp.land(integrator.u, integrator.t, integrator) && rp.land(integrator)
        rp.gd_model(integrator.u, integrator.t, integrator) && rp.gd_model(integrator)
        rp.connections(integrator)

        sum(rp.living) == 0 && continue
        li += 1
        n_clumps[li] = sum(rp.living)
        raft_com[li,:] .= com(integrator.u[:,rp.living])

        for i in (1:rp.n_clumps_max)[rp.living]
            trajs[i] = hcat(trajs[i], integrator.u[:,i])
            lifespans[li, i] = true
        end

        if showprogress
            t0, tend = integrator.sol.prob.tspan
            val = round(100*(t - t0)/(tend - t0), sigdigits = 3)
            print(WHITE_BG("Integrating: $(val)%   \r"))
            flush(stdout)
        end
    end

    return_raw && return integrator.sol

    l_idx = 1:li
    ts = collect(ts)[l_idx]
    trajs = [Trajectory(permutedims(xy2sph(trajs[i][:,2:end])), ts[lifespans[l_idx,i]]) for i = 1:rp.n_clumps_max if size(trajs[i], 2) > 1]
    trajs = Dict(1:length(trajs) .=> trajs)
    raft_com = raft_com[l_idx,:] |> permutedims |> xy2sph |> permutedims |> r -> Trajectory(r, ts)

    return RaftTrajectory(trajectories = trajs, n_clumps = n_clumps[l_idx], com = raft_com)
end