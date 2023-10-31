include(joinpath(@__DIR__, "SargassumBOMB.jl"))

# using Optimization 
# using OptimizationOptimJL, OptimizationBBO, OptimizationCMAEvolutionStrategy, OptimizationNLopt
using Surrogates

########################################################################

function integrate_water(
    α::Real, 
    β::Real;
    initial_time::NTuple{2, Integer} = (2018, 3), 
    final_time::NTuple{2, Integer} = (2018, 4),
    t_extra::Real = 7,
    seed::Integer = 1234)

    seed!(seed)

    tstart = Day(DateTime(initial_time...) - DateTime(yearmonth(water_itp.time_start)...)).value |> float
    tend = tstart + Day(DateTime(final_time...) - DateTime(initial_time...)).value + t_extra
    tspan = (tstart, tend)

    dists = SargassumDistribution(joinpath(@__DIR__, "..", "..", "SargassumFromAFAI.jl", "data", "dist-2018.nc"))
    dist = dists[initial_time]
    ics = initial_conditions(dist, [1], 2, "levels", ref_itp)

    cp_default = ClumpParameters(ref_itp) 
    cp = ClumpParameters(ref_itp, α, cp_default.τ, cp_default.R, cp_default.f, β)

    sp = SpringParameters(k -> 0.0, 0.0)

    nw_type = "none"
    icons = form_connections(ics, nw_type)

    gdm = ImmortalModel()

    land = Land(verbose = false)

    rp = RaftParameters(
        ics = ics,
        clumps = cp,
        springs = sp,
        connections = icons,
        t0 = first(tspan),
        gd_model = gdm
    )

    prob = ODEProblem(WaterWind!, rp.ics, tspan, rp)

    sol = solve(
        prob, 
        Tsit5(), abstol = 1e-6, reltol = 1e-6,
        callback = CallbackSet(
            cb_update(showprogress = false), 
            callback(land), 
            callback(gdm), 
            cb_connections(network_type = nw_type))
        )

    # return (sol, tstart, tend)
    return (RaftTrajectory(sol, rp, ref_itp, dt = 0.1), tstart, tend)
end

function loss_water(
    α::Real, 
    β::Real;
    initial_time::NTuple{2, Integer} = (2018, 3), 
    final_time::NTuple{2, Integer} = (2018, 4))

    target = dists[final_time].sargassum[:,:,1]
    target = target/sum(target)

    rtr, tstart, tend = integrate_water(α, β, initial_time = initial_time, final_time = final_time)
    rtr = time_slice(rtr, (tend - 8, tend))
    data = bins(rtr, dists[final_time])
    data = data/sum(data)

    return sum((data - target) .^ 2)
end

#################################################################
# OPTIMIZING

### USING OPTIMIZATION.JL
# loss_opt(u, p) = loss_water(u[1], u[2])
# u0 = [ClumpParameters(ref_itp).α, ClumpParameters(ref_itp).β]
# lb = [0.0, 0.0]
# ub = [0.05, 0.05]

# prob = OptimizationProblem(loss_opt, u0, lb = lb, ub = ub)

# @info "Optimizing optim."
# @time sol = solve(prob, Optim.ParticleSwarm(lower = prob.lb, upper = prob.ub, n_particles = 100))
# @time sol = solve(prob, BBO_adaptive_de_rand_1_bin_radiuslimited())
# @time sol = solve(prob, CMAEvolutionStrategyOpt())
# @time sol = solve(prob, BBO_dxnes())
# @time sol_opt = solve(prob, NLopt.LN_NELDERMEAD())
# @show sol_opt

### USING SURROGATES.JL
loss_opt(u) = loss_water(u[1], u[2])

n_samples_sur = 100         # default 100
maxiters_opt = 50           # default 50
lower_bound = [0.0, 0.0]    # [α, β] lower
upper_bound = [0.05, 0.05]  # [α, β] upper

@info "Computing surrogate"
xys = Surrogates.sample(n_samples_sur, lower_bound, upper_bound, SobolSample())
@time zs = loss_opt.(xys)

radial_basis = RadialBasis(xys, zs, lower_bound, upper_bound)

@info "Optimizing surrogate"
@time sol_sur = surrogate_optimize(loss_opt, DYCORS(), lower_bound, upper_bound, 
                                    radial_basis, SobolSample(), maxiters = maxiters_opt)


default_loss = loss_opt([ClumpParameters(ref_itp).α, ClumpParameters(ref_itp).β])
optimized_loss = loss_opt([sol_sur[1][1], sol_sur[1][2]])

@info "Default loss: $(default_loss)"
@info "Optmzed loss: $(optimized_loss)"

#################################################################
### PLOTTING

fig = Figure(
    # resolution = (1920, 1080), 
    resolution = (2220, 1920),
    fontsize = 50,
    figure_padding = (5, 100, 5, 5))

limits = (-100, -40, 5, 35)

### AFAI
# initial distribution (AFAI)
ax = geo_axis(fig[1, 1], limits = limits, title = "AFAI initial (March week 1)")
SFA_plot!(ax, (2018, 3), 1)
land!(ax)

# final distribution (AFAI)
ax = geo_axis(fig[1, 2], limits = limits, title = "AFAI final (April week 1)")
SFA_plot!(ax, (2018, 4), 1)
land!(ax)

### UNOPTIMIZED
# initial distribution (RAFT, unoptimized)
ax = geo_axis(fig[2, 1], limits = limits, title = "RAFT initial [default] (March week 1)")
rtr_dt, tstart, tend = integrate_water(ClumpParameters(ref_itp).α, ClumpParameters(ref_itp).β)
dist = dists[(2018, 3)]
rtr_dt_initial = time_slice(rtr_dt, (first(rtr_dt.t), first(rtr_dt.t)))
trajectory_hist!(ax, rtr_dt_initial, dist)
land!(ax)

# final distribution (RAFT, unoptimized)
ax = geo_axis(fig[2, 2], limits = limits, title = "RAFT final [default] (April week 1)")
rtr_final = time_slice(rtr_dt, (tend - 8, tend))
trajectory_hist!(ax, rtr_final, dist)
land!(ax)

### OPTIMIZED
# initial distribution (RAFT, unoptimized)
ax = geo_axis(fig[3, 1], limits = limits, title = "RAFT initial [optim] (March week 1)")
rtr_dt, tstart, tend = integrate_water(sol_sur[1][1], sol_sur[1][2])
dist = dists[(2018, 3)]
rtr_dt_initial = time_slice(rtr_dt, (first(rtr_dt.t), first(rtr_dt.t)))
trajectory_hist!(ax, rtr_dt_initial, dist)
land!(ax)

# final distribution (RAFT, unoptimized)
ax = geo_axis(fig[3, 2], limits = limits, title = "RAFT final [optim] (April week 1)")
rtr_final = time_slice(rtr_dt, (tend - 8, tend))
trajectory_hist!(ax, rtr_final, dist)
land!(ax)

fig[0,:] = Label(fig, "Default: $(default_loss), Optimized: $(optimized_loss)")

fig
