include(joinpath(@__DIR__, "SargassumBOMB.jl"))

using Optimization 
using OptimizationOptimJL, OptimizationBBO, OptimizationCMAEvolutionStrategy, OptimizationNLopt

########################################################################

function integrate_water(
    α::Real, 
    τ::Real;
    initial_time::NTuple{2, Integer} = (2018, 3), 
    final_time::NTuple{2, Integer} = (2018, 4),
    t_extra::Real = 7)

    tstart = Day(DateTime(initial_time...) - DateTime(yearmonth(water_itp.time_start)...)).value |> float
    tend = tstart + Day(DateTime(final_time...) - DateTime(initial_time...)).value + t_extra
    tspan = (tstart, tend)

    dists = SargassumDistribution(joinpath(@__DIR__, "..", "..", "SargassumFromAFAI.jl", "data", "dist-2018.nc"))
    dist = dists[initial_time]
    ics = initial_conditions(dist, [1], 1, "uniform", ref_itp)

    cp_default = ClumpParameters(ref_itp) 
    cp = ClumpParameters(ref_itp, α, τ, cp_default.R, cp_default.f, cp_default.β)

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

    return RaftTrajectory(sol, rp, ref_itp, dt = 0.1), tstart, tend
end

function loss_water(
    α::Real, 
    τ::Real;
    initial_time::NTuple{2, Integer} = (2018, 3), 
    final_time::NTuple{2, Integer} = (2018, 4))

    target = dists[final_time].sargassum[:,:,1]
    target = target/sum(target)

    rtr, tstart, tend = integrate_water(α, τ, initial_time = initial_time, final_time = final_time)
    rtr = time_slice(rtr, (tend - 8, tend))
    data = bins(rtr, dists[final_time])
    data = data/sum(data)

    return sum((data - target) .^ 2)
end

#################################################################

u0 = [ClumpParameters(ref_itp).α, ClumpParameters(ref_itp).τ]

fig = Figure(
    resolution = (1920, 800), 
    fontsize = 50,
    figure_padding = (5, 100, 5, 5));

dist = dists[(2018, 4)]
lons = dist.lon
δx = (lons[2] - lons[1])/2
x_bins = range(lons[1] - δx, stop = dist.lon[end] + δx, length = length(lons) + 1)

lats = dist.lat
δy = (lats[2] - lats[1])/2
y_bins = range(lats[1] - δy, stop = dist.lat[end] + δy, length = length(lats) + 1)  

ax = geo_axis(fig[1, 1], title = "Default", limits = (-90, -38, -5, 22))
rtr, tstart, tend = integrate_water(u0[1], u0[2])
rtr = time_slice(rtr, (tend - 8, tend))
# rtr = time_slice(rtr, (tstart, tstart + 0.1))
trajectory_hist!(ax, rtr, x_bins, y_bins)
land!(ax)
# fig

# objective function (to be minimized)
loss_opt(u, p) = -loss_water(u[1], u[2])
# u0 = [ClumpParameters(ref_itp).α, ClumpParameters(ref_itp).τ]
lb = [0.0, 0.0]
ub = [0.05, 0.1]

prob = OptimizationProblem(loss_opt, u0, lb = lb, ub = ub)

@info "Optimizing."
# @time sol = solve(prob, Optim.ParticleSwarm(lower = prob.lb, upper = prob.ub, n_particles = 100))
# @time sol = solve(prob, BBO_adaptive_de_rand_1_bin_radiuslimited())
# @time sol = solve(prob, CMAEvolutionStrategyOpt())
# @time sol = solve(prob, BBO_dxnes())
# @time sol = solve(prob, NLopt.LN_NELDERMEAD())

# ax = geo_axis(fig[1, 2], title = "Opt", limits = (-90, -38, -5, 22))
# rtr, tstart, tend = integrate_water(sol.u[1], sol.u[2])
# rtr = time_slice(rtr, (tend - 8, tend))
# trajectory_hist!(ax, rtr, x_bins, y_bins)
# land!(ax)
# fig
