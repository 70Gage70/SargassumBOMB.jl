include(joinpath(@__DIR__, "SargassumBOMB.jl"))

using Optimization 
using OptimizationOptimJL, OptimizationBBO, OptimizationCMAEvolutionStrategy, OptimizationNLopt
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

    return RaftTrajectory(sol, rp, ref_itp, dt = 0.1), tstart, tend
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

u0 = [ClumpParameters(ref_itp).α, ClumpParameters(ref_itp).β]
# u0 = [ClumpParameters(ref_itp).α, 0.02]

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
fig

### USING OPTIMIZATION.JL
loss_opt(u, p) = loss_water(u[1], u[2])
u0 = [ClumpParameters(ref_itp).α, ClumpParameters(ref_itp).β]
lb = [0.0, 0.0]
ub = [0.05, 0.05]

prob = OptimizationProblem(loss_opt, u0, lb = lb, ub = ub)

@info "Optimizing optim."
# @time sol = solve(prob, Optim.ParticleSwarm(lower = prob.lb, upper = prob.ub, n_particles = 100))
@time sol = solve(prob, BBO_adaptive_de_rand_1_bin_radiuslimited())
# @time sol = solve(prob, CMAEvolutionStrategyOpt())
# @time sol = solve(prob, BBO_dxnes())
# @time sol_opt = solve(prob, NLopt.LN_NELDERMEAD())
@show sol_opt

### USING SURROGATES.JL

# loss_opt(u) = loss_water(u[1], u[2])

# n_samples = 100
# lower_bound = [0.0, 0.0]
# upper_bound = [0.05, 0.05]

# @info "Computing surrogate"
# xys = Surrogates.sample(n_samples, lower_bound, upper_bound, GoldenSample())
# @time zs = loss_opt.(xys)

# radial_basis = RadialBasis(xys, zs, lower_bound, upper_bound)

# @info "Optimizing surrogate"
# @time sol_sur = surrogate_optimize(loss_opt, SRBF(), lower_bound, upper_bound, radial_basis, UniformSample(), maxiters=50)
# @show sol_sur

### PLOTTING

# ax = geo_axis(fig[1, 2], title = "Opt", limits = (-90, -38, -5, 22))
# rtr, tstart, tend = integrate_water(sol.u[1], sol.u[2])
# rtr = time_slice(rtr, (tend - 8, tend))
# trajectory_hist!(ax, rtr, x_bins, y_bins)
# land!(ax)
# fig
