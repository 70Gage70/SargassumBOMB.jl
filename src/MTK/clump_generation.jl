include("models.jl")
include(joinpath(@__DIR__, "..", "plotting.jl"))

###########################################

@info "Generating model."

xy0 = [0.0, 0.0] # arbitrary
params = ClumpParameters(ref_itp)
@named clump = Clump(xy0, params)
clump = structural_simplify(clump)

t_range = (0.0, 100.0)

prob = ODEProblem(
    clump, 
    [],
    t_range, 
    []
)

# callbacks, so computations end when a particle exits the water
condition(u, t, integrator) = (abs(water_itp.u(u..., t)) < 0.1) && (abs(water_itp.v(u..., t)) < 0.1)
affect!(integrator) = terminate!(integrator)
cb = DiscreteCallback(condition, affect!)

# ensembles
n_traj_x = 100
n_traj_y = 100
x_range = range(start = -100.0 + 0.5, length = n_traj_x, stop = -50.0 - 0.5)
y_range = range(start = 5.0 + 0.5, length = n_traj_y, stop = 35.0 - 0.5)
x_range, y_range = sph2xy(x_range, y_range, ref_itp)
xy0 = [[x, y] for x in x_range for y in y_range]
n_traj = length(xy0)

function prob_func(prob, i, repeat)
    remake(prob, u0 = xy0[i])
end

Eprob = EnsembleProblem(
    prob, 
    prob_func = prob_func
)

@info "Solving model."

@time sim = solve(Eprob, EnsembleThreads(), trajectories = n_traj, saveat = 5.0, callback = cb);


# @info "Plotting results."
# traj = xy2sph(sol.u, ref_itp)
# lon_traj, lat_traj = (traj[:,1], traj[:, 2])
# times = sol.t

# fig, ax = geo_fig()
# ln = trajectory(ax, lon_traj, lat_traj, times)
# colorbar(fig, ln)
# fig


