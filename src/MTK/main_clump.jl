include("models.jl")
include(joinpath(@__DIR__, "..", "plotting.jl"))

###########################################

@info "Loading interpolants."

# ref = EquirectangularReference(lon0 = -75.0, lat0 = 10.0)

# construct_wind_itp(wind_file_default, ref)
# construct_water_itp(water_file_default, ref)

# @load "water_itp.jld"
# @load "wind_itp.jld" 

@info "Generating model."

# xy0 = sph2xy(-64, 14, ref_itp) # loop current, use t_range = (0.0, 200.0)
xy0 = sph2xy(-75, 20, ref_itp) # loop current, use t_range = (0.0, 200.0)
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

@info "Solving model."

condition(u, t, integrator) = (abs(water_itp.u(u..., t)) < 0.1) && (abs(water_itp.v(u..., t)) < 0.1)
affect!(integrator) = terminate!(integrator)
cb = DiscreteCallback(condition, affect!)
sol = solve(prob, saveat = 5.0, callback = cb)

@info "Plotting results."
traj = xy2sph(sol.u, ref_itp)
lon_traj, lat_traj = (traj[:,1], traj[:, 2])
times = sol.t

fig, ax = geo_fig()
ln = trajectory(ax, lon_traj, lat_traj, times)
colorbar(fig, ln)
fig

# CAN GET OBSERVABLES with 
# @variables Fx(t)
# sol[Fx(t)]

# Ensembles

# n_traj = 10000
# xy0_ens = [xy0 .+ rand() for i = 1:n_traj]

# function prob_func(prob, i, repeat)
#     remake(prob, u0 = xy0_ens[i])
# end

# Eprob = EnsembleProblem(
#     prob, 
#     prob_func = prob_func
# )

# @time sim = solve(Eprob, EnsembleThreads(), trajectories = n_traj, saveat = 5.0);