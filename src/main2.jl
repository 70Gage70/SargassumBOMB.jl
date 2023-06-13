using GeoMakie, CairoMakie

include("models.jl")
include("vector-field-files.jl")

###########################################

@info "Loading interpolants."

ref = EquirectangularReference(lon0 = -75.0, lat0 = 10.0)

# construct_wind_itp(wind_file_default, ref)
# construct_water_itp(water_file_default, ref)

@load "water_itp.jld"
@load "wind_itp.jld" 

@info "Generating model."

# ics = sph2xy(-79.5, 25.5, ref) # gulf stream, use t_range = (0.0, 10.0)
xy0 = sph2xy(-64, 14, ref) # loop current, use t_range = (0.0, 200.0)
params = BOMParameters(ref)
@named clump = Clump(xy0, params)
clump = structural_simplify(clump)

t_range = (0.0, 100.0)

prob = ODEProblem(
    clump, 
    [],
    t_range, 
    []
)

function prob_func(prob, i, repeat)
    remake(prob, u0 = xy0 .+ rand())
end

Eprob = EnsembleProblem(
    prob, 
    prob_func = prob_func
)

sim(n_traj) = solve(Eprob, Tsit5(), EnsembleThreads(), trajectories = n_traj)

# @info "Solving model."

sol = solve(prob)
traj = xy2sph(sol.u, ref)
lon_traj, lat_traj = (traj[:,1], traj[:, 2]) #
times = sol.t

# @info "Plotting results."

# fig = Figure(resolution = (1920, 1080))

# ga(fig, row, col, title) = GeoAxis(
#     fig[row, col],
#     dest = "+proj=eqc", # https://proj.org/en/9.2/operations/projections/eqc.html
#     lonlims = (-100, -50),
#     latlims = (5, 35),
#     coastlines = true,
#     title = title
# )


# ln = lines!(ga(fig, 1, 1, "traj"), lon_traj, lat_traj; color = times, linewidth = 4)
# Colorbar(fig[1,2], ln, label = "Days")

# CAN GET OBSERVABLES with 
# @variables Fx(t)
# sol[Fx(t)]