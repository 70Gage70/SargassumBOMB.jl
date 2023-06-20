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

x0, y0 = sph2xy(-64, 14, ref_itp) 
x_range = range(start = x0 - 5, length = 5, stop = x0 + 5)
y_range = range(start = y0 - 5, length = 5, stop = y0 + 5)
clump_parameters = ClumpParameters(ref_itp)
spring_parameters = SpringParameters(k -> 20, step(x_range))

@named RRaft = RectangularRaft(x_range, y_range, clump_parameters, spring_parameters)
RRaft = structural_simplify(RRaft)

# t_range = (0.0, 150.0)
t_range = (0.0, 10.0)

@info "Generating problem."

prob = ODEProblem(
    RRaft, 
    [],
    t_range, 
    [],
    jac = true,
    sparse = true
)

# @time solve(prob)

n_traj = 10
# xy0_ens = [xy0 .+ rand() for i = 1:n_traj]

function prob_func(prob, i, repeat)
    remake(prob, u0 = prob.u0 .+ 2*rand())
end

Eprob = EnsembleProblem(
    prob, 
    prob_func = prob_func
)

@time sim = solve(Eprob, EnsembleThreads(), trajectories = n_traj, saveat = 5.0);


# @info "Solving model."

# sol = solve(prob)

# @variables COM(t)[1:2]
# com = [sol[COM[1]] ;; sol[COM[2]]]
# traj = xy2sph(com, ref_itp)
# lon_traj, lat_traj = (traj[:, 1], traj[:, 2]) 
# times = sol.t

# @info "Plotting results."


