using ProgressMeter
using MAT

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
# save_positions ensures we don't save before or after the condition is satisfied, avoiding partial trajectories
cb = DiscreteCallback(condition, affect!, save_positions = (false, false))  

# generating initial conditions for each run
t_long = 100.0
t_short = 5.0
n_steps = Integer(t_long/t_short) + 1
n_traj_x = 100
n_traj_y = 100
x_range = range(start = -100.0 + 0.5, length = n_traj_x, stop = -50.0 - 0.5)
y_range = range(start = 5.0 + 0.5, length = n_traj_y, stop = 35.0 - 0.5)
x_range, y_range = sph2xy(x_range, y_range, ref_itp)
tspans = [(i, i + t_long) for i = 0.0:1.0:260.0]
xyt0 = [
    [[x, y] + [step(x_range)*(rand() - 1/2), step(y_range)*(rand() - 1/2)], t]  # random wiggle in x0, y0 to probe space more
    for x in x_range for y in y_range for t in tspans
] 

n_traj = length(xyt0)

# prob_fun tells the ensemble how to "start over", in this case the i'th run is the problem remade with the i'th initial condition and time range.
function prob_func(prob, i, repeat)
    remake(prob, u0 = xyt0[i][1], tspan = xyt0[i][2])
end

# output_func generally handles what is done with the solution at each iteration; here we only use it to keep track of progress
progress = Progress(n_traj)

function output_func(sol, i)
    next!(progress)
    sol, false
end

Eprob = EnsembleProblem(
    prob, 
    prob_func = prob_func,
    output_func = output_func
)

@info "Solving model."

@time sim = solve(Eprob, EnsembleThreads(), trajectories = n_traj, saveat = t_short, callback = cb);

finish!(progress)

@info "Plotting results."

fig, ax = geo_fig(title = L"\mathrm{BOM} \, 100 \, \mathrm{day} \,  \mathrm{trajectories}")

for sol in sim.u
    if length(sol.t) == n_steps
        ln = trajectory(ax, sol, ref_itp)
    end
end

# colorbar(fig, ln)
fig

@info "Writing results."

x0_long = Float64[]
y0_long = Float64[]
xT_long = Float64[]
yT_long = Float64[]

x0_short = Float64[]
y0_short = Float64[]
xT_short = Float64[]
yT_short = Float64[]

for sol in sim.u
    if length(sol.t) == n_steps
        xy_sol = xy2sph(sol.u, ref_itp)

        push!(x0_long, xy_sol[1, 1])
        push!(y0_long, xy_sol[1, 2])
        push!(xT_long, xy_sol[end, 1])
        push!(yT_long, xy_sol[end, 2])
        
        for j = 1:n_steps - 1
            push!(x0_short, xy_sol[j, 1])
            push!(y0_short, xy_sol[j, 2])
            push!(xT_short, xy_sol[j + 1, 1])
            push!(yT_short, xy_sol[j + 1, 2])
        end
    end
end            

matwrite("x0xT_BOM_clump_5days.mat", Dict("x0" => x0_short, "y0" => y0_short, "xT" => xT_short, "yT" => yT_short))
matwrite("x0xT_BOM_clump_100days.mat", Dict("x0" => x0_long, "y0" => y0_long, "xT" => xT_long, "yT" => yT_long))
        
