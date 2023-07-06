using ProgressMeter

include("../models.jl")
include("../control.jl")
include("../../../../CustomMakie.jl/src/geo-methods.jl")

######################################################################

t_start = 0.0
t_long = 100.0
t_short = 5.0
n_traj_x = 100
n_traj_y = 100

xmin, xmax = -100.0, -50.0
ymin, ymax = 5.0, 35.0
delta = 0.5

@info "Generating model."

xy0 = [0.0, 0.0] # arbitrary
params = ClumpParameters(ref_itp)
@named clump = Clump(xy0, params)
clump = structural_simplify(clump)

prob = ODEProblem(
    clump, 
    [],
    (t_start, t_long), 
    []
)

@info "Generating ensemble."

n_steps = Integer(t_long/t_short) + 1
x_range = range(start = xmin + delta, length = n_traj_x, stop = xmax - delta)
y_range = range(start = ymin + delta, length = n_traj_y, stop = ymax - delta)
x_range, y_range = sph2xy(x_range, y_range, ref_itp)
xy0 = [
    [x, y] + [step(x_range)*(rand() - 1/2), step(y_range)*(rand() - 1/2)]  # random wiggle in x0, y0 to probe space more
    for x in x_range for y in y_range
] 

n_traj = length(xy0)

# prob_fun tells the ensemble how to "start over", in this case the i'th run is the problem remade with the i'th initial condition and time range.
function prob_func(prob, i, repeat)
    remake(prob, u0 = xy0[i])
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

sim = solve(Eprob, EnsembleThreads(), trajectories = n_traj, 
    saveat = t_short, 
    callback = avoid_shore(water_itp, tol = 0.1, save_positions = (true, false))
)

finish!(progress)

@info "Plotting results."

fig = default_fig()
ax = geo_axis(fig[1, 1], limits = (-100, -50, 5, 35), title = L"\mathrm{Clump}")


for sol in sim
    traj = xy2sph(sol.u, ref_itp)
    times = sol.t
    trajectory!(ax, traj, times)
end

land!(ax)

tticks = collect(range(start = minimum(times), stop = maximum(times), length = 5))
data_legend!(fig[1,2], L"\mathrm{Days}", ticks = tticks)

fig