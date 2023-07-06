include("../models.jl")
include("../control.jl")
include("../../../../CustomMakie.jl/src/geo-methods.jl")

######################################################################

@info "Generating model."

xy0 = sph2xy(-64, 14, ref_itp) # loop current, use t_range = (0.0, 200.0)
# xy0 = sph2xy(-75, 20, ref_itp) # loop current, use t_range = (0.0, 200.0)
params = ClumpParameters(ref_itp)
@named clump = Clump(xy0, params)
clump = structural_simplify(clump)

t_range = (0.0, 200.0)

prob = ODEProblem(
    clump, 
    [],
    t_range, 
    []
)

@info "Solving model."

sol = solve(prob, 
    saveat = 5.0, 
    callback = avoid_shore(water_itp, tol = 0.1, save_positions = (true, false))
)

@info "Plotting results."

fig = default_fig()
ax = geo_axis(fig[1, 1], limits = (-100, -50, 5, 35), title = L"\mathrm{Clump}")

traj = xy2sph(sol.u, ref_itp)
times = sol.t
trajectory!(ax, traj, times)
land!(ax)

tticks = collect(range(start = minimum(times), stop = maximum(times), length = 5))
data_legend!(fig[1,2], L"\mathrm{Days}", ticks = tticks)

fig