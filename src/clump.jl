include(joinpath(@__DIR__, "models.jl"))
include(joinpath(@__DIR__, "control.jl"))
include(joinpath(@__DIR__, "../../CustomMakie.jl/src/geo-methods.jl"))
include(joinpath(@__DIR__, "../../CustomMakie.jl/src/statistic-methods.jl"))

######################################################################

@info "Generating model."

# xy0 = sph2xy(-64, 14, ref_itp) 
xy0 = sph2xy(-80, 14, ref_itp) 
tspan = (0.0,1000.0)
cp = ClumpParameters(ref_itp)
clump_prob = ODEProblem(Clump!, xy0, tspan, cp)

@info "Solving model."

@time sol = solve(clump_prob, 
    saveat = 5.0,
    callback = die_land(land_itp)
)

@info "Plotting results."

### Trajectory

fig_traj2 = default_fig()
ax = geo_axis(fig_traj2[1, 1], limits = (-100, -50, 5, 35), title = L"\mathrm{Clump}")

traj = xy2sph(sol.u, ref_itp)
times = sol.t
trajectory!(ax, traj, times)
land!(ax)

tticks = collect(range(start = minimum(times), stop = maximum(times), length = 5))
data_legend!(fig_traj2[1,2], L"\mathrm{Days}", ticks = tticks)

fig_traj2