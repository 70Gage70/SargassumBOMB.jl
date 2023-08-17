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
    callback = die_land(land_itp)
)

tr = Trajectory(sol.u, sol.t, ref_itp)

@info "Plotting results."

### Trajectory

limits = (-100, -50, 5, 35)

fig_clump_traj = default_fig()
ax = geo_axis(fig_clump_traj[1, 1], limits = limits, title = L"\mathrm{Clump}")

trajectory!(ax, tr)
land!(ax)

tticks = collect(range(start = minimum(tr.t), stop = maximum(tr.t), length = 5))
data_legend!(fig_clump_traj[1,2], L"\mathrm{Days}", ticks = tticks)

fig_clump_traj