include("../models.jl")
include("../control.jl")
include("../../../../CustomMakie.jl/src/geo-methods.jl")
include("../../../../CustomMakie.jl/src/statistic-methods.jl")

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

@time sol = solve(prob, 
    saveat = 5.0, 
    callback = avoid_shore(water_itp, tol = 0.1, save_positions = (true, false))
)

@info "Plotting results."

### Trajectory

fig_traj = default_fig()
ax = geo_axis(fig_traj[1, 1], limits = (-100, -50, 5, 35), title = L"\mathrm{Clump}")

traj = xy2sph(sol.u, ref_itp)
times = sol.t
trajectory!(ax, traj, times)
land!(ax)

tticks = collect(range(start = minimum(times), stop = maximum(times), length = 5))
data_legend!(fig_traj[1,2], L"\mathrm{Days}", ticks = tticks)

fig_traj

### Forces

fig_force = default_fig()


force_U = stack(sol[clump.U], dims = 1)
force_τU = stack(sol[clump.τU], dims = 1)
force_ratio = log10.(abs.(force_U ./ force_τU))
fr_x = force_ratio[:,1]
fr_y = force_ratio[:,2]

tmin, tmax = extrema(sol.t)
ymin, ymax = min(minimum(fr_x), minimum(fr_y)), max(maximum(fr_x), maximum(fr_y))

ax = axis(fig_force[1, 1], 
    title = L"\mathrm{Clump}", 
    limits = (tmin, tmax, ymin, ymax), 
    xlabel = L"t \, \text{(days)}",
    ylabel = L"\log_{10}\left(\frac{F(u)}{F(u_\tau)}\right)",
    xticks = range(start = tmin, stop = tmax, length = 8),
    yticks = range(start = ymin, stop = ymax, length = 5))

xcolor = colorant"blue"
ycolor = colorant"red"

lines!(ax, sol.t, force_ratio[:,1], color = xcolor)
lines!(ax, sol.t, force_ratio[:,2], color = ycolor)

markers = [
    PolyElement(color = xcolor, strokecolor = :black, strokewidth = 1), 
    PolyElement(color = ycolor, strokecolor = :black, strokewidth = 1)
]

labels = [
    L"x",
    L"y"
]

Legend(fig_force[1,1], markers, labels; (halign = :right, valign = :bottom, tellheight = false, tellwidth = false)...)

fig_force