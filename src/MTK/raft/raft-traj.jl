include("../models.jl")
include("../control.jl")
include("../../../../CustomMakie.jl/src/geo-methods.jl")
include("../../../../CustomMakie.jl/src/statistic-methods.jl")



######################################################################

@info "Generating model."

x0, y0 = sph2xy(-64, 14, ref_itp) # loop current
# x0, y0 = sph2xy(-60, 25, ref_itp)
t_range = (0.0, 100.0)

x_range = range(start = x0 - 5, length = 5, stop = x0 + 5)
y_range = range(start = y0 - 5, length = 5, stop = y0 + 5)
clump_parameters = ClumpParameters(ref_itp)
spring_parameters = SpringParameters(k -> 20, step(x_range))

@named RRaft = RectangularRaft(x_range, y_range, clump_parameters, spring_parameters)
RRaft = structural_simplify(RRaft)

@info "Generating problem."

prob_raft = ODEProblem(
    RRaft, 
    [],
    t_range, 
    [],
    jac = true,
    sparse = true
)

@info "Solving model."

@time sol = solve(prob_raft)

# sol = solve(prob, 
#     saveat = 0.1, 
#     callback = avoid_shore(water_itp, tol = 0.1, save_positions = (true, false))
# )

@info "Generating reference clump."

@named clump = Clump([x0, y0], clump_parameters)
clump = structural_simplify(clump)

prob_clump = ODEProblem(
    clump, 
    [],
    t_range, 
    []
)

sol_clump = solve(prob_clump)
clump_traj = xy2sph(sol_clump.u, ref_itp)
clump_times = sol_clump.t

@info "Plotting results."

limits = (-100, -50, 5, 35)

### Trajectory COM

fig_COM = default_fig()
ax = geo_axis(fig_COM[1, 1], limits = limits, title = L"\mathrm{Raft COM}")

traj = xy2sph(sol[RRaft.COM], ref_itp)
times = sol.t
trajectory!(ax, traj, times)
land!(ax)

tticks = collect(range(start = minimum(times), stop = maximum(times), length = 5))
data_legend!(fig_COM[1,2], L"\mathrm{Days}", ticks = tticks)

fig_COM

### Trajectory NET

fig_NET = default_fig()
ax = geo_axis(fig_NET[1, 1], limits = limits, title = L"\mathrm{Raft NET}")

for i = 1:Integer(length(states(RRaft))/2)
    clump_x = Symbol("clump_$(i)₊x")
    clump_y = Symbol("clump_$(i)₊y")
    traj = xy2sph([sol[getproperty(RRaft, clump_x)] ;; sol[getproperty(RRaft, clump_y)]], ref_itp)
    times = sol.t
    trajectory!(ax, traj, times)
end

trajectory!(ax, clump_traj, clump_times)

land!(ax)

tticks = collect(range(start = minimum(times), stop = maximum(times), length = 5))
data_legend!(fig_NET[1,2], L"\mathrm{Days}", ticks = tticks)

fig_NET

### Forces

fig_forces = default_fig()
ax = Axis(fig_forces[1, 1])

Fsx = stack(sol[RRaft.clump_1₊τF], dims = 1)[:,1];
τUx = stack(sol[RRaft.clump_1₊τU], dims = 1)[:,1];
Ux = stack(sol[RRaft.clump_1₊U], dims = 1)[:,1];

lines!(ax, sol.t, Fsx)
lines!(ax, sol.t, τUx)
lines!(ax, sol.t, Ux)

fig_forces

# force_U = stack(sol[clump.U], dims = 1)
# force_τU = stack(sol[clump.τU], dims = 1)
# force_ratio = log10.(abs.(force_U ./ force_τU))
# fr_x = force_ratio[:,1]
# fr_y = force_ratio[:,2]

# tmin, tmax = extrema(sol.t)
# ymin, ymax = min(minimum(fr_x), minimum(fr_y)), max(maximum(fr_x), maximum(fr_y))

# ax = axis(fig[1, 1], 
#     title = L"\mathrm{Clump}", 
#     limits = (tmin, tmax, ymin, ymax), 
#     xlabel = L"t \, \text{(days)}",
#     ylabel = L"\log_{10}\left(\frac{F(u)}{F(u_\tau)}\right)",
#     xticks = range(start = tmin, stop = tmax, length = 10),
#     yticks = range(start = ymin, stop = ymax, length = 5))

# lines!(ax, sol.t, force_ratio[:,1], color = :blue)
# lines!(ax, sol.t, force_ratio[:,2], color = :red)

# fig