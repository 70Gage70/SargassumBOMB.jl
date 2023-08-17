include(joinpath(@__DIR__, "models.jl"))
include(joinpath(@__DIR__, "control.jl"))
include(joinpath(@__DIR__, "../../CustomMakie.jl/src/geo-methods.jl"))
include(joinpath(@__DIR__, "../../CustomMakie.jl/src/statistic-methods.jl"))

######################################################################

@info "Generating model."

# x0, y0 = -64, 14
x0, y0 = -55, 10

x_range, y_range = sph2xy(range(x0 - 1, x0 + 1, step = 0.5), range(y0 - 1, y0 + 1, step = 0.5), ref_itp)
tspan = (0.0, 200.0)
cp = ClumpParameters(ref_itp)
spring_k = x -> 20
rp = RaftParameters(x_range, y_range, cp, spring_k)
prob_raft = ODEProblem(Raft!, rp.xy0, tspan, rp)

@info "Solving model."

# @time sol_no_shore = solve(prob_raft, Tsit5(), reltol = 1e-6, abstol = 1e-6)

@time sol_raft = solve(prob_raft, 
    Tsit5(),
    #callback = die_land(water_itp, tol = 0.1)
    callback = CallbackSet(die_land(land_itp), grow_test([5.0, 50.0]))
)

rtr = RaftTrajectory(sol_raft, rp, ref_itp)

@info "Generating reference clump."

xy0 = sph2xy(x0, y0, ref_itp) 
clump_prob = ODEProblem(Clump!, xy0, tspan, cp)

sol_clump = solve(clump_prob)

ctr = Trajectory(sol_clump.u, sol_clump.t, ref_itp)

@info "Plotting results."

limits = (-100, -50, 5, 35)

### Trajectory COM

fig_COM = default_fig()
ax = geo_axis(fig_COM[1, 1], limits = limits, title = L"\mathrm{Raft COM}")

# raft
trajectory!(ax, rtr)

# COM
trajectory!(ax, rtr.com, 
    opts = (linestyle = :dot, color = rtr.com.t, linewidth = 5)
) 

# clump
trajectory!(ax, ctr, 
    opts = (color = ctr.t, colormap = :heat, linewidth = 2)
) 

land!(ax)

tticks = collect(range(start = minimum(rtr.t), stop = maximum(rtr.t), length = 5))
data_legend!(fig_COM[1,2], L"\mathrm{Days}", ticks = tticks)

fig_COM

# ### Trajectory NET

# fig_NET = default_fig()
# ax = geo_axis(fig_NET[1, 1], limits = limits, title = L"\mathrm{Raft NET}")

# for i = 1:Integer(length(rp.xy0)/2)
#     times = sol.t
#     traj = xy2sph(Rafti(sol, i), ref_itp)
#     trajectory!(ax, traj, times)
# end



# land!(ax)

# tticks = collect(range(start = minimum(times), stop = maximum(times), length = 5))
# data_legend!(fig_NET[1,2], L"\mathrm{Days}", ticks = tticks)

# fig_NET