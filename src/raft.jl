include(joinpath(@__DIR__, "geography.jl"))
include(joinpath(@__DIR__, "physics.jl"))
include(joinpath(@__DIR__, "biology.jl"))
include(joinpath(@__DIR__, "trajectories.jl"))
include(joinpath(@__DIR__, "control.jl"))
include(joinpath(@__DIR__, "../../CustomMakie.jl/src/geo-methods.jl"))
include(joinpath(@__DIR__, "../../CustomMakie.jl/src/statistic-methods.jl"))

######################################################################

@info "Generating model."

# x0, y0 = -64, 14
x0, y0 = -55, 10

x_range, y_range = sph2xy(range(x0 - 1, x0 + 1, step = 0.5), range(y0 - 1, y0 + 1, step = 0.5), ref_itp)
# x_range, y_range = sph2xy(range(x0 - 2, x0 + 2, step = 0.3), range(y0 - 2, y0 + 2, step = 0.3), ref_itp)
tspan = (0.0, 200.0)
# tspan = (100.0, 200)
cp = ClumpParameters(ref_itp)

# spring_k = x -> 20
# rp = RaftParameters(x_range, y_range, cp, spring_k, first(tspan), "nearest")
k10 = 2*step(x_range)
spring_k = x -> 20 * (5/k10) * x * exp(1 - (5/k10)*x) # A(5/k10) * x e^(1 - (5/k10)x)

gd_model = ImmortalModel()
rp = RaftParameters(x_range, y_range, cp, spring_k, first(tspan), "full", gd_model)

prob_raft = ODEProblem(Raft!, rp.ics, tspan, rp)

@info "Solving model."

land = Land(land_itp)
# bm = BrooksModel(BrooksModelParameters(temp_itp, no3_itp))

@time sol_raft = solve(prob_raft, 
    Tsit5(),
    callback = CallbackSet(cb_loc2label(), callback(land), callback(gd_model))
    # callback = CallbackSet(cb_loc2label(), callback(land), grow_test([20.0, 40.0, 60.0, 120.0]))
);

# @time sol_raft = solve(prob_raft, 
#     Tsit5(),
#     callback = die_land(land_itp)
#     # callback = CallbackSet(die_land(land_itp), growth_death_temperature(temp_itp, t_lag = 5.0, T_opt = 26.0, T_thresh = 0.5))
#     # callback = DiscreteCallback(bmp, bmp)
# );

rtr = RaftTrajectory(sol_raft, rp, ref_itp)

@info "Generating reference clump."

xy0 = sph2xy(x0, y0, ref_itp) 
clump_prob = ODEProblem(Clump!, xy0, tspan, cp)

sol_clump = solve(clump_prob, Tsit5())

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