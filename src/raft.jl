include(joinpath(@__DIR__, "geography.jl"))
include(joinpath(@__DIR__, "physics.jl"))
include(joinpath(@__DIR__, "biology.jl"))
include(joinpath(@__DIR__, "trajectories.jl"))
include(joinpath(@__DIR__, "control.jl"))
include(joinpath(@__DIR__, "../../CustomMakie.jl/src/geo-methods.jl"))
include(joinpath(@__DIR__, "../../CustomMakie.jl/src/statistic-methods.jl"))

######################################################################

@info "Generating model."

# x0, y0 = -90, 23 # GoM
# x0, y0 = -64, 14
x0, y0 = -55, 10

# x_range, y_range = sph2xy(range(x0 - 1, x0 + 1, step = 0.5), range(y0 - 1, y0 + 1, step = 0.5), ref_itp)
x_range, y_range = sph2xy(range(x0 - 2, x0 + 2, step = 0.3), range(y0 - 2, y0 + 2, step = 0.3), ref_itp)
# tspan = (0.0, 200.0)
tspan = (151.0, 300) # May 1 - Sep. 27
cp = ClumpParameters(ref_itp)

# spring_k = x -> 20
# rp = RaftParameters(x_range, y_range, cp, spring_k, first(tspan), "nearest")
# k10 = 2*step(x_range)
# spring_k = x -> 5 * (5/k10) * x * exp(1 - (5/k10)*x) # A(5/k10) * x e^(1 - (5/k10)x)

function spring_k(x::Real; k10::Real = 2*step(x_range))
    return 5 * (5/k10) * x * exp(1 - (5/k10)*x)
end

# gd_model = ImmortalModel()
gd_model = BrooksModel(verbose = true)
rp = RaftParameters(x_range, y_range, cp, spring_k, first(tspan), "full", gd_model)

prob_raft = ODEProblem(Raft!, rp.ics, tspan, rp)

@info "Solving model."

land = Land(verbose = true)

@time sol_raft = solve(prob_raft, 
    Tsit5(), abstol = 1e-6, reltol = 1e-6,
    callback = CallbackSet(cb_loc2label(), callback(land), callback(gd_model))
);

rtr = RaftTrajectory(sol_raft, rp, ref_itp)

@info "Generating reference clump."

gd_model = ImmortalModel()
land = Land(verbose = true)
rp_1c = OneClumpRaftParameters(sph2xy(x0, y0, ref_itp)..., cp, first(tspan), gd_model)
prob_raft_1c = ODEProblem(Raft!, rp_1c.ics, tspan, rp_1c)

@time sol_clump = solve(prob_raft_1c, 
    Tsit5(), abstol = 1e-6, reltol = 1e-6,
    callback = CallbackSet(cb_loc2label(), callback(land), callback(gd_model))
);

ctr = RaftTrajectory(sol_clump, rp_1c, ref_itp)

@info "Plotting results."

# limits = (-100, -50, 5, 35)
limits = (-70, -50, 5, 25)

### Trajectory COM

fig_COM = default_fig()
ax = geo_axis(fig_COM[1, 1], limits = limits, title = L"\mathrm{Raft COM}")

# raft
# trajectory!(ax, rtr)

# COM
trajectory!(ax, rtr.com, 
    opts = (linestyle = :dot, color = rtr.com.t, colormap = :heat, linewidth = 5)
) 

# COM 5
rtr_5 = RaftTrajectory(sol_raft, rp, ref_itp, dt = 5.0)

# trajectory!(ax, rtr_5.com, 
#     opts = (linestyle = :dot, color = rtr_5.com.t, colormap = :viridis, linewidth = 5)
# ) 

scatter!(ax, rtr_5.com.xy[:,1], rtr_5.com.xy[:,2])

# clump
# trajectory!(ax, ctr, 
#     opts = (color = ctr.t, colormap = :heat, linewidth = 2)
# ) 

land!(ax)

tticks = collect(range(start = minimum(rtr.t), stop = maximum(rtr.t), length = 5))
data_legend!(fig_COM[1,2], L"\mathrm{Days}", ticks = tticks)

fig_COM