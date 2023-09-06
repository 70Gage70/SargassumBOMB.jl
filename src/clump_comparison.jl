include(joinpath(@__DIR__, "geography.jl"))
include(joinpath(@__DIR__, "physics.jl"))
include(joinpath(@__DIR__, "biology.jl"))
include(joinpath(@__DIR__, "trajectories.jl"))
include(joinpath(@__DIR__, "control.jl"))
include(joinpath(@__DIR__, "../../CustomMakie.jl/src/geo-methods.jl"))
include(joinpath(@__DIR__, "../../CustomMakie.jl/src/statistic-methods.jl"))

######################################################################

# Clump!

# x0, y0 = -55, 10
x0, y0 = -64, 14
tspan = (151.0, 300) # May 1 - Sep. 27
cp = ClumpParameters(ref_itp)

xy0 = sph2xy(x0, y0, ref_itp) 
clump_prob = ODEProblem(Clump!, xy0, tspan, cp)

# @btime sol_clump = solve(clump_prob, Tsit5())
sol_clump = solve(clump_prob, Tsit5(), abstol = 1e-6, reltol = 1e-6)
ctr = Trajectory(sol_clump.u, sol_clump.t, ref_itp)

# xy02 = sph2xy(x0 + 0.01*rand(), y0 + 0.01*rand(), ref_itp)
xy02 = [1205.378748125225, 445.358753207288] 
clump_prob2 = ODEProblem(Clump!, xy02, tspan, cp)
sol_clump2 = solve(clump_prob2, Tsit5(), abstol = 1e-6, reltol = 1e-6)
ctr2 = Trajectory(sol_clump2.u, sol_clump2.t, ref_itp)


# Raft!
gd_model = ImmortalModel()
rp = RPClumpTest(xy0[1], xy0[2], cp, first(tspan), gd_model)
prob_raft = ODEProblem(Raft!, rp.ics, tspan, rp)
land = Land(verbose = true)

# @btime sol_raft = solve(prob_raft, 
#     Tsit5(),
#     callback = CallbackSet(cb_loc2label(), callback(land), callback(gd_model))
# );

sol_raft = solve(prob_raft, 
    Tsit5(), abstol = 1e-6, reltol = 1e-6,
    callback = CallbackSet(cb_loc2label(), callback(land), callback(gd_model))
);

rtr = RaftTrajectory(sol_raft, rp, ref_itp)


####################

fig_COM = default_fig()
limits = (-100, -50, 5, 35)
ax = geo_axis(fig_COM[1, 1], limits = limits, title = L"\mathrm{Raft COM}")

# raft
trajectory!(ax, rtr)

# clump
trajectory!(ax, ctr, 
    opts = (color = ctr.t, colormap = :heat, linewidth = 2)
) 

trajectory!(ax, ctr2, 
    opts = (color = ctr2.t, colormap = :heat, linewidth = 2)
) 

land!(ax)

tticks = collect(range(start = minimum(rtr.t), stop = maximum(rtr.t), length = 5))
data_legend!(fig_COM[1,2], L"\mathrm{Days}", ticks = tticks)

fig_COM
